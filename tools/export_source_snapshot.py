#!/usr/bin/env python3
"""Export retained research source without altering the original workspace.

Standard library only. This tool does not execute research code, download data,
change experiment parameters, create GitHub repositories, or push commits.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
from collections import Counter


SOURCE_EXTENSIONS = {
    ".py", ".pyw", ".ipynb", ".sh", ".ps1", ".bat", ".cmd", ".r", ".m",
    ".jl", ".js", ".mjs", ".ts", ".tsx", ".jsx", ".c", ".cpp", ".h",
    ".hpp", ".yaml", ".yml", ".toml",
}
TEXT_EXTENSIONS = SOURCE_EXTENSIONS | {".md", ".json", ".csv", ".tsv", ".txt"}
NO_WALK = {
    ".git", ".venv-dpas", ".venv", "venv", "__MACOSX", "__pycache__",
    "node_modules", "data", "artifacts", "tmp", "git-publication",
    "archives", "download_attempts", "pyg_wheels",
}
ROOT_FILES = {
    "README.md", ".gitignore", ".gitattributes", "PROJECT_INDEX.md",
    "PROJECT_PLAN.md", "FILE_PROVENANCE.md", "ST2P_TASK_SCOPE.md",
    "retained_cleanup_manifest.json", "AGENTS.md",
}
TREES = ("ST2P", "SCP-NET-work", "results", "figures", "tools")
SECRET_PATTERNS = {
    "github_token": re.compile(rb"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b"),
    "private_key": re.compile(rb"-----BEGIN (?:[A-Z0-9 ]+ )?PRIVATE KEY-----"),
    "aws_access_key": re.compile(rb"\bAKIA[0-9A-Z]{16}\b"),
    "openai_style_key": re.compile(rb"\bsk-(?:proj-)?[A-Za-z0-9_-]{32,}\b"),
    # A shell reference such as $HF_TOKEN is not a literal credential.
    "credential_url": re.compile(rb"https?://[^\s/:<>]+:(?!\$)[^\s/@<>]+@"),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def retained_files(root: Path):
    for dirname, dirs, names in os.walk(root, followlinks=False):
        dirs[:] = [d for d in dirs if d not in NO_WALK and not d.startswith(".venv")]
        for name in sorted(names):
            p = Path(dirname) / name
            if p.is_symlink() or name.startswith(("thread_", "GIT_UPLOAD_RECEIPT_")):
                continue
            if p.suffix.lower() in TEXT_EXTENSIONS or name in {".gitignore", ".gitattributes", "Makefile", "Dockerfile", "LICENSE", "COPYING"}:
                if not name.endswith(".inspect.ndjson"):
                    yield p


def upstream_records(root: Path):
    primary = json.loads((root / "upstream/SOURCE_MANIFEST.json").read_text(encoding="utf-8-sig"))
    projects = {item["name"]: item for item in primary["repositories"]}
    for name, filename in [("NicheTrans", "NicheTrans_SOURCE_MANIFEST.json"), ("SPINE", "SPINE_SOURCE_MANIFEST.json")]:
        record = json.loads((root / "upstream" / filename).read_text(encoding="utf-8-sig"))
        projects[name] = {
            "name": name, "commit": record["commit"],
            "official_repository": record["repository"], "files": record["files"],
        }
    return projects


def export(source: Path, output: Path):
    source = source.resolve(strict=True)
    output = output.resolve()
    if output == source or source in output.parents and output.name != "git-publication":
        raise ValueError("Use a separate export directory, or the dedicated git-publication directory")
    if output in source.parents:
        raise ValueError("The output must not contain the original workspace")
    output.mkdir(parents=True, exist_ok=True)
    files = []
    credentials = []
    checked_python = 0
    checked_notebooks = 0
    syntax_failures = []

    def write_file(rel: str, data: bytes, provenance: str, mode: str = "100644"):
        nonlocal checked_python, checked_notebooks
        rel_path = Path(rel)
        if rel_path.is_absolute() or ".." in rel_path.parts or ".git" in rel_path.parts:
            raise ValueError(f"Unsafe archive path: {rel}")
        dest = output / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        if dest.read_bytes() != data:
            raise RuntimeError(f"Copy verification failed: {rel}")
        if len(data) >= 100 * 1024 * 1024:
            raise RuntimeError(f"File exceeds ordinary GitHub file limit: {rel}")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(data):
                credentials.append({"path": rel, "pattern": label})
        if rel_path.suffix.lower() in {".py", ".pyw"}:
            try:
                ast.parse(data, filename=rel)
                checked_python += 1
            except (SyntaxError, UnicodeError) as exc:
                syntax_failures.append({"path": rel, "error": str(exc)})
        if rel_path.suffix.lower() == ".ipynb":
            json.loads(data.decode("utf-8-sig"))
            checked_notebooks += 1
        files.append({"path": rel, "bytes": len(data), "sha256": sha256(data), "provenance": provenance, "mode": mode})

    for name in sorted(ROOT_FILES):
        path = source / name
        if path.is_file():
            write_file(name, path.read_bytes(), "workspace_copy")
    for tree in TREES:
        if not (source / tree).exists():
            continue
        for p in retained_files(source / tree):
            rel = p.relative_to(source).as_posix()
            write_file(rel, p.read_bytes(), "workspace_copy", "100755" if p.suffix.lower() == ".sh" else "100644")

    # Preserve official tracked files, even when the original Git tree tracked
    # bytecode. Restore those blobs from the fixed Git objects, not new caches.
    upstream = []
    for name, record in upstream_records(source).items():
        local = source / "upstream" / name
        count = 0
        for item in record["files"]:
            rel = item["path"]
            if (local / ".git").exists():
                data = subprocess.check_output(["git", "-C", str(local), "show", f"{record['commit']}:{rel}"])
            else:
                data = (local / rel).read_bytes()
            if sha256(data) != item["sha256"]:
                raise RuntimeError(f"Upstream checksum mismatch: {name}/{rel}")
            write_file(f"upstream/{name}/{rel}", data, f"upstream:{record['commit']}", item.get("mode", "100644"))
            count += 1
        upstream.append({"name": name, "repository": record["official_repository"], "commit": record["commit"], "verified_tracked_files": count})
    for p in sorted((source / "upstream").iterdir()):
        if p.is_file() and p.suffix.lower() in {".md", ".json"}:
            write_file(p.relative_to(source).as_posix(), p.read_bytes(), "workspace_provenance_copy")

    # An independent inventory prevents accidentally dropping model/helper code.
    coverage = []
    for tree in TREES:
        if not (source / tree).exists():
            continue
        for p in retained_files(source / tree):
            if p.suffix.lower() in SOURCE_EXTENSIONS:
                rel = p.relative_to(source).as_posix()
                dest = output / rel
                if not dest.is_file() or sha256(dest.read_bytes()) != sha256(p.read_bytes()):
                    coverage.append(rel)
    if coverage:
        raise RuntimeError("Missing or changed source files: " + json.dumps(coverage))
    if credentials:
        # Report paths and detector labels only; never print matched values.
        raise RuntimeError("Review potential credentials before staging: " + json.dumps(credentials))
    if syntax_failures:
        raise RuntimeError("Static Python syntax review required: " + json.dumps(syntax_failures))

    manifest = {
        "archive_date": "2026-11-12", "repository": "xk271521-droid/DGAT",
        "visibility_requested": "private", "source_root": str(source),
        "policy": "Research files are copied unchanged; fixed upstream tracked files are verified against SHA-256 manifests. Data, weights, prediction arrays, local environments, caches and secrets are excluded except original upstream tracked cache blobs.",
        "file_count": len(files), "total_bytes": sum(f["bytes"] for f in files),
        "counts_by_top_directory": dict(sorted(Counter(f["path"].split("/")[0] for f in files).items())),
        "upstream": upstream,
        "verification": {"copy_hashes": "passed", "retained_source_coverage": "passed", "python_syntax_files": checked_python, "notebook_json_files": checked_notebooks, "credential_pattern_matches": 0, "research_code_executed": False, "training_rerun": False},
        "files": sorted(files, key=lambda f: f["path"]),
    }
    (output / "GIT_SOURCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    # Raw NUL-delimited path list is local-only and used for precise git staging.
    # Original author-tracked bytecode is intentionally force-added by this list.
    stage = b"".join((f["path"].encode("utf-8") + b"\0") for f in manifest["files"])
    stage += b"GIT_SOURCE_MANIFEST.json\0"
    (output / ".git-stage-paths").write_bytes(stage)
    print(json.dumps({k: v for k, v in manifest.items() if k != "files"}, ensure_ascii=True, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    export(args.source_root, args.output_root)


if __name__ == "__main__":
    main()

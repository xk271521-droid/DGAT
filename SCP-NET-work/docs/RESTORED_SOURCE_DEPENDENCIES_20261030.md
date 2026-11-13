# Restored SCP-NET source dependencies

On 2026-10-30 the missing helper files were recovered from the user's supplied
WeChat file directory and copied into `ST2P/src` with their contents unchanged.
The parenthesized WeChat filenames were normalized to the module names used by
the retained code.

| Project path | Supplied source | SHA-256 |
|---|---|---|
| `ST2P/src/triscale_center_protein_v1.py` | `triscale_center_protein_v1(1).py` | `ADA7E1ABC59E21AE928BC5EEE5E43FA2D6721786492CD5DBA1E1CCEA0E362671` |
| `ST2P/src/triscale_center_protein_v5_program_decoder.py` | `triscale_center_protein_v5_program_decoder.py` | `34BBFF2178A35D33CBFEE4DEBA7FC23AC66329D3F920F4C16125DCCF8883A098` |
| `ST2P/src/triscale_center_protein_v7_anchored_cross_attention.py` | `triscale_center_protein_v7_anchored_cross_attention(1).py` | `12C9C2A098C74B43EBDEB80974CF43C1B2E3FC53C5CF8294ECFC9C35EC35C286` |

The files form one dependency chain: V2 imports V1; V5 and V7 import V2; and
the retained V9-compatible evaluator imports V2, V5, and V7. All restored
modules and the bounded-anchor candidate pass Python compilation and import
checks in the `st2p` environment.

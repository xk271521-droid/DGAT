import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle

OUT = r'E:\DGAT\figures'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5})
fig, ax = plt.subplots(figsize=(18,9.5), dpi=220)
ax.set_xlim(0,18); ax.set_ylim(0,10); ax.axis('off')
fig.patch.set_facecolor('white'); ax.set_facecolor('white')

navy='#1F3A5F'; blue='#DCEBFA'; teal='#DCEFEA'; teal2='#4E9587'
orange='#F7E1C1'; orange2='#C97620'; purple='#E8DDF3'; purple2='#7A559D'
red='#F6D6D6'; red2='#AA4141'; dark='#1F2933'; muted='#5B6570'

def box(x,y,w,h,text,fc='white',ec=dark,lw=1.35,fs=9.4,weight='normal',ls='-'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.03,rounding_size=0.11',fc=fc,ec=ec,lw=lw,linestyle=ls))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color=dark,weight=weight,wrap=True)

def arrow(x1,y1,x2,y2,color=dark,lw=1.8,dash=False,rad=0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,
        linewidth=lw,color=color,linestyle='--' if dash else '-',connectionstyle=f'arc3,rad={rad}'))

ax.text(9,9.62,'ProteinConvex for spatial RNA-to-ADT prediction',ha='center',va='center',fontsize=16,weight='bold',color=navy)
ax.text(9,9.30,'Overall framework: main prediction network → checkpoint selection → frozen output calibration',ha='center',va='center',fontsize=9.5,color=muted)

panels=[(0.35,1.05,3.35,7.85,'A  Data & context',blue,navy),
        (3.95,1.05,6.45,7.85,'B  Main predictor',teal,teal2),
        (10.65,1.05,3.0,7.85,'C  Selection',orange,orange2),
        (13.90,1.05,3.75,7.85,'D  Calibration & evaluation',purple,purple2)]
for x,y,w,h,t,fc,ec in panels:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=0.14',fc=fc,ec=ec,lw=1.9))
    ax.text(x+0.16,y+h-0.23,t,ha='left',va='center',fontsize=10.5,weight='bold',color=ec)

# A
box(0.70,7.30,2.65,0.92,'RNA input\n9384 common genes',ec=navy,fs=9.7,weight='bold')
for i in range(4):
    for j in range(6):
        c=['#C9DDF2','#8EB4DC','#E9F1F8'][(i*3+j)%3]
        ax.add_patch(Rectangle((0.94+j*0.13,7.55+i*0.12),0.105,0.095,fc=c,ec='white',lw=.35))
box(0.70,5.95,2.65,0.82,'Spatial coordinates',ec=navy,fs=9.7,weight='bold')
nodes=[(1.02,6.18),(1.44,6.41),(1.82,6.12),(1.37,5.85),(2.18,6.32)]
for a,b in [(0,1),(1,2),(0,3),(1,4),(2,4),(3,2)]:
    ax.plot([nodes[a][0],nodes[b][0]],[nodes[a][1],nodes[b][1]],color=navy,lw=1.0)
for x,y in nodes: ax.add_patch(Circle((x,y),0.075,fc='#AFC9E8',ec=navy,lw=.9))
box(0.70,4.55,2.65,0.82,'Spatial KNN\ncenter / local 4 / regional 4',ec=navy,fs=9.2,weight='bold')
box(0.70,2.98,2.65,0.92,'Slice context\nRNA summary + spatial scale',ec=navy,fs=9.2,weight='bold')
box(0.70,1.67,2.65,0.68,'Slice canonicalization',ec=navy,fs=9.2,weight='bold')
arrow(2.02,7.30,2.02,6.80,navy); arrow(2.02,5.95,2.02,5.42,navy); arrow(2.02,4.55,2.02,3.94,navy); arrow(2.02,2.98,2.02,2.65,navy)

# B
box(4.38,7.55,5.58,0.66,'RNA encoder + relative position',ec=teal2,fs=9.8,weight='bold')
box(4.38,6.34,5.58,0.70,'Protein-specific neighbour attention',ec=teal2,fs=9.8,weight='bold')
box(4.38,5.10,5.58,0.76,'ProteinConvex: 3-scale routing\ncenter / local / regional',ec=orange2,fs=9.7,weight='bold')
box(4.38,3.82,5.58,0.76,'Cognate RNA anchor\ncorresponding gene → protein',ec=purple2,fs=9.6,weight='bold')
box(4.38,2.54,5.58,0.76,'31 protein-specific heads\n→ 31 ADT values per spot',ec=teal2,fs=9.7,weight='bold')
box(4.38,1.70,5.58,0.42,'Slice context conditions routing and decoder adaptation',ec=teal2,lw=1.0,fs=8.2)
for y in [7.55,6.34,5.10,3.82]: arrow(7.17,y,7.17,y-0.45,teal2)
arrow(3.35,2.10,4.33,2.92,navy,rad=-0.08); arrow(3.35,3.44,4.33,7.88,navy,rad=-0.27)

# C
box(11.00,7.55,2.30,0.70,'Source fit\nbackpropagation',ec=orange2,fs=9.4,weight='bold')
box(11.00,6.30,2.30,0.70,'Source validation\nscore',ec=orange2,fs=9.4,weight='bold')
box(11.00,5.05,2.30,0.70,'Target RNA + space\nroughness',ec=orange2,fs=9.2,weight='bold')
box(11.00,3.52,2.30,0.86,'20 epochs\nchoose best checkpoint',ec=orange2,fs=9.4,weight='bold')
box(11.00,1.80,2.30,0.76,'Target ADT withheld\nuntil final evaluation',ec=red2,fs=8.8,weight='bold')
arrow(9.98,2.92,10.96,7.88,teal2,rad=-0.28)
arrow(12.15,7.55,12.15,7.05,orange2); arrow(12.15,6.30,12.15,5.76,orange2); arrow(12.15,5.05,12.15,4.48,orange2)
arrow(11.00,3.95,10.60,3.95,orange2,rad=0.2); arrow(10.60,3.95,10.60,7.88,orange2); arrow(10.60,7.88,10.96,7.88,orange2)

# D
box(14.25,7.55,3.05,0.70,'Freeze best checkpoint',ec=purple2,fs=9.7,weight='bold')
box(14.25,6.24,3.05,0.78,'Strict RNA\nRNA / spatial neighbour calibration',ec=purple2,fs=9.0,weight='bold')
box(14.25,4.86,3.05,0.78,'Source Affine\nprotein-wise scale + offset',ec=purple2,fs=9.0,weight='bold')
box(14.25,3.48,3.05,0.78,'Final target prediction\n31 ADT values per spot',ec=purple2,fs=9.0,weight='bold')
box(14.25,1.85,3.05,0.76,'Target ADT\nfinal metrics only',ec=red2,fs=9.2,weight='bold',ls='--')
arrow(13.64,3.95,14.20,7.90,orange2,rad=-0.18)
arrow(15.78,7.55,15.78,7.04,purple2); arrow(15.78,6.24,15.78,5.64,purple2); arrow(15.78,4.86,15.78,4.28,purple2); arrow(15.78,3.48,15.78,2.62,purple2); arrow(15.78,2.62,15.78,2.28,red2,dash=True)

# compact legend and workflow cues
ax.text(3.66,8.62,'forward pass',ha='center',va='center',fontsize=8.4,color=navy)
ax.text(10.48,8.62,'epoch-level choice',ha='center',va='center',fontsize=8.4,color=orange2)
ax.text(13.68,8.62,'frozen post-processing',ha='center',va='center',fontsize=8.4,color=purple2)
ax.add_line(plt.Line2D([0.50,0.86],[0.55,0.55],color=navy,lw=2)); ax.text(0.96,0.55,'forward / data flow',va='center',fontsize=8.3,color=muted)
ax.add_line(plt.Line2D([3.18,3.54],[0.55,0.55],color=red2,lw=1.6,linestyle='--')); ax.text(3.64,0.55,'target ADT only at final evaluation',va='center',fontsize=8.3,color=muted)
ax.text(17.55,0.55,'Overview · draft 2',ha='right',va='center',fontsize=8.2,color='#6B7280')
ax.add_patch(Rectangle((0.18,0.20),17.64,9.15,fill=False,ec='#AAB4BE',lw=0.8))

fig.savefig(os.path.join(OUT,'proteinconvex_overview_draft2.png'),dpi=300,bbox_inches='tight',facecolor='white')
fig.savefig(os.path.join(OUT,'proteinconvex_overview_draft2.svg'),bbox_inches='tight',facecolor='white')
plt.close(fig)
print('created')

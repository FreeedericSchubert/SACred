import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
F=Path('/mnt/data/hackathon_report_assets');F.mkdir(exist_ok=True)
NAVY='#10253C';BLUE='#316FA2';TEAL='#149C90';AMBER='#E39B35';CORAL='#DA7060';SLATE='#65758A';LIGHT='#D8E2EA'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def style(ax,axis='y'):
 ax.grid(axis=axis,color='#E6EDF3',linewidth=.7)
 ax.set_axisbelow(True)
 for name in ['left','bottom']:ax.spines[name].set_color('#B8C5D2')
 ax.tick_params(colors='#53647A')
def save(fig,path):fig.savefig(F/path,dpi=220,bbox_inches='tight',pad_inches=.15);plt.close(fig)
# Training: known epoch 0-7 only
E=[9.37,7.19,5.54,4.86,4.29,4.22,4.06,4.21,4.22]
R=[246.23,95.31,81.14,75.30,71.66,69.35,67.68,66.84,65.94]
fig,axs=plt.subplots(1,2,figsize=(11.0,3.6),layout='constrained')
for ax,vals,c,title,unit in [(axs[0],R,TEAL,'Validation force RMSE','meV/Å'),(axs[1],E,BLUE,'Validation energy RMSE','meV/atom')]:
 style(ax);ax.plot(range(-1,8),vals,'-o',lw=2.4,ms=5.7,color=c,zorder=5)
 ax.scatter([-1],[vals[0]],color=AMBER,s=55,zorder=6)
 ax.set_xticks([-1,0,1,2,3,4,5,6,7],['init','0','1','2','3','4','5','6','7']);ax.set_xlabel('Fine-tuning epoch');ax.set_ylabel(unit);ax.set_title(title,loc='left',fontweight='bold',color=NAVY)
 ax.annotate(f'{vals[-1]:.2f}',(7,vals[-1]),xytext=(-9,13),textcoords='offset points',ha='right',fontsize=10,fontweight='bold',color=c)
 ax.annotate(f'{vals[0]:.1f}',(-1,vals[0]),xytext=(8,-20),textcoords='offset points',fontsize=9,color=SLATE)
save(fig,'train_validation_curve.png')

# test bars errors force role
roles=['Graphene C\nfar','Graphene C\nnear','C near\nmetal','H','N','O','S','Metal']
vals=[53.70,113.02,222.77,75.79,121.93,145.27,111.07,334.23]
colors=[BLUE,BLUE,AMBER,SLATE,TEAL,TEAL,TEAL,CORAL]
fig,ax=plt.subplots(figsize=(11.3,4.15),layout='constrained');style(ax)
bars=ax.bar(roles,vals,color=colors,width=.7)
for p,v in zip(bars,vals):ax.text(p.get_x()+p.get_width()/2,v+8,f'{v:.0f}',ha='center',fontsize=9,fontweight='bold',color=NAVY)
ax.set_ylim(0,395);ax.set_ylabel('Force RMSE (meV/Å)');ax.set_title('Held-out test  •  error by atomic role',loc='left',fontweight='bold',color=NAVY)
save(fig,'test_errors_roles.png')

# test metals
metals=['Cu','Co','Ni','Fe','Zn'];v=[194.61,203.08,320.28,364.23,421.22];counts=[72,27,70,121,81]
fig,ax=plt.subplots(figsize=(9.7,3.6),layout='constrained');style(ax)
colors=[BLUE,TEAL,AMBER,CORAL,'#B74F62'];bars=ax.bar(metals,v,color=colors,width=.55)
for p,a,n in zip(bars,v,counts):ax.text(p.get_x()+p.get_width()/2,a+10,f'{a:.0f}',ha='center',fontweight='bold',color=NAVY);ax.text(p.get_x()+p.get_width()/2,15,f'n={n}',ha='center',fontsize=9,color='white')
ax.set_ylim(0,475);ax.set_ylabel('Force RMSE on metal (meV/Å)');ax.set_title('Averaged error hides metal-specific failures',loc='left',fontweight='bold',color=NAVY)
save(fig,'test_errors_metals.png')

# DFT vs xTB along xTB NEB images
d=pd.read_csv('/mnt/data/mace_1500K_ready/md_start/escape_paths/escape_path_energies.csv')
fig,axs=plt.subplots(1,2,figsize=(11.1,3.75),layout='constrained',sharex=True)
for ax,(s,title) in zip(axs,[('Zn-NS2:SV','Zn–NS₂'),('Ni-NS2:SV','Ni–NS₂')]):
 b=d.loc[d.site==s].sort_values('image');style(ax)
 ax.plot(b.image,b.E_xtb_rel,marker='o',lw=2.2,color=TEAL,label='GFN2-xTB NEB')
 ax.plot(b.image,b.E_dft_rel,marker='s',lw=2.2,ms=5,color=CORAL,label='DFT single points')
 ax.set_xticks(list(range(9)));ax.set_xlabel('Image along xTB NEB path');ax.set_ylabel('Relative energy (eV)');ax.set_title(title,loc='left',color=NAVY,fontweight='bold');ax.legend(frameon=False,fontsize=8.8,loc='upper left')
save(fig,'xtb_dft_paths.png')

# separate MD coordination (real-frame values and log observations)
dist=pd.read_csv('/mnt/data/pyrolysis_analysis/metal_donor_distances.csv');fig,axs=plt.subplots(1,3,figsize=(13,4.0),layout='constrained')
for ax in axs:style(ax);ax.set_ylim(-.2,4.5);ax.set_yticks([0,1,2,3,4]);ax.set_xlabel('MD time (ps)')
zn=dist[dist.site.str.startswith('Zn')].groupby('time_ps').bonded.sum(); ni=dist[dist.site.str.startswith('Ni')].groupby('time_ps').bonded.sum()
axs[0].step(zn.index,zn.values,where='post',lw=2.5,color=CORAL);axs[0].axvline(.15,ls='--',lw=1.2,color=AMBER);axs[0].axvline(1.15,ls=':',lw=1.2,color=SLATE);axs[0].set_title('Zn–NS₂  |  xTB bonded start',loc='left',fontweight='bold',color=NAVY);axs[0].set_xlim(0,3.15)
axs[0].text(0.24,3.78,'loss at 0.15 ps',color=CORAL,fontsize=9)
axs[1].axvspan(0,3.1,color='#F2F5F8');axs[1].plot([0.25,.5,1,1.5,2,2.5,3],[3]*7,'o',color=BLUE,ms=4,alpha=.7);axs[1].step(ni.index,ni.values,where='post',lw=2.5,color=BLUE);axs[1].set_xlim(0,20);axs[1].set_title('Ni–NS₂  |  log + XYZ',loc='left',fontweight='bold',color=NAVY);axs[1].text(5,3.78,'3 → 2 donors',color=BLUE,fontsize=9)
coT=[.25,.5,1,1.5,2,2.5,3,3.5,4,4.5,5];axs[2].scatter(coT,[4]*len(coT),s=31,color=TEAL);axs[2].set_xlim(0,5.15);axs[2].set_title('Co–N₂OC  |  logs only',loc='left',fontweight='bold',color=NAVY);axs[2].text(.5,3.43,'4 / 4 retained',color=TEAL,fontsize=9)
axs[0].set_ylabel('Donor contacts within cutoffs')
save(fig,'md_coordination_three.png')

# Metal donor bond distance 2 panels; Co not included
fig,axs=plt.subplots(1,2,figsize=(11.4,4.25),layout='constrained')
for ax,(site,title,tmax) in zip(axs,[('Zn','Zn–NS₂ · 0–3.15 ps',3.15),('Ni','Ni–NS₂ · 3.1–20 ps',20.0)]):
 style(ax);ds=dist[dist.site.str.startswith(site)]
 for ix,(idx,part) in enumerate(ds.groupby('donor_index')):
  elem=part.donor_element.iloc[0];c=[BLUE,AMBER,CORAL][ix]
  ax.plot(part.time_ps,part.distance_A,lw=1.8,color=c,label=f'{elem} ({idx})');ax.axhline(part.cutoff_A.iloc[0],color=c,alpha=.5,lw=.8,ls='--')
 ax.set_xlim(0,tmax);ax.set_ylim(1.8, min(18,ds.distance_A.max()+1));ax.set_title(title,loc='left',fontweight='bold',color=NAVY);ax.set_xlabel('Simulation time (ps)');ax.set_ylabel('Periodic metal–donor distance (Å)');ax.legend(frameon=False,fontsize=8,ncol=3)
save(fig,'md_donor_distances.png')

# quick content schematic: circles, not real atomic positions labelled illustrative
from matplotlib.patches import FancyBboxPatch
fig,ax=plt.subplots(figsize=(10.7,1.7));ax.set_xlim(0,12);ax.set_ylim(0,2);ax.axis('off')
cols=[NAVY,'#1A5276',BLUE,TEAL,'#168676',CORAL]
labels=['PDB motif\nmining','Graphene\nsite building','xTB / DFT\nscreening','DFT trajectories\nE and forces','MACE v2\nfine-tuning','1500 K\nexploratory MD']
for i,l in enumerate(labels):
 x=i*2.06;ax.add_patch(FancyBboxPatch((x,.4),1.72,1.2,boxstyle='round,pad=.02,rounding_size=.16',facecolor=cols[i],edgecolor='none'));ax.text(x+.86,1,l,color='white',ha='center',va='center',fontsize=9,fontweight='bold')
 if i<5:ax.annotate('',xy=(x+2.00,1),xytext=(x+1.77,1),arrowprops={'arrowstyle':'-|>','lw':1.3,'color':'#8BA0B3'})
save(fig,'pipeline.png')
print('Created plots:',[p.name for p in F.glob('*.png')])

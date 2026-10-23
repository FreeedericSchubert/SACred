#!/usr/bin/env python3
"""Parse uploaded joined extended XYZ trajectories and generate transparent 1500 K MD plots.
The uploaded 'pyrolysis_results' file contains complete Zn MD frames, then a
fragmentary Ni frame, then complete Ni frames for 3.10-20.00 ps; Co absent.
"""
import re, zipfile, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ARCHIVE=Path('/mnt/data/res.zip')
DEST=Path('/mnt/data/pyrolysis_analysis')
DEST.mkdir(exist_ok=True)

SITE_NAMES={'Zn-NS2:SV':'Zn–NS₂', 'Ni-NS2:SV':'Ni–NS₂','Co-N2OC:DV:NONC':'Co–N₂OC'}
CUTOFF={'N':3.1,'O':3.1,'S':3.3,'C':2.8}

def parse(data):
    lines=data.decode('utf-8').splitlines()
    results=[]; skipped=0; i=0
    while i<len(lines):
        if not lines[i].strip().isdigit() or i+1>=len(lines) or not lines[i+1].startswith('Lattice='):
            skipped+=1; i+=1;continue
        count=int(lines[i])
        if i+count+2>len(lines):break
        header=lines[i+1]
        m=re.search(r'motif=([^\s]+)',header)
        tt=re.search(r'md_time_ps=([\d.]+)',header)
        metal=re.search(r'metal_index=([\d]+)',header)
        donor=re.search(r'donor_indices="([\d,\s]+)"',header)
        cell=re.search(r'Lattice="([^"]+)"',header)
        if not all((m,tt,metal,donor,cell)):
            skipped+=1; i+=1;continue
        body=[row.split() for row in lines[i+2:i+2+count]]
        if not all(len(row)>=4 for row in body):
            skipped+=1;i+=1;continue
        try:
            species=[row[0] for row in body]
            positions=np.array([[float(x) for x in row[1:4]] for row in body])
            c=np.array([float(x) for x in cell.group(1).split()]).reshape(3,3)
            mid=int(metal.group(1)); ds=[int(x) for x in donor.group(1).replace(',',' ').split()]
            t=float(tt.group(1))
        except (ValueError,IndexError):
            skipped+=1;i+=1;continue
        # Exhaustive MIC for the small skewed periodic graphene cell: enumerate
        # nearby lattice translations around the nearest fractional offset.
        ds_vectors=positions[ds]-positions[mid]
        u=ds_vectors@np.linalg.inv(c)
        base=-np.rint(u).astype(int)
        offsets=np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1],indexing='ij')).reshape(3,-1).T
        shifts=base[:,None,:]+offsets[None,:,:]
        vectors=ds_vectors[:,None,:]+np.matmul(shifts,c)
        dist=np.linalg.norm(vectors,axis=-1).min(axis=1)
        cuts=[CUTOFF.get(species[d],3.1) for d in ds]
        results.append(dict(motif=m.group(1),time_ps=t,n_atoms=count,metal_index=mid,
                            donors=ds,species=[species[d] for d in ds],dist=dist.tolist(),
                            cuts=cuts, bound=int(sum(v<c for v,c in zip(dist,cuts))),
                            metal_z=float(positions[mid,2]),source='extxyz'))
        i+=2+count
    return results,skipped

with zipfile.ZipFile(ARCHIVE) as archive:
    original,skipped=parse(archive.read('pyrolysis_results'))
    bonded,skipped2=parse(archive.read('zn_bonded_results'))

# The original Zn trajectory started detached, so exclude it from the scientific plot.
ni=[r for r in original if r['motif']=='Ni-NS2:SV']
zn=[r for r in bonded if r['motif']=='Zn-NS2:SV']
original_zn=[r for r in original if r['motif']=='Zn-NS2:SV']
ni.sort(key=lambda r:r['time_ps']);zn.sort(key=lambda r:r['time_ps']);original_zn.sort(key=lambda r:r['time_ps'])
print('Parsed',len(original),'orig,',len(bonded),'bonded. Ignored unparsable lines',skipped,skipped2)
for name,rows in [('bonded Zn',zn),('Ni',ni),('original detached Zn',original_zn)]:
 if rows:
  print(name,'frames',len(rows),'t',rows[0]['time_ps'],rows[-1]['time_ps'],'first/last d',np.round(rows[0]['dist'],2),np.round(rows[-1]['dist'],2),'first/last count',rows[0]['bound'],rows[-1]['bound'])

# CSV of actual trajectory-extracted distances only.
with (DEST/'metal_donor_distances.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['site','time_ps','metal_index','donor_index','donor_element','distance_A','cutoff_A','bonded','source'])
    for site, rows in [('Zn–NS2 (bonded NEB start)',zn),('Ni–NS2 (partial trajectory)',ni)]:
        for r in rows:
            for idx,el,dd,cut in zip(r['donors'],r['species'],r['dist'],r['cuts']):
                w.writerow([site,r['time_ps'],r['metal_index'],idx,el,round(dd,6),cut,int(dd<cut),r['source']])

# Useful log reports from previous shared run, not from Co trajectory.
co_log_times=[0.25,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
ni_log_times=[0.25,0.5,1,1.5,2,2.5,3]
zn_event=0.15
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10.5,'axes.titlesize':11.5,'xtick.labelsize':9,'ytick.labelsize':9,'savefig.facecolor':'white','figure.facecolor':'white'})
fig=plt.figure(figsize=(12.3,9.2),layout='constrained')
gs=fig.add_gridspec(3,2,width_ratios=[1.33,1],height_ratios=[1,1,0.85],hspace=0.1)
axzn=fig.add_subplot(gs[0,0]);axni=fig.add_subplot(gs[1,0]);axco=fig.add_subplot(gs[2,0]);axcount=fig.add_subplot(gs[:,1])
fig.suptitle('1500 K pyrolysis MD  |  coordination persistence',fontsize=17,fontweight='bold',x=0.02,ha='left',y=1.04)
colors=['#3866A3','#D98329','#7B58A4','#2A9D8F']

def style(ax):
    ax.grid(axis='y',color='#DEE3E9',linewidth=0.6)
    ax.spines[['top','right']].set_visible(False)
    ax.spines['left'].set_color('#BEC6CE')
    ax.spines['bottom'].set_color('#BEC6CE')
    ax.set_axisbelow(True)
    ax.tick_params(colors='#394450')

for ax,rows,name,maxt in [(axzn,zn,'Zn–NS₂  |  bonded xTB NEB start',3.15),(axni,ni,'Ni–NS₂  |  3.1–20 ps available',20)]:
    style(ax)
    x=[r['time_ps'] for r in rows]
    for j,(donor,element,cut) in enumerate(zip(rows[0]['donors'],rows[0]['species'],rows[0]['cuts'])):
        y=[r['dist'][j] for r in rows]
        ax.plot(x,y,lw=1.65,color=colors[j],label=f'{element}{donor} (limit {cut:.1f} Å)')
        ax.axhline(cut,color=colors[j],alpha=0.4,ls=(0,(4,3)),lw=0.95)
    ax.set_xlim(0,maxt)
    ax.set_ylim(0,max(7.5, min(14, max(r['dist'][j] for r in rows for j in range(len(r['dist'])))+1)))
    ax.set_ylabel('Metal–donor distance (Å)')
    ax.set_title(name,loc='left',pad=8,weight='bold')
    ax.legend(loc='upper right',fontsize=7.6,ncol=3,frameon=False,columnspacing=0.6,handlelength=1.5)
    if ax is axni:
        ax.axvspan(0,3.1,facecolor='#F1F2F4',alpha=.8,zorder=-1)
        ax.text(1.55,ax.get_ylim()[1]*0.5,'Coordinates\nnot in archive',ha='center',va='center',size=9,color='#68727F')
        ax.set_xlabel('Simulation time (ps)')
    else:
        ax.axvline(zn_event,color='#BA4A34',lw=1.5,ls=':',label='event')
        ax.text(zn_event+0.12,ax.get_ylim()[1]*0.67,'<2 contacts at\n0.15 ps',color='#BA4A34',size=9,va='center')

style(axco)
axco.set_facecolor('#F9FAFC')
axco.set_title('Co–N₂OC  |  trajectory missing',loc='left',pad=8,weight='bold')
axco.set_xlim(0,5);axco.set_ylim(0,1)
axco.set_xticks([]);axco.set_yticks([])
axco.text(2.5,0.55,'Metal–donor distances unavailable',ha='center',fontsize=10,color='#465364')
axco.text(2.5,0.27,'Log: 4/4 contacts retained through 5 ps',ha='center',fontsize=9.5,color='#2A806A')

style(axcount)
axcount.set_title('Coordination contacts retained',loc='left',pad=8,weight='bold')
axcount.set_xlim(-0.3,20.5);axcount.set_ylim(-0.25,4.65)
axcount.set_yticks([0,1,2,3,4]);axcount.set_xlabel('Simulation time (ps)');axcount.set_ylabel('Donor contacts (distance cutoffs)')
# actual frame-derived Zn and Ni
axcount.step([r['time_ps'] for r in zn],[r['bound'] for r in zn],where='post',color='#D98329',lw=2.15,label='Zn–NS₂ (bonded start)')
axcount.step([r['time_ps'] for r in ni],[r['bound'] for r in ni],where='post',color='#3866A3',lw=2.15,label='Ni–NS₂ (XYZ ≥3.1 ps)')
# Explicitly show log-only intervals with discrete markers, no interpolated trajectory
axcount.plot(ni_log_times,[3]*len(ni_log_times),'o',ms=4.2,alpha=.6,color='#3866A3',label='Ni log samples (<3.1 ps)')
axcount.plot(co_log_times,[4]*len(co_log_times),'s',ms=4.7,color='#2A9D8F',label='Co log samples (no XYZ)')
axcount.axvline(zn_event,ls=':',color='#D98329',alpha=.7)
axcount.text(0.7,0.35,'Zn escape criterion\nmet at 0.15 ps',fontsize=9,color='#AB6215',va='bottom')
axcount.legend(frameon=False,loc='center right',fontsize=8.5)
fig.text(.015,-0.01,'Source: uploaded res.zip (atom positions for bonded Zn and Ni 3.1–20 ps); Co contact counts and early Ni counts from previously shared Modal logs.\nDistances use 3D periodic minimum-image convention; contacts defined by original script thresholds (N/O 3.1 Å, S 3.3 Å, C 2.8 Å). One seed per site; exploratory MLIP MD.',fontsize=8.5,color='#556070')
fig.savefig(DEST/'pyrolysis_1500K_coordination.png',dpi=230,bbox_inches='tight',pad_inches=.24)
fig.savefig(DEST/'pyrolysis_1500K_coordination.svg',bbox_inches='tight',pad_inches=.24)
plt.close(fig)

summary={
  'source_file':'res.zip',
  'bonded_zn':{'frames':len(zn),'range_ps':[zn[0]['time_ps'],zn[-1]['time_ps']], 'initial_distances_A':zn[0]['dist'],'initial_contacts':zn[0]['bound'], 'final_contacts':zn[-1]['bound'],'escape_start_ps_from_log':0.15,'escape_confirm_ps_from_log':1.15,'zero_contacts_at_ps_from_saved_frames':0.25},
  'ni':{'frames':len(ni),'range_ps':[ni[0]['time_ps'],ni[-1]['time_ps']],'first_contacts':ni[0]['bound'],'final_contacts':ni[-1]['bound'],'missing_early_frames':True},
  'co':{'xyz_in_uploaded_archive':False,'log_contacts':4,'log_time_range_ps':[0.25,5.0]},
  'original_zn':{'frames':len(original_zn),'comment':'Original Zn DFT starting geometry already detached; excluded from comparison.'},
  'limitations':['Single simulation seed per site','Zn starts at bonded xTB NEB image, not DFT-minimized minimum','Ni coordinates before 3.1 ps missing','Co complete trajectory not included in ZIP','MACE trained mostly on relaxation data; metal escape energetics not validated']
}
(DEST/'analysis_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2,ensure_ascii=False))
print('Files',[(p.name,p.stat().st_size) for p in DEST.iterdir()])

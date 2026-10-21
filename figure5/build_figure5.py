"""Rebuild Figure 5 from the supplied aggregate results, without model refitting.
Requires Python, matplotlib, numpy and scipy. Run from any working directory.
"""
from pathlib import Path
import csv, json, math, hashlib
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'source_data'
def read(name):
    with (DATA/name).open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f))
a = read('Fig5A_ADNC_proportions.csv')
am = read('Fig5A_ADNC_model_results.csv')[0]
b = read('Fig5B_34protein_results.csv')
bs = read('Fig5B_replication_summary.csv')[0]
c = read('Fig5C_MMSE_predicted_trajectories.csv')
cn = read('Fig5C_MMSE_counts_by_stage.csv')
cs = [r for r in read('Fig5C_MMSE_stage_slopes.csv') if r['set']=='Primary']
cg = read('Fig5C_MMSE_global_test.csv')[0]
x = np.array([float(r['KROAD_beta']) for r in b])
y = np.array([float(r['pittsburgh_beta']) for r in b])
p = np.array([float(r['P']) for r in b])
q = stats.false_discovery_control(p, method='bh')
con = np.sign(x)==np.sign(y)
rep = con & (q<0.05)
rho, rp = stats.spearmanr(x,y)
assert len(b)==34 and sum(con)==27 and sum(rep)==16
assert np.allclose(q,[float(r['q']) for r in b],rtol=1e-10)
assert abs(rho-float(bs['Spearman_rho']))<1e-12
assert abs(rp-float(bs['Spearman_P']))<1e-12
assert sum(int(r['participants']) for r in cn)==int(cg['n_participants'])==323
assert sum(int(r['assessments']) for r in cn)==int(cg['n_assessments'])==1386
assert sum(int(r['n']) for r in a)==399
assert float(am['P'])<0.001 and float(cg['P'])<0.001
for s in range(5):
    cc=[r for r in c if r['stage']==f'Stage {s}']
    slope=(float(cc[-1]['predicted_MMSE'])-float(cc[0]['predicted_MMSE']))/4
    assert abs(slope-float(cs[s]['annual_MMSE_slope']))<1e-10
check={'all_34_spearman_rho':float(rho),'all_34_spearman_P':float(rp),
 'directional_concordance':int(sum(con)),'BH_FDR_concordant_replicated':int(sum(rep)),
 'one_sided_exact_binomial_P':stats.binomtest(int(sum(con)),34,0.5,alternative='greater').pvalue,
 'replicated_proteins':[r['protein'] for r,k in zip(b,rep) if k],
 'ADNC_plot_n':399,'ADNC_model_n':int(am['n']),'MMSE_participants':323,'MMSE_assessments':1386,
 'checks':'BH adjustment, coefficient sign concordance, correlation, stage counts and trajectory slopes independently reconciled.'}
(ROOT/'verification.json').write_text(json.dumps(check,indent=2),encoding='utf-8')

plt.rcParams.update({'font.family':'Arial','font.size':7.5,'axes.labelsize':8,
 'axes.titlesize':9,'axes.linewidth':0.6,'xtick.labelsize':7,'ytick.labelsize':7,
 'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':3,'ytick.major.size':3,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'dejavusans',
 'savefig.facecolor':'white'})
def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(direction='out',pad=3)
def panel_a(ax):
    clean(ax)
    bottom=np.zeros(5)
    for cat,col in zip(['Not/Low','Intermediate','High'],['white','0.65','0.18']):
        rr=[next(r for r in a if r['stage']==f'Stage {s}' and r['ADNC']==cat) for s in range(5)]
        vals=np.array([100*int(r['n'])/int(r['stage_n_with_ADNC']) for r in rr])
        ax.bar(range(5),vals,bottom=bottom,width=.69,color=col,edgecolor='black',linewidth=.5,label=cat)
        bottom+=vals
    ax.set(ylim=(0,100),xlim=(-.6,4.6),ylabel='Participants (%)',xlabel='Plasma stage')
    ax.set_yticks([0,25,50,75,100])
    ns=[next(r['stage_n_with_ADNC'] for r in a if r['stage']==f'Stage {s}') for s in range(5)]
    ax.set_xticks(range(5),[f'{s}\n$n$ = {n}' for s,n in enumerate(ns)])
    ax.text(0,1.19,'Adjusted OR per stage = 4.73\n95% CI 3.71–6.05; $P$ < 0.001',transform=ax.transAxes,va='bottom',linespacing=1.6)
    ax.legend(title='ADNC',loc='lower left',bbox_to_anchor=(-.025,1.02),ncol=3,frameon=False,
              handlelength=1,handletextpad=.4,columnspacing=.8,borderpad=0,title_fontsize=7,fontsize=7)
def panel_b(ax):
    clean(ax)
    ax.axhline(0,color='0.82',lw=.5,zorder=0)
    ax.axvline(0,color='0.82',lw=.5,zorder=0)
    ax.plot([-.14,.49],[-.14,.49],ls=(0,(4,3)),color='0.6',lw=.7,zorder=0)
    for mask,color,label in [(rep,'black','Replicated'),(con & ~rep,'0.65','Concordant, $q$ ≥ 0.05'),(~con,'white','Opposite direction')]:
        ax.scatter(x[mask],y[mask],s=19,facecolors=color,edgecolors='black',linewidths=.5,label=label,zorder=3)
    offsets={'GFAP':(-6,7),'ACHE':(7,-4),'GAP43':(8,10),'NfL':(6,-10),'BACE1':(-35,10),'SPP1':(8,-20)}
    for name,(dx,dy) in offsets.items():
        i=next(i for i,r in enumerate(b) if r['protein']==name)
        ax.annotate(name,(x[i],y[i]),xytext=(dx,dy),textcoords='offset points',fontsize=7,
         ha='right' if name=='GFAP' else 'left',va='center',
         arrowprops=dict(arrowstyle='-',color='0.35',lw=.4,shrinkA=2,shrinkB=3))
    ax.set(xlim=(-.14,.49),ylim=(-.14,.49),xlabel='K-ROAD stage β (SD per stage)',ylabel='Pittsburgh stage β (SD per stage)')
    ax.set_xticks([-.1,0,.1,.2,.3,.4]); ax.set_yticks([-.1,0,.1,.2,.3,.4])
    ax.text(0,1.07,'Spearman ρ = 0.67, $P$ < 0.001\n27/34 concordant (79.4%)\n16/34 replicated (BH-FDR $q$ < 0.05)',
            transform=ax.transAxes,va='bottom',linespacing=1.45)
    ax.legend(loc='lower right',frameon=False,fontsize=6.8,handletextpad=.2,labelspacing=.45,borderpad=.1)
def panel_c(ax):
    clean(ax)
    colors=['0.60','0.45','0.30','0.18','0.0']
    styles=['solid',(0,(5,2)),(0,(4,1.5,1,1.5)),(0,(1,1.6)),(0,(7,2,2,2))]
    ax.axvline(0,color='.72',ls=(0,(3,3)),lw=.65)
    for s in range(5):
        rr=[r for r in c if r['stage']==f'Stage {s}']
        t=np.array([float(r['time_years']) for r in rr])
        mu=np.array([float(r['predicted_MMSE']) for r in rr])
        lo=np.array([float(r['ci_low']) for r in rr]); hi=np.array([float(r['ci_high']) for r in rr])
        ax.fill_between(t,lo,hi,color=colors[s],alpha=.065,lw=0)
        ax.plot(t,mu,color=colors[s],ls=styles[s],lw=1.2,label=f'Stage {s} ($n$ = {cn[s]["participants"]})')
    ax.set(xlim=(-1,3),ylim=(5,30),xlabel='Years from index plasma',ylabel='Predicted MMSE')
    ax.set_xticks([-1,0,1,2,3]); ax.set_yticks([5,10,15,20,25,30])
    ax.text(0,1.06,'Global stage × time: $P$ < 0.001',transform=ax.transAxes,va='bottom')
    ax.legend(loc='upper left',bbox_to_anchor=(1.03,1),frameon=False,handlelength=3.3,fontsize=7.5,borderaxespad=0,labelspacing=.9)
    ax.text(1.03,.41,'323 participants\n1,386 assessments',transform=ax.transAxes,va='top',linespacing=1.6)

def save(fig,stem):
    fig.savefig(ROOT/f'{stem}.pdf')
    fig.savefig(ROOT/f'{stem}.svg')
    fig.savefig(ROOT/f'{stem}.png',dpi=600)

fig=plt.figure(figsize=(183/25.4,178/25.4))
aa=fig.add_axes([.095,.565,.315,.285]); bb=fig.add_axes([.57,.565,.39,.285]); cc=fig.add_axes([.095,.09,.66,.31])
panel_a(aa);panel_b(bb);panel_c(cc)
for text,title,xx,yy in [('a','Neuropathological validation',.025,.98),('b','Protein replication',.48,.98),('c','Longitudinal cognition',.025,.465)]:
    fig.text(xx,yy,text,fontweight='bold',fontsize=12,va='top')
    fig.text(xx+.035,yy-.001,title,fontsize=9,va='top')
save(fig,'Figure5_final')
plt.close(fig)
for tag,fn in [('A',panel_a),('B',panel_b),('C',panel_c)]:
    f=plt.figure(figsize=(4.6,4.5) if tag!='C' else (6.7,3.8))
    ax=f.add_axes([.16,.18,.77,.56] if tag!='C' else [.11,.17,.64,.65])
    fn(ax); save(f,f'Figure5{tag}');plt.close(f)
print(json.dumps(check,indent=2))

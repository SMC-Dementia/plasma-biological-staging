"""Generate a synthetic K-ROAD style demo dataset (150 participants, 30 per stage)."""
import argparse,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from biological_staging.staging import assign_stage
def main():
 a=argparse.ArgumentParser();a.add_argument('outdir');args=a.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
 rng=np.random.default_rng(20261021);n=150;stage=np.repeat(np.arange(5),30);latent=stage+rng.normal(0,.2,n)
 b=pd.DataFrame({'participant_id':['SYN%04d'%i for i in range(n)],'stage':stage,'age':rng.normal(70,7,n),'sex':np.where(rng.random(n)<.6,'F','M'),'assay_phase':np.where(np.arange(n)%2,'I','II'),'centiloid':latent*25+rng.normal(0,10,n),'education_years':rng.normal(12,3,n),'apoe4_carrier':rng.integers(0,2,n),'diagnosis':np.where(stage<2,'CU',np.where(stage<4,'MCI','dementia')),'mmse':30-stage*2,'cdr_sb':stage.astype(float)})
 b['BD-p-tau217']=np.where(stage==0,9.,10.5)+rng.uniform(-.1,.1,n)
 b['eMTBR-tau']=np.choose(stage,[9.,9.,10.,11.,12.])+rng.uniform(-.1,.1,n)
 for site in [181,212,231]: b['BD-p-tau'+str(site)]=latent+rng.normal(0,.6,n)
 for site in [181,205,212,217,231]: b['p-tau'+str(site)]=latent+rng.normal(0,.8,n)
 b['BD-p-tau205']=latent+rng.normal(0,.5,n)
 selected=pd.read_csv(ROOT/'metadata/stage_associated_34.csv').protein.tolist()
 for protein in (ROOT/'metadata/non_tau_197.txt').read_text(encoding='utf-8').splitlines(): b[protein]=rng.normal(0,1,n)+(.5*stage if protein in selected else 0)
 b['amyloid_pet_positive']=(b.centiloid>24).astype(int);b.to_csv(out/'baseline.csv',index=False)
 t=b.copy();t['tau_temporal_meta_suvr']=np.choose(stage,[1.,1.,1.5,1.5,1.5]);t['tau_neoc_suvr']=np.choose(stage,[1.,1.,1.1,1.4,1.8]);t.to_csv(out/'tau.csv',index=False)
 intercept=rng.normal(0,.4,n);slope=rng.normal(0,.12,n);cs=[];ss=[]
 for time in [0.,1.,2.,3.]:
  c=b[['participant_id','stage']].copy();c['years_from_baseline']=time;c['cdr_sb']=stage+intercept+(.2+.25*stage+slope)*time+rng.normal(0,.2,n);cs.append(c)
 for time in [0.,1.,2.]:
  s=b.copy();s['years_from_baseline']=time;s['eMTBR-tau']+=time*np.choose(np.arange(n)%3,[.7,0.,-.7])
  for protein in selected: s[protein]+=time*(.06+rng.normal(0,.04,n))+rng.normal(0,.15,n)
  s=assign_stage(s,'BD-p-tau217','eMTBR-tau');ss.append(s)
 pd.concat(cs).to_csv(out/'clinical.csv',index=False);pd.concat(cs).to_csv(out/'coupling_clinical.csv',index=False);pd.concat(ss).to_csv(out/'serial.csv',index=False)
 print('Synthetic data: 150 participants, 600 clinical assessments, 450 serial draws')
if __name__=='__main__':main()

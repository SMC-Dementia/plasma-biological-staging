"""Generate a synthetic PITT-ADRC style demo dataset (300 participants, 1,200 MMSE assessments)."""
from pathlib import Path
import random,csv,re,json,datetime
ROOT=Path(__file__).resolve().parent
rng=random.Random(20261021)
src=(ROOT.parent/'02_34protein_replication.R').read_text(encoding='utf-8')
block=src.split('KROAD_BETA <- c(',1)[1].split('\n)',1)[0]
betas={a.strip().strip('`'):float(b) for a,b in re.findall(r'([^,\n]+?)=([-0-9.]+)',block)}
plasma=[];clinical=[]
for i in range(300):
 s=i//60;pid='SYNTHETIC_%04d'%(i+1);dt=datetime.date(2018,1,1)+datetime.timedelta(days=i)
 age=65+18*rng.random();sex='Female' if i%2 else 'Male'
 p=9.0+rng.random()*.4 if s==0 else 9.8+rng.random()*.5
 e=[9.1,9.1,9.9,10.8,11.9][s]+rng.random()*.2
 latent=.5*s+rng.gauss(0,1.4)
 row={'PTID':pid,'sample_id':pid+'_P1','draw_date':dt.isoformat(),'DateDeath':(dt+datetime.timedelta(days=730+i%500)).isoformat(),
      'neuro220_BD_pTau_217':p,'multiplex_MTBR':e,'blood_draw_age':age,'sex':sex,
      'neuro220_plateID':'SYN_PLATE_'+str(i%3),
      'ADNC':'Not/Low' if latent<.5 else 'Intermediate' if latent<1.5 else 'High'}
 for name,b in betas.items():row[name]=10+b*s+rng.gauss(0,.8)
 plasma.append(row)
 intercept=27-1.3*s+rng.gauss(0,1.3);slope=-.3-.4*s+rng.gauss(0,.35)
 for days in [-180,0,365,730]:
  clinical.append({'PTID':pid,'MMSE_date':(dt+datetime.timedelta(days=days)).isoformat(),'MMSE':round(max(0,min(30,intercept+slope*days/365.25+rng.gauss(0,.6))),3)})
for name,rows in [('synthetic_plasma.csv',plasma),('synthetic_mmse.csv',clinical)]:
 with (ROOT/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(ROOT/'synthetic_design.json').write_text(json.dumps({'seed':20261021,'participants':300,'stages':{'Stage '+str(i):60 for i in range(5)},'MMSE_assessments':1200,'proteins':34,'data_origin':'Synthetic'},indent=2),encoding='utf-8')
print('Created synthetic demo: 300 participants; 1200 MMSE assessments; 34 proteins.')

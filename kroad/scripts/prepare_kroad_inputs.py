"""Map the K-ROAD database export files to the canonical input schema (kroad/metadata/schema_template.csv).
Inputs: baseline export (one row per participant), serial-draw export, clinical export, coupling-clinical export
and, optionally, the tau PET workbook. The output contains participant-level data; keep it outside the repository.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from biological_staging.io import canonicalize_columns
from biological_staging.staging import assign_stage
MAPPING={'RID': 'participant_id', 'BD_pTau_217': 'BD-p-tau217', 'BD_pTau_205': 'BD-p-tau205', 'BD_pTau_T205': 'BD-p-tau205', 'eMTBR_Tau': 'eMTBR-tau', 'Age': 'age', 'Sex': 'sex', 'Phase': 'assay_phase', 'LEC_yn': 'anti_amyloid_treatment', 'Klunk_Centiloid': 'centiloid', 'Diagnosis': 'diagnosis', 'MMSE': 'mmse', 'CDR.SOB': 'cdr_sb', 'CDR_SOB': 'cdr_sb', 'Education': 'education_years', 'e4': 'apoe4_carrier', 'time_yr': 'years_from_baseline', 'time_yr_from_draw1': 'years_from_baseline', 'Sample_ID': 'draw_id', 'pRAB10_T73': 'pRAB10-T73', 'NT_proBNP': 'NT-proBNP', 'pLRRK2_S1292': 'pLRRK2-S1292', 'pQ_HTT': 'pQ-HTT', 'pSNCA_129': 'pSNCA-129', 'Oligo_SNCA': 'Oligo-SNCA', 'Ab40': 'Aβ40', 'pTDP43_409': 'pTDP43-409', 'BD_pTau_181': 'BD-p-tau181', 'BD_pTau_T212': 'BD-p-tau212', 'BD_pTau_231': 'BD-p-tau231', 'pTau_181': 'p-tau181', 'pTau_212': 'p-tau212', 'pTau_217': 'p-tau217', 'pTau_231': 'p-tau231', 'pTau_205': 'p-tau205'}
def prepare(baseline,serial,clinical,coupling,out,tau_workbook=None):
 out=Path(out)
 if out.exists() and any(out.iterdir()): raise ValueError('Use a new or empty output directory')
 out.mkdir(parents=True,exist_ok=True)
 def read(path):
  d=pd.read_csv(path,low_memory=False)
  # Analytes exported twice (rounded merge copy, then the full-precision panel copy): keep the panel copy.
  for col in [c for c in d.columns if c.endswith('.1') and c[:-2] in d.columns]:
   base=col[:-2]
   if not np.allclose(pd.to_numeric(d[base],errors='coerce'),pd.to_numeric(d[col],errors='coerce'),atol=1e-6,equal_nan=True): raise ValueError('Conflicting duplicate column '+base)
   d[base]=d[col]; d=d.drop(columns=[col])
  d=canonicalize_columns(d.rename(columns=MAPPING)).copy()
  if 'apoe4_carrier' in d:
   g=d.apoe4_carrier.astype('string');known=g.isin(['E2/E2','E2/E3','E2/E4','E3/E3','E3/E4','E4/E4'])
   d['apoe4_carrier']=g.str.contains('E4',regex=False).where(known).astype('Int64')
  if 'centiloid' in d: d['amyloid_pet_positive']=(d.centiloid>24).where(d.centiloid.notna()).astype('Int64')
  return d
 b=assign_stage(read(baseline),'BD-p-tau217','eMTBR-tau')
 if b.participant_id.duplicated().any(): raise ValueError('Baseline IDs must be unique')
 s=assign_stage(read(serial),'BD-p-tau217','eMTBR-tau')
 c=read(clinical).drop(columns=['stage'],errors='ignore').merge(b[['participant_id','stage']],on='participant_id',validate='many_to_one')
 c=c[c.anti_amyloid_treatment.eq(0)&c.stage.isin(range(5))]
 tables={'baseline':b,'serial':s,'clinical':c,'coupling_clinical':read(coupling)}
 if tau_workbook:
  t=pd.read_excel(tau_workbook,sheet_name='TauPET')
  t=b.merge(t,left_on='participant_id',right_on='RID',how='inner')
  t['distance']=(pd.to_datetime(t.Tau_PET_date)-pd.to_datetime(t.Blood_Collection_date)).dt.days.abs()
  if t.distance.isna().any(): raise ValueError('Missing PET/plasma date')
  t=t.sort_values(['participant_id','distance','Tau_PET_date']).drop_duplicates('participant_id')
  tables['tau']=t.rename(columns={'metaT_SMC':'tau_temporal_meta_suvr','NeoC':'tau_neoc_suvr'})
 for name,d in tables.items(): d.to_csv(out/(name+'.csv'),index=False)
 (out/'input_summary.json').write_text(json.dumps({k:{'rows':len(v),'participants':v.participant_id.nunique()} for k,v in tables.items()},indent=2),encoding='utf-8')
if __name__=='__main__':
 a=argparse.ArgumentParser()
 a.add_argument('--baseline',required=True);a.add_argument('--serial',required=True);a.add_argument('--clinical',required=True);a.add_argument('--coupling',required=True)
 a.add_argument('--tau-workbook');a.add_argument('--out',required=True)
 args=a.parse_args();prepare(args.baseline,args.serial,args.clinical,args.coupling,args.out,args.tau_workbook)

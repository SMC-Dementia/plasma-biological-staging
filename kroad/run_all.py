"""Run the K-ROAD analysis scripts on a directory of canonical input CSVs."""
import argparse,subprocess,sys,time,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data-dir',required=True);ap.add_argument('--outdir',required=True);ap.add_argument('--demo',action='store_true');ap.add_argument('--figures',action='store_true');a=ap.parse_args()
 d=Path(a.data_dir).resolve();o=Path(a.outdir).resolve()
 if o.exists() and any(o.iterdir()): raise ValueError('Use a new or empty output directory')
 o.mkdir(parents=True,exist_ok=True);tables=o/'tables';tables.mkdir();env=os.environ.copy();env['BIOSTAGING_MODEL_LOG']=str(o/'model_diagnostics.jsonl')
 jobs=[('03_stage_assignment.py',[str(d/'baseline.csv'),str(o/'baseline_staged.csv')]),('04_table1_baseline.py',[str(d/'baseline.csv'),str(tables/'Table1_baseline.csv')])]
 specs=[('01_T1_selection.py',['baseline.csv']),('02_T2_selection.py',['tau.csv']),('05_clinical_longitudinal.py',['clinical.csv']),('06_stage_transition.py',['serial.csv']),('08_crosssectional_proteome.py',['baseline.csv',str(ROOT/'metadata/non_tau_197.txt')]),('09_adjacent_stage_proteome.py',['baseline.csv']),('10_longitudinal_proteome.py',['serial.csv','baseline.csv']),('11_domain_analysis.py',['serial.csv','baseline.csv'])]
 for script,files in specs:
  args=[str(d/f) if not Path(f).is_absolute() else f for f in files]+['--outdir',str(tables)]
  if a.demo and script.startswith('01'):args+=['--n-boot','30']
  if a.demo and script.startswith('02'):args+=['--n-boot-auc','30','--n-boot-delta','30']
  jobs.append((script,args))
 jobs.append(('07_transition_clinical_coupling.py',[str(d/'coupling_clinical.csv'),str(tables/'Figure3B_participant_transitions.csv'),'--outdir',str(tables)]))
 if a.figures:jobs.append(('91_generate_figures.py',['--table-dir',str(tables),'--figure-dir',str(o/'figures')]))
 records=[]
 for script,args in jobs:
  t=time.monotonic();r=subprocess.run([sys.executable,str(ROOT/'scripts'/script),*args],env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
  (o/(script+'.log')).write_text(r.stdout+'\n'+r.stderr,encoding='utf-8');records.append({'script':script,'returncode':r.returncode,'seconds':round(time.monotonic()-t,3)})
  (o/'execution.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
  if r.returncode: raise RuntimeError('Analysis failed: '+script+'; see saved log')
 import shutil
 shutil.copy2(o/'stage_counts.csv',tables/'stage_counts.csv')
 if not a.demo:
  r=subprocess.run([sys.executable,str(ROOT/'scripts/99_reproducibility_QC.py'),'--table-dir',str(tables),'--output',str(o/'QC_REPORT.csv')],env=env)
  sys.exit(r.returncode)
 print('Demo completed.')
if __name__=='__main__':main()

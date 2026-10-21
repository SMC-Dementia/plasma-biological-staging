"""Run the three PITT-ADRC validation scripts with the fixed K-ROAD thresholds."""
import argparse, csv, json, pathlib, shutil, subprocess, sys, time

def read_csv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True,type=pathlib.Path)
    ap.add_argument('--clinical',required=True,type=pathlib.Path)
    ap.add_argument('--outdir',required=True,type=pathlib.Path)
    ap.add_argument('--rscript',default='Rscript')
    args,extra=ap.parse_known_args()
    allowed={'--sheet','--clinical-sheet','--clinical-id-col','--mmse-date-col','--mmse-col'}
    allowed.update('--'+v+'-col' for v in ['id','sample-id','sample-date','death-date','p217','emtbr','age','sex','plate','adnc'])
    if len(extra)%2 or any(extra[i] not in allowed for i in range(0,len(extra),2)):
        ap.error('Only column and sheet overrides are accepted.')
    for f in [args.input,args.clinical]:
        if not f.is_file():ap.error('Input file does not exist: '+str(f))
    rs=shutil.which(args.rscript)
    if not rs:ap.error('Rscript not found; supply --rscript with its full path.')
    out=args.outdir.resolve()
    if out.exists() and any(out.iterdir()):ap.error('Use a new or empty output directory.')
    out.mkdir(parents=True,exist_ok=True)
    src=pathlib.Path(__file__).resolve().parent
    start=time.perf_counter()
    expr='cat(R.version.string,"\\n"); print(Sys.info()[c("sysname","release","machine")]); print(sapply(c("MASS","nlme"),function(x)as.character(packageVersion(x)))); sessionInfo()'
    env=subprocess.run([rs,'--vanilla','-e',expr],capture_output=True,text=True,encoding='utf-8',errors='replace')
    (out/'sessionInfo.txt').write_text(env.stdout+env.stderr,encoding='utf-8')
    if env.returncode:raise RuntimeError('R dependency check failed; see sessionInfo.txt')
    for script in ['01_ADNC_validation.R','02_34protein_replication.R','03_MMSE_validation.R']:
        cmd=[rs,'--vanilla',str(src/script),'--input',str(args.input.resolve()),'--clinical',str(args.clinical.resolve()),'--outdir',str(out),*extra]
        t=time.perf_counter(); r=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace')
        (out/(script+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
        if r.returncode:raise RuntimeError(script+' failed; see its log')
        print(script+' completed in %.2f seconds'%(time.perf_counter()-t))
    b=read_csv(out/'Fig5B_34protein_results.csv')
    if len(b)!=34:raise RuntimeError('Primary manuscript analysis requires all 34 protein models; got '+str(len(b)))
    a=read_csv(out/'Fig5A_ADNC_model_results.csv')[0]
    g=read_csv(out/'Fig5C_MMSE_global_test.csv')[0]
    if g['stages_present'].split(';') != ['Stage '+str(i) for i in range(5)]:
        raise RuntimeError('The manuscript primary MMSE comparison requires all five stages.')
    if 'time' not in g['random_effects']:raise RuntimeError('Random-intercept fallback detected; the manuscript model has random slopes.')
    s=next(x for x in read_csv(out/'Fig5B_replication_summary.csv') if x['set']=='all analysed proteins')
    report={'status':'completed','protein_family':34,'ADNC_model_n':int(a['n']),
      'protein_concordant':int(s['n_same_direction']),'protein_replicated':int(s['n_replicated_q_lt_0.05_same_direction']),
      'MMSE_participants':int(g['n_participants']),'MMSE_assessments':int(g['n_assessments']),
      'MMSE_random_effects':g['random_effects'],'elapsed_seconds':round(time.perf_counter()-start,3)}
    (out/'run_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()

import json, numpy as np, sys
def stats(R, nominal=8604.8):
    errs=[]; Ms=[]; raw=[]; chis=[]; dt=[]
    for r in R:
        ev=r['events']; k=np.array([e['k'] for e in ev]); M=np.array([e['M'] for e in ev]); est=np.array([e['est'] for e in ev])
        err=est-M; raw+=list(err); ref=(k%4==0)
        p=np.polyfit(k[ref],err[ref],1); c=err-np.polyval(p,k)
        errs+=list(c[~ref]); Ms+=list(M[~ref]); chis+=[e['chi'] for e in ev]
        dt+=[e['te']-e['tau'] for e in ev if 'tau' in e]
    e=np.array(errs); M=np.array(Ms); raw=np.array(raw)
    mis=((M+e>nominal)!=(M>nominal)).mean()
    out=dict(n=int(e.size), raw_mean=float(raw.mean()), raw_sd=float(raw.std()), mean=float(e.mean()), sd=float(e.std()),
             robust_sd=float(1.4826*np.median(np.abs(e-np.median(e)))), mis=float(mis), frac_gt40=float(np.mean(np.abs(e)>40)))
    if dt: out['tension_err_sd']=float(np.std(dt))
    return out, e, M
if __name__=='__main__':
    for tag in sys.argv[1:]:
        R=json.load(open(f'tension_{tag}.json'))
        o,_,_=stats(R)
        print(tag, len(R), 'tubes', {k:(round(v,4) if isinstance(v,float) else v) for k,v in o.items()})

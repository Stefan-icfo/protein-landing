import json, numpy as np
R=json.load(open('batch_win05.json'))
ev=[dict(e,tube=r['seed']) for r in R for e in r['events']]
err=np.array([e['est']-e['M'] for e in ev]); M=np.array([e['M'] for e in ev]); est=np.array([e['est'] for e in ev])
k=np.array([e['k'] for e in ev]); tube=np.array([e['tube'] for e in ev])
print('slit window +-0.05 L, 20 tubes x 30 landings')
print(f' raw: mean {err.mean():+.1f}  SD {err.std():.1f}  robust SD {1.4826*np.median(np.abs(err-np.median(err))):.1f}  misassigned {((est>8604.8)!=(M>8604.8)).mean():.3f}')
# drift correction: per tube, linear trend vs landing index estimated from interleaved references.
# references = every 4th landing (index%4==0) treated as known-mass calibrants
corr=np.empty_like(err); used=np.zeros(len(err),bool)
for t in np.unique(tube):
    m=tube==t; kk=k[m]; ee=err[m]
    ref=(kk%4==0)
    p=np.polyfit(kk[ref], ee[ref], 1)
    corr[m]=ee-np.polyval(p,kk)
    used[np.where(m)[0][~ref]]=True       # evaluate only on non-reference landings
e2=corr[used]; est2=(M+corr)[used]; M2=M[used]
print(f' drift-corrected (1 in 4 landings as reference): n={e2.size} mean {e2.mean():+.1f} SD {e2.std():.1f} robust SD {1.4826*np.median(np.abs(e2-np.median(e2))):.1f} misassigned {((est2>8604.8)!=(M2>8604.8)).mean():.3f}  |err|>40: {np.mean(np.abs(e2)>40):.3f}')
# add protein-specific terms (assumed, Gaussian): hydrogen/adducts 7 Da, scale calibration 5 Da
rng=np.random.default_rng(0)
tot=e2+rng.normal(0,np.hypot(7,5),e2.size)
print(f' + protein-specific terms (7 Da, 5 Da): SD {tot.std():.1f}  misassigned {(((M2+tot)>8604.8)!=(M2>8604.8)).mean():.3f}')
# histogram data for the figure (drift-corrected, mechanical + protein terms)
bins=np.arange(8480,8740,10)
h_ub,_=np.histogram((M2+tot)[M2<8600],bins); h_p,_=np.histogram((M2+tot)[M2>8600],bins)
json.dump(dict(bins=bins.tolist(), ub=h_ub.tolist(), pub=h_p.tolist(), n_ub=int((M2<8600).sum()), n_pub=int((M2>8600).sum()),
               sd_mech=float(e2.std()), sd_tot=float(tot.std()), mis=float((((M2+tot)>8604.8)!=(M2>8604.8)).mean())), open('hist.json','w'))
print(' histogram counts Ub :', h_ub.tolist()); print(' histogram counts pUb:', h_p.tolist())
# whole-tube random landings (earlier batch): usable molecules per tube at <=15 Da
R2=json.load(open('batch_mildA.json'))+json.load(open('batch_mildB.json'))
ev2=[e for r in R2 for e in r['events']]
xe=np.array([e['xe'] for e in ev2]); kk=np.array([e['k'] for e in ev2])
print(f' random landings: central fraction {np.mean(np.abs(xe-0.5)<=0.05):.3f}; usable (central, first 5 landings) per tube: {np.sum((np.abs(xe-0.5)<=0.05)&(kk<5))/len(R2):.2f}')

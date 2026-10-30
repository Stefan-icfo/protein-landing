"""Case A (random tubes, slit +-0.05 L, tension steps 3e-5) at worse frequency stability."""
import json, sys, time
import batch_random as br
sigma_y = float(sys.argv[1]); tag = sys.argv[2]; seeds = [int(a) for a in sys.argv[3].split(',')]
res = []; t = time.time()
for sd in seeds:
    res.append(br.tube(sd, sigma_tau=3e-5, sigma_y=sigma_y))
    json.dump(res, open(f'noise_{tag}.json', 'w'))
    print(f'{tag} {sd} {time.time() - t:.0f}s', flush=True)

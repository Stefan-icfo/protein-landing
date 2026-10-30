"""Case A conditions with a chosen landing window and frequency stability (20 random tubes, seeds 600-619)."""
import json, sys, time
import batch_random as br
sigma_y = float(sys.argv[1]); window = float(sys.argv[2]); tag = sys.argv[3]
res = []; t = time.time()
for sd in range(600, 620):
    res.append(br.tube(sd, sigma_tau=3e-5, sigma_y=sigma_y, window=window))
    json.dump(res, open(f'wn_{tag}.json', 'w'))
    print(f'{tag} {sd} {time.time() - t:.0f}s', flush=True)

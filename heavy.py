"""Same protocol as batch.py for a heavier protein pair (RNase A 13,682 Da and +80 Da), 20 landings per tube (~30 % of the tube mass)."""
import json, sys, time
import batch
batch.UB, batch.PUB, batch.NOM = 13682.0, 13762.0, 13722.0
seeds = [int(a) for a in sys.argv[1].split(',')]
res = []; t = time.time()
for sd in seeds:
    res.append(batch.tube(sd, kind='mild', n_land=20, window=0.05))
    json.dump(res, open('batch_heavy.json', 'w'))
    print(f'heavy {sd} {time.time() - t:.0f}s', flush=True)

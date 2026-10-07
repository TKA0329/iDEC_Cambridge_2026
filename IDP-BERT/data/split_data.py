import os
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, PowerTransformer


def main(name='split0', oversample=True):
    with open('./data/data.pkl', 'rb') as f:
        raw = pickle.load(f)

    # scaler = StandardScaler()
    scaler = PowerTransformer()
    raw['rog'] = scaler.fit_transform(raw['rog'].reshape(-1, 1)).reshape(-1)
    raw['cv'] = scaler.fit_transform(raw['cv'].reshape(-1, 1)).reshape(-1)
    raw['tau'] = scaler.fit_transform(raw['tau'].reshape(-1, 1)).reshape(-1)

    data = [(
        raw['seqs'][i], raw['rog'][i],
        raw['cv'][i], raw['tau'][i]
    ) for i in range(len(raw['seqs']))]

    buckets = {i+116: [] for i in range(39, 1199, 116)}
    for d in data:
        for b in buckets:
            if len(d[0]) <= b:
                buckets[b].append(d)
                break

    train_buckets, val_buckets, test_buckets = {}, {}, {}
    tv_buckets = {}
    for b in buckets:
        tv_buckets[b], test_buckets[b] = train_test_split(buckets[b], test_size=0.2)
        train_buckets[b], val_buckets[b] = train_test_split(tv_buckets[b], test_size=0.25)

    if oversample:
        max_bucket_size = max([len(buckets[b]) for b in buckets.keys()])
        for b in train_buckets:
            multiplier = round(max_bucket_size/len(train_buckets[b]))
            train_buckets[b] = train_buckets[b]*multiplier

    data_tr = []
    for b in train_buckets:
        data_tr.extend(train_buckets[b])

    data_vl = []
    for b in val_buckets:
        data_vl.extend(val_buckets[b])

    data_te = []
    for b in test_buckets:
        data_te.extend(test_buckets[b])

    train_X, train_y = [d[0] for d in data_tr], [d[1:] for d in data_tr]
    val_X, val_y = [d[0] for d in data_vl], [d[1:] for d in data_vl]
    test_X, test_y = [d[0] for d in data_te], [d[1:] for d in data_te]

    os.makedirs(f'./data/{name}', exist_ok=True)

    np.savez(f'./data/{name}/train.npz', X=train_X, y=train_y)
    np.savez(f'./data/{name}/val.npz', X=val_X, y=val_y)
    np.savez(f'./data/{name}/test.npz', X=test_X, y=test_y)

main()

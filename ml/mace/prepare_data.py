#!/usr/bin/env python3
"""Prepare motif-disjoint graphene SAC training splits from dft.zip.
No dependencies beyond Python standard library.
"""
import argparse
import collections
import csv
import random
import re
import zipfile
from pathlib import Path


def metadata(header, key):
    m = re.search(rb'(?:^|\s)' + key.encode() + rb'=([^\s]+)', header)
    return m.group(1).decode('utf-8') if m else ''


def read_blocks(f):
    while True:
        line = f.readline()
        if not line:
            break
        if not line.strip():
            continue
        n = int(line.strip())
        head = f.readline()
        block = [line, head]
        for _ in range(n):
            v = f.readline()
            if not v:
                raise ValueError('Truncated extxyz frame')
            block.append(v)
        yield head, b''.join(block)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('zipfile', type=Path, help='Original dft.zip')
    ap.add_argument('--output', type=Path, default=Path('mlff_data'))
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    groups = collections.defaultdict(list)
    stats = collections.Counter()
    with zipfile.ZipFile(args.zipfile) as z:
        with z.open('graphene_relax_gbrv.extxyz') as f:
            for head, frame in read_blocks(f):
                motif = metadata(head, 'motif')
                state = metadata(head, 'state')
                run = metadata(head, 'run')
                if not motif or not run or not state:
                    raise ValueError(f'Missing motif/run/state in {head[:150]!r}')
                if b'energy=' not in head or b'forces:R:3' not in head:
                    raise ValueError(f'Missing training labels in {head[:150]!r}')
                groups[motif].append((run, state, frame))
                stats[state] += 1

    metals = collections.defaultdict(list)
    for motif in groups:
        metal = motif.split('-')[0]
        metals[metal].append(motif)

    rng = random.Random(args.seed)
    split_for_motif = {}
    for metal, motifs in sorted(metals.items()):
        motifs = sorted(motifs)
        rng.shuffle(motifs)
        # Rare elements stay in train; never claim OOD validation for them.
        if len(motifs) < 4:
            nval = ntest = 0
        else:
            nval = max(1, round(len(motifs) * 0.10))
            ntest = max(1, round(len(motifs) * 0.10))
        for i, motif in enumerate(motifs):
            split_for_motif[motif] = 'valid' if i < nval else ('test' if i < nval + ntest else 'train')

    files = {k: (args.output / f'{k}.xyz').open('wb') for k in ('train','valid','test')}
    manifest = []
    frames = collections.Counter()
    try:
        for motif in sorted(groups):
            split = split_for_motif[motif]
            for run, state, frame in groups[motif]:
                files[split].write(frame)
                frames[split] += 1
            manifest.append({
                'motif': motif,
                'metal': motif.split('-')[0],
                'split': split,
                'n_frames': len(groups[motif]),
                'n_runs': len(set(x[0] for x in groups[motif])),
                'states': ','.join(sorted(set(x[1] for x in groups[motif]))),
            })
    finally:
        for f in files.values():
            f.close()
    with (args.output / 'manifest.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['motif','metal','split','n_frames','n_runs','states'])
        writer.writeheader()
        writer.writerows(manifest)
    print('Frames by split:', dict(frames))
    print('Motifs by split:', dict(collections.Counter(v for v in split_for_motif.values())))
    print('States:', dict(stats))
    print('WARNING: W and other metals with <4 motifs are train-only; no held-out test for those metals.')
    print('Written to', args.output.resolve())

if __name__ == '__main__':
    main()

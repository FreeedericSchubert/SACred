#!/usr/bin/env python3
"""Quick audit of MACE train/validation/test extxyz, with special focus on 'ref'.
Usage: python audit_v2.py mlff_data_v2
No third-party dependencies.
"""
from collections import defaultdict, Counter
from pathlib import Path
import re
import sys


def meta(header, field):
    found = re.search(r'(?:^|\s)' + re.escape(field) + r'=([^\s]+)', header)
    return found.group(1).strip('"') if found else '?'


def frames(fname):
    with open(fname, 'r') as f:
        while True:
            first = f.readline()
            if not first:
                return
            if not first.strip():
                continue
            n = int(first.strip())
            header = f.readline().strip()
            species = Counter()
            for _ in range(n):
                row = f.readline()
                if not row:
                    raise ValueError(f'Truncated frame in {fname}')
                species[row.split()[0]] += 1
            yield n, header, species


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else 'mlff_data_v2')
    for name in ('train', 'valid', 'test'):
        path = root / (name + '.xyz')
        counts = Counter(); motifs=defaultdict(set); sizes=defaultdict(Counter); elems=defaultdict(Counter)
        samples=defaultdict(list); lack_labels=Counter()
        for n,header,sp in frames(path):
            state=meta(header, 'state'); motif=meta(header, 'motif')
            counts[state]+=1
            motifs[state].add(motif)
            sizes[state][n]+=1
            elems[state].update(sp)
            if 'energy=' not in header or 'forces:R:3' not in header:
                lack_labels[state] +=1
            if len(samples[state])<3:
                samples[state].append((motif,n,dict(sp)))
        print(f'\n{name.upper()} — {sum(counts.values())} frames')
        for state, count in sorted(counts.items()):
            print(f'  {state:16s} frames={count:5d} motifs={len(motifs[state]):3d} sizes={dict(sizes[state])} elements={dict(elems[state])}')
            if state=='ref':
                print('    ref examples:', *samples[state],sep='\n      ')
        if lack_labels: print('  WARNING missing labels:',dict(lack_labels))

if __name__ == '__main__':
    main()

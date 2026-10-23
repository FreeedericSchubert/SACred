#!/bin/bash
# Everything that needs the campus network, in one go: results from Curta, new DFT inputs to Curta and submitted.
#   bash tools/campus.sh          # pull DFT and xTB results, push and submit new run directories, show the queue
#   bash tools/campus.sh pull     # only pull (e.g. the last pull before a data freeze)
# New run directories (no relax.out on Curta yet) get 24 cores and 12 h, so they need no restart while the VPN is down.
set -e
cd "$(dirname "$0")/.."
ssh -o ConnectTimeout=10 -o BatchMode=yes curta true || { echo "Curta not reachable (campus network needed)"; exit 1; }
echo "== $(date +%H:%M) pull"
bash dft/sync.sh pull
bash dft/sync.sh xtb
[ "$1" = pull ] && exit 0
echo "== $(date +%H:%M) push and submit new run directories"
bash dft/sync.sh push
ssh curta 'cd ~/sac_dft && for d in $(bash submit.sh -n runs/graphene 2>/dev/null | grep "would submit" | cut -d: -f1); do
    [ -f $d/relax.out ] && continue
    printf "NTASKS=24\nMEM=24G\nHOURS=12\n" > $d/resources
    sed -i "s/max_seconds = .*/max_seconds = 41999/" $d/relax.in
    bash submit.sh $d
done; squeue -u $USER -h -o "%.9j %.8T" | cut -c1-14 | sort | uniq -c'

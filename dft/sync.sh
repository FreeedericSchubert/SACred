#!/bin/bash
# Copy run directories between this machine and Curta (~/sac_dft); tar over ssh, no rsync needed here.
#   bash dft/sync.sh push    # scripts and new inputs to Curta; files already there are kept (relax.in may
#                            # hold a restart geometry written by step.py)
#   bash dft/sync.sh pull    # pw.x outputs, the generated scf.in and meta.json back to dft/runs (git-ignored outputs)
#   bash dft/sync.sh xtb [m] # xTB results of ~/sac_dft/xtb back to data/screen (only files of the last m minutes)
set -e
cd "$(dirname "$0")"
REMOTE=sac_dft
case "$1" in
push)
    tar -C curta -cf - step.py submit.sh mos2_h.py mos2_fe.py | ssh curta "mkdir -p $REMOTE && tar -C $REMOTE -xf -"
    tar -C . -cf - --exclude='*.out' --exclude='*.out.*' runs | ssh curta "tar -C $REMOTE --skip-old-files -xf -" ;;
pull)
    ssh curta "cd $REMOTE && find runs \( -name '*.out' -o -name 'relax.out.*' -o -name scf.in -o -name meta.json \) -print0 | tar -cf - --null -T -" \
        | tar -C . -xf - ;;
xtb)
    ssh curta "cd $REMOTE/xtb && find data/screen -type f -mmin -${2:-1000000} -print0 | tar -cf - --null -T -" \
        | tar -C .. --keep-newer-files -xf - ;;   # a local retry (xtb_scan.py retry) is newer: kept
*)
    echo "usage: bash dft/sync.sh push|pull|xtb [minutes]"; exit 1 ;;
esac

#!/bin/bash
# Rebuild every result table, the database, the website and the 3D view from the data on this machine, in order.
#   bash dft/sync.sh pull && bash dft/sync.sh xtb    # on the campus network: DFT outputs and xTB results from Curta
#   bash tools/refresh.sh                            # then this (about 15 min; the database build is the slow step)
#   bash tools/refresh.sh export                     # and copy everything to the Nextcloud folder as well
# Each step reads what the steps before it wrote; a step that fails stops the script.
set -e
cd "$(dirname "$0")/.."
PY=.venv/bin/python
step() { echo "== $(date +%H:%M:%S) $*"; "$@"; }
step $PY dft/analyse.py > /dev/null
step $PY dft/orr.py > /dev/null
step $PY dft/conditions.py > /dev/null
step $PY dft/points.py > /dev/null
step $PY screen/xtb_scan.py table
step $PY screen/calibrate.py > /dev/null
step $PY screen/xtb_supports.py table > /dev/null
step $PY screen/xtb_reactions.py table > /dev/null
step $PY screen/xtb_stability.py table > /dev/null
step $PY screen/stability_descriptors.py > /dev/null
step $PY screen/geometry_mismatch.py > /dev/null
step $PY screen/xtb_free_complexes.py table > /dev/null
step $PY site/zoom_viewer.py
step $PY site/build.py
step $PY db/build.py > /dev/null
if [ "$1" = export ]; then
    step $PY db/export.py ~/Nextcloud/SAC_database
fi
echo "== $(date +%H:%M:%S) done; publish site/index.html, restart the dashboards, commit"

#!/bin/bash
# Submit the single-atom-site Quantum ESPRESSO runs on Curta (SLURM). Lives in ~/sac_dft next to runs/.
#   bash submit.sh runs/graphene                  # every run directory below that has no finished SCF
#   bash submit.sh runs/graphene/Ni-S4_DV         # one motif (all its states)
#   bash submit.sh -n runs/graphene               # dry run: only list what would be submitted
# One job per structure: Gamma-point relaxation (restarted from its last geometry when an earlier job ran
# out of time), then the 3x3x1 k-point SCF at the relaxed geometry that the energies come from.
# Cluster settings as in our earlier MoS2 runs: QE 7.1 foss (mpirun, its OpenMPI lacks SLURM PMI), OmniPath/
# PSM2 off, work files on node-local /dev/shm (Lustre writes hung pw.x), the same bad nodes excluded.
# 16 cores and 20 GB per job (seff: ~12.5 GB, plus ~3 GB of work files in /dev/shm during the SCF): no full
# nodes, and a 6 h limit so the backfill scheduler can start the jobs in gaps. Jobs already queued or
# running under the same name are left alone. A run directory may override these with a `resources` file
# (NTASKS=, MEM=, HOURS=; the MoS2 cells), and a `gamma_only` marker skips the SCF (energy from the relaxation).
set -e
cd "$(dirname "$0")"
ROOT=$PWD
QE=QuantumESPRESSO/7.1-foss-2022a
EXCLUDE=c061,c114,c116,c124
NTASKS0=16      # defaults; a `resources` file overrides them for its own directory only
MEM0=20G
HOURS0=6
dry=""
[ "$1" = -n ] && { dry=1; shift; }
for d in $(find "${@:-runs}" -name relax.in -printf "%h\n" | sort); do
    name=sac_$(echo "${d#runs/}" | tr / _)
    (cd "$d" && python3 "$ROOT/step.py" done) && { echo "$d: done"; continue; }
    id=$(squeue -h -u "$USER" -n "$name" -o %i | head -1)
    [ -n "$id" ] && { echo "$d: queued or running ($id)"; continue; }
    [ -n "$dry" ] && { echo "$d: would submit"; continue; }
    nk=1
    grep -q "nspin = 2" "$d/relax.in" && nk=2    # spin up/down as two pools
    NTASKS=$NTASKS0 MEM=$MEM0 HOURS=$HOURS0       # reset: an earlier directory's resources must not carry over
    [ -f "$d/resources" ] && eval "$(tr '\n' ' ' < "$d/resources")"
    id=$(sbatch --parsable --job-name="$name" --partition=main --qos=standard --exclude=$EXCLUDE \
        --nodes=1 --ntasks=$NTASKS --mem=$MEM --time=$HOURS:00:00 --chdir="$d" --output=slurm.out <<SCRIPT
#!/bin/bash -l
module load $QE
export OMP_NUM_THREADS=1 OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader OMPI_MCA_mtl=^psm2,ofi
export OMPI_MCA_osc=^ucx,rdma OMPI_MCA_io=romio321 HDF5_USE_FILE_LOCKING=FALSE
export ESPRESSO_TMPDIR=/dev/shm/\$USER/qe_\$SLURM_JOB_ID
mkdir -p \$ESPRESSO_TMPDIR; trap "rm -rf \$ESPRESSO_TMPDIR" EXIT
step="python3 $ROOT/step.py"
if ! \$step relaxed; then
    \$step restart
    mpirun pw.x -nk $nk -in relax.in > relax.out
fi
if [ ! -f gamma_only ]; then
    \$step relaxed && \$step scf && mpirun pw.x -nk $nk -in scf.in > scf.out
fi
SCRIPT
)
    echo "$d: $id"
done

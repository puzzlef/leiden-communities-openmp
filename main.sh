#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --job-name slurm
#SBATCH --output=slurm.out
# source scl_source enable gcc-toolset-11
# module load hpcx-2.7.0/hpcx-ompi
# module load openmpi/4.1.5
src="leiden-communities-openmp"
out="$HOME/Logs/$src$1.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
if [[ "$DOWNLOAD" != "0" ]]; then
  rm -rf $src
  git clone https://github.com/puzzlef/$src
  cd $src
fi

# Fixed config
# export OMP_STACKSIZE="4G"
: "${TYPE:=float}"
: "${MAX_THREADS:=64}"
: "${REPEAT_METHOD:=1}"
# Define macros (dont forget to add here)
DEFINES=(""
"-DTYPE=$TYPE"
"-DMAX_THREADS=$MAX_THREADS"
"-DREPEAT_METHOD=$REPEAT_METHOD"
)

# Compile
g++ ${DEFINES[*]} -std=c++17 -O3 -fopenmp main.cxx  # -Wall -Wno-unknown-pragmas

# Run on each graph
runEach() {
# stdbuf --output=L ./a.out ~/Data/web-Stanford.mtx   0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/indochina-2004.mtx  0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2002.mtx         0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/arabic-2005.mtx     0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2005.mtx         0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/webbase-2001.mtx    0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/it-2004.mtx         0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/sk-2005.mtx         0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-LiveJournal.mtx 1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Orkut.mtx       1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/asia_osm.mtx        1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/europe_osm.mtx      1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_A2a.mtx        1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_V1r.mtx        1 0 2>&1 | tee -a "$out"
}

# Run 5 times
for i in {1..5}; do
  runEach
done

# Signal completion
curl -X POST "https://maker.ifttt.com/trigger/puzzlef/with/key/${IFTTT_KEY}?value1=$src$1"

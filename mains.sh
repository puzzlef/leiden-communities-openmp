#!/usr/bin/env bash
src="leiden-communities-openmp"
ulimit -s unlimited

# Download program
if [[ "$DOWNLOAD" != "0" ]]; then
  rm -rf $src
  git clone https://github.com/puzzlef/$src
  cd $src
fi

# Don't need to download program again.
export DOWNLOAD="0"

# 1. Static Leiden
./main.sh

# Signal completion
curl -X POST "https://maker.ifttt.com/trigger/puzzlef/with/key/${IFTTT_KEY}?value1=$src"

#!/usr/bin/env bash
# The International Stroke Trial database (version 2), Sandercock et al. 2011, University of Edinburgh,
# https://doi.org/10.7488/ds/104, ODC-By licence. A copy ships in the RESIDE R package on GitHub.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/src
[ -d data/src/RESIDE ] || git clone -q --depth 1 https://github.com/hehta/RESIDE data/src/RESIDE
Rscript -e 'load("data/src/RESIDE/data/IST.rda"); write.csv(IST, "data/ist_raw.csv", row.names = FALSE); cat(nrow(IST), "patients written to data/ist_raw.csv\n")'

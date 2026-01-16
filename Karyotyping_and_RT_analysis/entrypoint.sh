#!/bin/bash
source /miniforge3/etc/profile.d/conda.sh
conda activate screpliseq
exec "$@"
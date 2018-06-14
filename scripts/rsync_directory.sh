#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Backup scratch folder to group directory
# ----------------------------------------------------------------------------------------

# Note that this ignores certain files to save room, including:
# - IMPUTE files with pt in their name
# - 1000 Genomes BAM files
# See ~/.cvsignore for full list

# Use -n for dry-run

cd $SCRATCH

rsync -avhC batwa-bakiga-exomes $GRP/batwa-bakiga-exomes_MIRROR

cd ..

exit

# To sync scratch directory with mirror (restoring things deleted on scratch):
# cd $SCRATCH
# rsync --dry-run -vhrlt --update $GRP/batwa-bakiga-exomes_MIRROR/batwa-bakiga-exomes .
# cd $SCRATCH/batwa-bakiga-exomes
# Create fake IMPUTE parts
# python scripts/create_fake_impute_parts.py
# Download 1000G BAM files and set their timestamp correctly
# sh scripts/restore_BAM_files.sh

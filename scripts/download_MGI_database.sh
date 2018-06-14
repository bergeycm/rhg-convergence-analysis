#!/bash/bin

# ----------------------------------------------------------------------------------------
# --- Download MGI (Mouse Genome Database, etc.) phenotype association data
# ----------------------------------------------------------------------------------------

# Download Vertebrate Homology:
# 4. Mouse/Human Orthology with Phenotype Annotations (tab-delimited)

# See here: http://www.informatics.jax.org/downloads/reports/index.html

mkdir -p data/MGI/

wget http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt \
    -O data/MGI/HMD_HumanPhenotype.rpt

exit

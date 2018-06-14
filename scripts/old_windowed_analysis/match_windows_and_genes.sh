
#!/bin/sh

module load bedtools/2.24.0

# ----------------------------------------------------------------------------------------
# --- Match windows to genes contained within
# ----------------------------------------------------------------------------------------

REFGENE=refGene/refGene.sort.simple.gtf
WINS=results/windows_100kb.bed 

# --- Simplify refGene file more to just gene info, not exons

sed -e "s/_exon[0-9]\+//" $REFGENE > ${REFGENE/.gtf/.genes.gtf}

# --- Sort genome file (Uses GNU's alphanumeric sort)
sort -k1,1V data/hg19.genome > data/hg19.sort.genome

# --- Sort windows file
sort -k1,1V -k2,2n $WINS > ${WINS/.bed/.sort.bed}

# --- Match genes to their windows

mapBed  -a ${WINS/.bed/.sort.bed} \
        -b ${REFGENE/.gtf/.genes.gtf} \
        -c 9 -o distinct \
        -g data/hg19.sort.genome \
        > ${WINS/.bed/.genes.bed}

#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Find windows that have insufficient exonic bases for our analyses
# ----------------------------------------------------------------------------------------

module load bedtools/2.24.0

# --- Make BED of windows of size 100kb

# MySQL blocked on cluster
#    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#        "select chrom, size from hg19.chromInfo" > data/hg19.genome

wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes -O data/hg19.genome

bedtools makewindows -g data/hg19.genome -w 100000 > results/windows_100kb.bed

# --- Find RefGene exon overlap with these windows

OUT_EXON_COV=results/windows_100kb.exon_cov.bed
echo -en "chr\twin_start\twin_end\t" > $OUT_EXON_COV
echo -e  "num_exons\tnum_exonic_bases\twin_length\tperc_covered" >> $OUT_EXON_COV

bedtools coverage -a results/windows_100kb.bed \
    -b refGene/refGene.sort.simple.gtf  >> $OUT_EXON_COV

# Find passing windows containing >= 500 exonic basis...
awk '{ if ($5 >= 500) { print $0 } }' $OUT_EXON_COV > ${OUT_EXON_COV/.bed/.pass.bed}

# ...and failing windows.
awk '{ if ($5 < 500) { print $0 } }' $OUT_EXON_COV > ${OUT_EXON_COV/.bed/.fail.bed}

#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download BAM files that are not restored by rsync
# ----------------------------------------------------------------------------------------

# Download BAM file and set modification date to just after individual's txt file

for pop in GBR LWK; do

    while read ind; do

        IND_TXT=data/1000genomes/BAMs/${ind}.txt
        IND_BAM=data/1000genomes/BAMs/${ind}.mapped.ILLUMINA.bwa.${pop}.exome.bam

        if [[ ! -f $IND_BAM ]] || [[ ! -f $IND_BAM.bai ]]; then
            sh scripts/download_1000G_BAMs.sh $ind
        else
            echo "Skipping download of BAM file for $ind as it already exists."
        fi

        touch -r $IND_TXT -d '+1 minute' $IND_BAM
        touch -r $IND_TXT -d '+1 minute' $IND_BAM.bai

    done < data/1000genomes/subset_${pop}.txt

done

exit

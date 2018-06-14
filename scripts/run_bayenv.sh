#!/bin/bash

# ========================================================================================
# --- Run Bayenv
# ========================================================================================

mkdir -p bayenv_run

cd bayenv_run

SNPSFILE=all.pops.merged.cleaned.bayenv.SNPfile
SNPSFILE_MINI=all.pops.merged.cleaned.thin.bayenv.SNPfile
NUMPOPS=4

cp ../data/$SNPSFILE .
cp ../data/$SNPSFILE_MINI .

# ----------------------------------------------------------------------------------------
# --- Estimate covariance matrix
# ----------------------------------------------------------------------------------------

MATRIXFULL=bayenv.matrix.out.full
~/bin/bayenv2 -i $SNPSFILE_MINI -p $NUMPOPS -k 100000 > $MATRIXFULL
MATRIXFILE=bayenv.matrix.out
tail -n 5 $MATRIXFULL | head -n4 > $MATRIXFILE

# ----------------------------------------------------------------------------------------
# --- Create environmental file
# ----------------------------------------------------------------------------------------

ENVIRONFILE=bayenv.environfile.txt

# Order is
# ANDA BAKIGA BATWA BR1
echo -e "-1.0\t1.0\t-1.0\t1.0" > $ENVIRONFILE

# To confirm order is correct:
# From VCF to Bayenv format conversion script's STDERR:
# chr1    9323910 rs6688832
# {'BATWA': '33,67', 'ANDA': '12,8', 'BAKIGA': '57,43', 'BR1': '16,4'}
# $ grep -n "rs6688832" ../data/all.pops.merged.cleaned.bayenv.SNPmap
# 1430:chr1       9323910 rs6688832
# 1430 * 2 = 2860
# $ head -n 2860 ../data/all.pops.merged.cleaned.bayenv.SNPfile | tail -n2
# 12      57      33      16
# 8       43      67      4

# ----------------------------------------------------------------------------------------
# --- Run Bayenv on single test SNPs
# ----------------------------------------------------------------------------------------

head -n2 $SNPSFILE | sed -e "s/\t$//" > test.SNPfile

~/bin/bayenv2 -i test.SNPfile -m $MATRIXFILE -e $ENVIRONFILE -p $NUMPOPS \
    -k 100000 -n 1 -t -r 1234

rm bf_environ.bayenv.environfile.txt

# ----------------------------------------------------------------------------------------
# --- Run Bayenv on all SNPs
# ----------------------------------------------------------------------------------------

LINE_COUNT=`wc -l $SNPSFILE | cut -f1 -d' '`

#   # Process 20,000 lines, or 10,000 SNPs at a time
#   HUNK_SIZE=20000

# Process 300 lines, or 150 SNPs at a time
HUNK_SIZE=300

rm ../results/bayenv_BF.txt

for start in `seq 1 $HUNK_SIZE $((LINE_COUNT + $HUNK_SIZE))`; do

    rm -r snp_batch_*

    end=$((start + $HUNK_SIZE))
    end=$((end - 1))
    SNPSFILE_SUBSET=${SNPSFILE}_${start}_${end}

    sed -n "${start},${end}p" $SNPSFILE > $SNPSFILE_SUBSET

    split -a 10 -l 2 $SNPSFILE_SUBSET snp_batch_

    # Can really be deleted
    ### rm $SNPSFILE_SUBSET

    module load parallel/20150122

    export ENVIRONFILE=$ENVIRONFILE
    export MATRIXFILE=$MATRIXFILE
    export NUMPOPS=$NUMPOPS

    function process_snp {

        SNP=$1

        ~/bin/bayenv2 -i $SNP -e $ENVIRONFILE -m $MATRIXFILE \
            -k 100000 -r $RANDOM -p $NUMPOPS -n 1 -t -o ${SNP}_OUTPUT

        rm $SNP
        rm $SNP.freqs

    }

    export -f process_snp

    ls snp_batch_* | \
        parallel --progress --jobs 30 process_snp

    # Back-up copy
    cat `ls -v *OUTPUT*` > bayenv_BF_${start}_${end}.txt

    cat `ls -v *OUTPUT*` >> ../results/bayenv_BF.txt

    rm -r snp_batch_*

done

cd ..
rm -r bayenv_run

# ----------------------------------------------------------------------------------------
# --- Combine results with original SNP info
# ----------------------------------------------------------------------------------------

paste data/${SNPSFILE/SNPfile/SNPmap} results/bayenv_BF.txt | \
    cut -f 1-3,5 > results/bayenv_BF_full.txt

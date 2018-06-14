#!/bin/sh

# ========================================================================================
# --- Impute missing SNPs by piping jobs to qsub
# ========================================================================================

IMPUTE2=$HOME/bin/impute_v2.3.2_x86_64_static/impute2
DIR_1000G=/gpfs/cyberstar/ghp3/Bergey/1000g_phase3/1000GP_Phase3

# ----------------------------------------------------------------------------------------
# --- Impute missing haplotypes with and without reference panel
# ----------------------------------------------------------------------------------------

for USE_REF in 1 0; do

    QSUB_OPTS="-l nodes=1:ppn=1,walltime=4:00:00 -m abe -M cxb585@psu.edu -N impute"
    if [ $USE_REF -eq 1 ]; then
        QSUB_OPTS="-l nodes=1:ppn=1,walltime=12:00:00,mem=20gb "
        QSUB_OPTS="$QSUB_OPTS -m abe -M cxb585@psu.edu -N impute_1kg"
    fi

    REF_SUFFIX=
    if [ $USE_REF -eq 1 ]; then
        REF_SUFFIX=.1000g
    fi

    for chr in `seq 1 22`; do

        # Find length of chromosome
        CHR_LEN=`grep -P "^chr${chr}\t" data/hg19.genome | cut -f 2`

        c=0

        for chunk_start in `seq 0 5000000 $CHR_LEN`; do

            chunk_end=$((chunk_start + 4999999))

            OUT_PREFIX=results/impute/AGRHUM_EASTERN_100x267251.1M
            OUT_PREFIX=${OUT_PREFIX}.${chr}${REF_SUFFIX}.gen.impute2_pt${c}

            CMD="$IMPUTE2 \
                -m ${DIR_1000G}/genetic_map_chr${chr}_combined_b37.txt \
                -g results/impute/AGRHUM_EASTERN_100x267251.1M.${chr}.gen \
                -int $chunk_start $chunk_end \
                -Ne 20000 \
                -o $OUT_PREFIX"

            # Add reference panel option if used
            if [ $USE_REF -eq 1 ]; then
                CMD="${CMD} \
                    -phase \
                    -h ${DIR_1000G}/1000GP_Phase3_chr${chr}.hap.gz \
                    -l ${DIR_1000G}/1000GP_Phase3_chr${chr}.legend.gz"
            fi

            # Run job via piping to qsub
            echo $CMD
            echo "cd $PWD; $CMD" | qsub $QSUB_OPTS

            sleep 2

            c=$((c + 1))
        done
    done
done

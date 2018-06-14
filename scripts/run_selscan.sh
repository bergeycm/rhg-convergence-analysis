#!/bin/sh

# ========================================================================================
# --- Compute iHS with selScan
# ========================================================================================

module load glib/2.48.2
module load gcc

SELSCAN=/storage/home/cxb585/bin/selscan/src/selscan
PREDICT_GMAP=/storage/home/cxb585/bin/predictGMAP/predictGMAP

DIR_1000G=/gpfs/cyberstar/ghp3/Bergey/1000g_phase3/1000GP_Phase3

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

chr=$1

for pop in kiga twa; do

    IMPUTE_DIR=results/impute/
    IN_HAP=${IMPUTE_DIR}/AGRHUM_EASTERN_100x267251.1M.$chr.1000g.gen.impute2_${pop}_haps

    MAP_1000G=$DIR_1000G/genetic_map_chr${chr}_combined_b37.txt

    OUT_FILE=results/selscan/`basename ${IN_HAP/.gen.*_haps/.${pop}.selscan}`

    # ------------------------------------------------------------------------------------
    # --- Reformat input data
    # ------------------------------------------------------------------------------------

    # Remove info columns from hap file
    cut -d' ' -f 6- $IN_HAP > $IN_HAP.fix

    # Reformat MAP file into format desired by selscan
    grep -v "^pos" $MAP_1000G | \
        awk  -v chr="$chr" 'BEGIN { OFS="\t" } { print chr,chr":"$1,$3,$1 }' > \
        $MAP_1000G.fix

    # ------------------------------------------------------------------------------------
    # --- Interpolate map info for SNPs missing from 1000G genetic map
    # ------------------------------------------------------------------------------------

    # Get list of SNPs in dataset to reduce map file
    cut -d' ' -f3 $IN_HAP > $IN_HAP.snplist_tmp

    # Interpolate to find missing SNP locations with predictGMAP
    # Max gap size is massive to include all SNPs
    $PREDICT_GMAP \
        --ref $MAP_1000G.fix \
        --query $IN_HAP.snplist_tmp \
        --out $MAP_1000G.predict.map \
        --max-gap 1000000

    sed -i_BACKUP -e "s/chr//" $MAP_1000G.predict.map

    # Check that all SNPs made it through
    PREDICTED_COUNT=`wc -l $MAP_1000G.predict.map | cut -f 1 -d' '`
    ORIG_COUNT=`wc -l $IN_HAP.snplist_tmp | cut -f 1 -d ' '`

    echo "Number of SNPs interpolated: [$PREDICTED_COUNT]"
    echo "Total number of SNPs: [$ORIG_COUNT]"

    if [ $PREDICTED_COUNT -ne $ORIG_COUNT ]; then

        NUM_MISSING=$(($ORIG_COUNT - $PREDICTED_COUNT))

        echo "Not all SNPs interpolated (N=$NUM_MISSING). Reducing haps file."

        # --- Figure out which SNPs did not get interpolated and are
        # --- missing from the map file and remove them from the haps file

        module load R

        Rscript scripts/remove_loci_in_haps.R \
            $IN_HAP.snplist_tmp \
            $MAP_1000G.predict.map \
            $IN_HAP.fix

        module unload R

        mv $IN_HAP.fix.reduced $IN_HAP.fix
        rm $IN_HAP.snplist_tmp

    fi

    # Remove MAP file now that interpolated one exists and temporary file for sed
    rm $MAP_1000G.fix
    rm $MAP_1000G.predict.map_BACKUP

    # ------------------------------------------------------------------------------------
    # --- Transpose haps input file
    # ------------------------------------------------------------------------------------

    # awk script to transpose file taken from:
    # http://stackoverflow.com/a/1729980
    awk '{
        for (i=1; i<=NF; i++)  {
            a[NR,i] = $i
        }
    }
    NF>p { p = NF }
    END {
        for(j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str" "a[i,j];
            }
            print str
        }
    }' $IN_HAP.fix > $IN_HAP.transpose

    # ------------------------------------------------------------------------------------
    # --- Run selscan
    # ------------------------------------------------------------------------------------

    $SELSCAN --ihs --hap $IN_HAP.transpose --map $MAP_1000G.predict.map --out $OUT_FILE

    # Delete haplotype file in selscan format
    rm $IN_HAP.fix

done

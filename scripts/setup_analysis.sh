
# --- Get scratch dir ready to go

# Download repo to /gpfs/scratch/cxb585/batwa_bakiga_exomes

# Make symlinks
ln -s /gpfs/cyberstar/ghp3/Bergey/batwa_bakiga_exomes/East_AGR_RHG.txt \
    data/East_AGR_RHG.txt
ln -s /gpfs/cyberstar/ghp3/Bergey/batwa_bakiga_exomes/AGRHUM_EASTERN_100/AGRHUM_EASTERN_100x267251.vcf.gz \
    data/AGRHUM_EASTERN_100x267251.vcf.gz

mkdir data/Perry_etal_2014-PNAS-1M/
for end in bed bim fam; do

	file_base=Batwa_Kiga.913651pos.230samples.PNAS2014.${end}

	orig_file=/gpfs/cyberstar/ghp3/Bergey/batwa_bakiga_exomes
	orig_file=${orig_file}/Perry_etal_2014-PNAS-1M/${file_base}

	link_file=data/Perry_etal_2014-PNAS-1M/${file_base}

	ln -s ${orig_file} ${link_file}

done

# Unpack cram files to bam format

qsub -t 1-100 pbs/cram2bam.pbs

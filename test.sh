snakemake -np --jobs 8 --cluster "qsub -V -M cxb585@psu.edu -m abe -l nodes=1:ppn={threads},walltime={params.runtime}:00:00{params.mem}"

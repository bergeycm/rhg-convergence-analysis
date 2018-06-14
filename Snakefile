
# snakemake --jobs 8 --cluster "qsub -V -M cxb585@psu.edu -m abe -l nodes=1:ppn={threads},walltime={params.runtime}:00:00{params.mem}"

# ----------------------------------------------------------------------------------------
# --- Variables
# ----------------------------------------------------------------------------------------

PERRY_2014_BEDS = [
    "GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc",
    "BayeScan_outlier_regions/Perry_etal_2014-S3-baka-vs-nzebi-bzime_BayeScan_outliers",
    "BayeScan_outlier_regions/Perry_etal_2014-S3-batwa-vs-bakiga_BayeScan_outliers",
    "iHS_outlier_regions/Perry_etal_2014-S3-baka_iHS_outliers",
    "iHS_outlier_regions/Perry_etal_2014-S3-bakiga_iHS_outliers",
    "iHS_outlier_regions/Perry_etal_2014-S3-batwa_iHS_outliers",
    "iHS_outlier_regions/Perry_etal_2014-S3-nzebi-nzime_iHS_outliers"
]

# ----------------------------------------------------------------------------------------

IDS_1000G, = glob_wildcards("data/1000genomes/BAMs/{ind_id}.txt")

inds_GBR = open("data/1000genomes/subset_GBR.txt", "r").readlines()
inds_GBR = [ind.rstrip('\n') for ind in inds_GBR]
inds_LWK = open("data/1000genomes/subset_LWK.txt", "r").readlines()
inds_LWK = [ind.rstrip('\n') for ind in inds_LWK]

# ----------------------------------------------------------------------------------------

GENE_LIST_ROOTS = [
    "eAGR_eRHG.weir.fst",
    "ihs.twa.txt",    "ihs.kiga.txt",
    "pbs.Batwa.GBR",  "pbs.Bakiga.GBR",
    "pbs.Anda.LWK",   "pbs.BR1.LWK",
    "bayenv_BF_pass",
]

GENE_LIST_ROOTS_PBS = list(filter(lambda x:'pbs' in x, GENE_LIST_ROOTS))

GENE_LIST_ROOTS_PBS_IHS = GENE_LIST_ROOTS_PBS[:]
GENE_LIST_ROOTS_PBS_IHS.append("ihs.twa.txt")
GENE_LIST_ROOTS_PBS_IHS.append("ihs.kiga.txt")

STATS = ['mean', 'max']
ONTS = ["BP", "MF", "CC"]
ONTS_BPMF = ["BP", "MF", "BPMF"]

# ----------------------------------------------------------------------------------------

import re

pattern = re.compile("^chr([0-9]+)$")

IMPUTE_OUTPUT = []
IMPUTE_OUTPUT_BY_CHR = {}

with open("data/hg19.genome") as gen:
    for line in gen:
        chr = line.split()[0]
        m = pattern.match(chr)
        if (m):
            chr_num = m.group(1)
            chr_len = int(line.split()[1])
            prefix = "results/impute/AGRHUM_EASTERN_100x267251.1M." + str(chr_num) + ".1000g.gen.impute2_pt"
            IMPUTE_OUTPUT_BY_CHR[chr_num] = []
            for i in range(0, 1 + chr_len // 5000000):
                IMPUTE_OUTPUT.append(prefix + str(i))
                IMPUTE_OUTPUT_BY_CHR[chr_num].append(prefix + str(i))

# ----------------------------------------------------------------------------------------

# Count of GO terms with pruning of terms with nodeSize=50
MAX_GO = {}
MAX_GO['BP'] = 5513
MAX_GO['MF'] = 1036
MAX_GO['CC'] =  389

GO_RANGES = []

for ont in ONTS:
    for start in range(1, MAX_GO[ont], 100):
        end = start + 99
        if (end > MAX_GO[ont]):
            end = MAX_GO[ont]
        GO_RANGES.append(ont + '_' + str(start) + '_' + str(end))

# ----------------------------------------------------------------------------------------
# --- Make all
# ----------------------------------------------------------------------------------------

rule all:
    input:
        # get_hg19_genome:
        "genomes/hg19/hg19.fa",
        "genomes/hg19/hg19.dict",
        "genomes/hg19/hg19.fa.fai",
        # get_refGene:
        "refGene/refGene.sort.simple.gtf",
        "refGene/refGene.sort.simple.justGenes.gtf",
        # make_sample_lists:
        "data/eRHG_IDs.txt",
        "data/eAGR_IDs.txt",
        # liftover_perry2014
        expand("data/Perry_etal_2014/{bed}.hg19.bed", bed=PERRY_2014_BEDS),
        # get_perry2014_genes:
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.genes.bed",
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.justGenes.txt",
        # download_mgi:
        "data/MGI/HMD_HumanPhenotype.rpt",
        # combine_growth_genes
        "data/all_growth_genes.txt",
        # download_andamanese:
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.vcf.gz",
        "data/Mondal_indiv_list.txt",
        "data/Mondal_indiv_list_Anda.txt",
        "data/Mondal_indiv_list_BR1.txt",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.SNP.locs.vcf",
        # annotate:
        "results/AGRHUM_EASTERN_100x267251.hg19_multianno.vcf",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.hg19_multianno.vcf",
        "reports/AGRHUM_EASTERN_100x267251_SNP_class_counts.txt",
        "reports/GreatApe_Indiafinal_Austosome.sm.recode_SNP_class_counts.txt",
        "reports/AGRHUM_EASTERN_100x267251_Polyphen_prediction_counts.txt",
        "reports/GreatApe_Indiafinal_Austosome.sm.recode_Polyphen_prediction_counts.txt",
        "results/AGRHUM_EASTERN_100x267251.syn.vcf.gz",
        "results/AGRHUM_EASTERN_100x267251.nonsyn.vcf.gz",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.syn.vcf.gz",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.nonsyn.vcf.gz",
        "results/AGRHUM_EASTERN_100x267251_stopgain.genes.txt",
        "results/GreatApe_Indiafinal_Austosome.sm.recode_stopgain.genes.txt",
        # find_pure_inds:
        "data/eAGR_IDs.iHS.txt",
        "data/eRHG_IDs.iHS.txt",
        "data/eAGR_IDs.iHS.VCF.txt",
        "data/eRHG_IDs.iHS.VCF.txt",
        # find_paralogs:
        expand("results/all_het_sites.e{pop}.txt", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe.filtered.bed", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe.filtered.justGenes.txt", pop=['AGR', 'RHG']),
        # fst:
        "results/eAGR_eRHG.weir.fst",
        "results/eAGR_eRHG.windowed.weir.fst",
        # merge_1M:
        "results/AGRHUM_EASTERN_100x267251.1M.bed",
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.gen", chr=range(1,23)),
        # run_impute:
        IMPUTE_OUTPUT,
        # combine_impute:
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2", chr=range(1,23)),
        # make_pop_haps:
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_kiga_haps", chr=range(1,23)),
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_twa_haps", chr=range(1,23)),
        # compute_ihs:
        expand("results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.twa.selscan.ihs.out", chr=range(1,23)),
        expand("results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.kiga.selscan.ihs.out", chr=range(1,23)),
        # combine_normalize_iHS:
        "results/ihs.twa.txt",
        "results/ihs.kiga.txt",
        # prep_for_1000G_SNP_calling:
        "genomes/hs37d5/hs37d5.fa",
        "genomes/hs37d5/hs37d5.dict",
        "genomes/hs37d5/hs37d5.fa.fai",
        "data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf",
        # download_1000G_indiv_info:
        "data/1000genomes/integrated_call_samples.20130502.seq.ped",
        # randomly_select_1000G_pops:
        "data/1000genomes/subset_GBR.txt",
        "data/1000genomes/subset_LWK.txt",
        # download_1000G_BAMs_GBR:
        expand("data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.GBR.exome.bam", id=inds_GBR),
        # download_1000G_BAMs_LWK:
        expand("data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.LWK.exome.bam", id=inds_LWK),
        # call_1000G_SNPs_GBR and call_1000G_SNPs_LWK:
        expand("data/1000genomes/hs37d5_snps/{pop}.chr{chr}.pass.snp.vcf", pop=['GBR', 'LWK'], chr=range(1,23)),
        # combine_1000G_SNPs_GBR and combine_1000G_SNPs_LWK:
        expand("data/1000genomes/hs37d5_snps/hs37d5.{pop}.pass.snp.vcf.gz", pop=['GBR', 'LWK']),
        # merge_1000G_batwa_bakiga_SNPs_GBR:
        "data/AGRHUM_EASTERN_100x267251.GBR.vcf",
        # merge_1000G_batwa_bakiga_SNPs_LWK:
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf",
        # run_admixture:
        expand("results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf.{k}.Q",
            k=range(2,7)),
        expand("results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.{k}.Q",
            k=range(2,7)),
        # plot_admixture:
        expand("results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf.{k}.pdf",
            k=range(2,7)),
        expand("results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.{k}.pdf",
            k=range(2,7)),
        # compute_Fst_for_PBS_GBR:
        "results/PBS_Bakiga_Batwa.weir.fst",
        "results/PBS_Bakiga_GBR.weir.fst",
        "results/PBS_Batwa_GBR.weir.fst",
        # compute_Fst_for_PBS_LWK:
        "results/PBS_Anda_BR1.weir.fst",
        "results/PBS_BR1_LWK.weir.fst",
        "results/PBS_Anda_LWK.weir.fst",
        # compute_PBS_GBR:
        "results/pbs.Batwa.GBR.bed",
        "results/pbs.Bakiga.GBR.bed",
        "reports/PBS_outliers.GBR.pdf",
        "reports/PBS_Batwa_vs_Fst.pdf",
        "reports/PBS_Bakiga_vs_Fst.pdf",
        "results/GBR.FST_for_PBS.txt",
        # compute_PBS_LWK:
        "results/pbs.Anda.LWK.bed",
        "results/pbs.BR1.LWK.bed",
        "reports/PBS_outliers.LWK.pdf",
        "reports/PBS_Anda_vs_Fst.pdf",
        "reports/PBS_BR1_vs_Fst.pdf",
        "results/LWK.FST_for_PBS.txt",
        # compute_pi:
        expand("results/pi_{pop}.sites.pi", pop=['Anda', 'Bakiga', 'Batwa', 'BR1',
            'GBR', 'GBR_all_pops', 'LWK', 'LWK_all_pops_BR1']),
        # merge_all_pops:
        "data/all.pops.merged.vcf.gz",
        "data/all.pops.merged.clean.recode.vcf",
        "data/all.pops.merged.cleaned.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.bayenv.SNPmap",
        "data/all.pops.merged.clean.thin.recode.vcf",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPmap",
        # run_bayenv:
        "results/bayenv_BF_full.txt",
        # filter_bayenv_results:
        "results/all.pops.merged.clean.overly.missing.uniq.txt",
        # plot_bayenv:
        "reports/bayenv_BFs.pdf",
        "results/bayenv_BF_pass.txt",
        "results/bayenv_BF_pass.bed",
        # plot_bayenv_geno:
        "results/bayenv_outlier_plotting_sentinel.txt",
        # make_gene_lists:
        expand("results/{root}.{stat}.txt", root=GENE_LIST_ROOTS, stat=STATS),
        expand("results/ihs.{pop}.txt.high.txt", pop=['twa', 'kiga']),
        expand("results/{root}.anno.bed", root=GENE_LIST_ROOTS),
        # split_stats_into_N_S:
        expand("results/AGRHUM_EASTERN_100x267251.{group}.bed",
            group=['syn', 'nonsyn']),
        expand("results/GreatApe_Indiafinal_Austosome.sm.recode.{group}.bed",
            group=['syn', 'nonsyn']),
        expand("results/{root}.anno.{group}.bed",
            root=GENE_LIST_ROOTS, group=['syn', 'nonsyn']),
        # find_outlier_snps:
        expand("reports/{root}.anno.paralog.flt.pdf",
            root=GENE_LIST_ROOTS_PBS_IHS),
        expand("reports/{root}.anno.extreme.vals{flt_suf}{all_suf}.txt",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt',''], all_suf=['.all','']),
        expand("reports/{root}.anno.extreme.vals{flt_suf}.tex",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt','']),
        expand("reports/{root}.anno.gwas.region.extreme.vals{flt_suf}.{end}",
            root=GENE_LIST_ROOTS_PBS_IHS,
            end=['txt','tex'], flt_suf=['.flt','']),
        expand("reports/{root}.anno.wood.region.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt','']),
        expand("reports/{root}.anno.MGI.region.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt','']),
        expand("reports/{root}.anno.lui.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt','']),
        expand("reports/{root}.anno.fisher.overrep{flt_suf}.{end}",
            root=GENE_LIST_ROOTS_PBS_IHS,
            end=['txt','tex'], flt_suf=['.flt','']),
        expand("reports/{root}.anno.wilcoxshift{flt_suf}.txt",
            root=GENE_LIST_ROOTS_PBS_IHS,
            flt_suf=['.flt','']),
        # compile_latex_outlier_snp:
        "reports/extreme_SNP_results_tables.pdf",
        # compute_joint_pval_snps:
        "results/joint_SNP_pvalues.txt",
        "reports/joint_SNP_pvalues.tex",
        "reports/joint_SNP_pvalues.deleterious.tex",
        # get_rsid_joint_pval_snps:
        "results/joint_SNP_pvalues.rsid.txt",
        # permute_to_get_pval:
        expand("results/{root}.pvals.txt",
            root=GENE_LIST_ROOTS),
        # permute_to_get_pval_size_corrected:
        expand("results_sizeCor/{root}.pvals.txt",
            root=GENE_LIST_ROOTS),
        # permute_to_get_pval_maf_corrected:
        expand("results_mafCor/{root}.pvals.txt",
            root=GENE_LIST_ROOTS),
        # permute_to_get_joint_pval:
        "results/pbs.Batwa.Anda.jointp.txt",
        # find_outlier_genes:
        expand("reports/{root}.pvals.paralog.flt.pdf",
            root=GENE_LIST_ROOTS),
        expand("reports/{root}.pvals.extreme.vals{flt_suf}{all_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt',''], all_suf=['.all','']),
        expand("reports/{root}.pvals.extreme.vals{flt_suf}.tex",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.gwas.region.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.wood.region.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.MGI.region.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.lui.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.growthgo.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.ophid_igf1.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.ophid_gh1.extreme.vals{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.fisher.overrep{flt_suf}.{end}", end=['txt','tex'],
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.pvals.wilcoxshift{flt_suf}.txt",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        # compute_joint_pval_genes:
        "results/joint_gene_pvalues.txt",
        "reports/joint_gene_pvalues.tex",
        # do_overrep:
        # - P-values - All SNPs
        expand("reports/{root}.pvals.{ont}.{testtype}.tex",
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - P-values - All SNPs (simple text version)
        expand("results/{root}.pvals.{ont}.{testtype}.results.txt",
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports/{root}.pvals.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("reports/{root}.pvals.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - Significant genes in significant GO terms - p-values - All SNPs
        expand("reports/{root}.pvals.sig.GO.sig.genes.{ont}.tex",
            root=GENE_LIST_ROOTS,
            ont=ONTS),
        # do_overrep_corrected:
        # - P-values - All SNPs
        expand("reports_{whatCor}Cor/{root}.pvals.{ont}.{testtype}.tex",
            whatCor=['size', 'maf'],
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - P-values - All SNPs (simple text version)
        expand("results_{whatCor}Cor/{root}.pvals.{ont}.{testtype}.results.txt",
            whatCor=['size', 'maf'],
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports_{whatCor}Cor/{root}.pvals.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            whatCor=['size', 'maf'],
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("reports_{whatCor}Cor/{root}.pvals.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            whatCor=['size', 'maf'],
            root=GENE_LIST_ROOTS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - Significant genes in significant GO terms - p-values - All SNPs
        expand("reports_{whatCor}Cor/{root}.pvals.sig.GO.sig.genes.{ont}.tex",
            whatCor=['size', 'maf'],
            root=GENE_LIST_ROOTS,
            ont=ONTS),
        # do_overrep_jointp:
        # - P-values - All SNPs
        expand("reports/pbs.{pops}.jointp.{ont}.{testtype}.tex",
            pops=['Batwa.Anda'],
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - Joint p-values - All SNPs (simple text version)
        expand("results/pbs.{pops}.jointp.{ont}.{testtype}.results.txt",
            pops=['Batwa.Anda'],
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # - P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports/pbs.{pops}.jointp.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            pops=['Batwa.Anda'],
            ont=ONTS, testtype=['overrep', 'enrich']),
        expand("reports/pbs.{pops}.jointp.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            pops=['Batwa.Anda'],
            ont=ONTS, testtype=['overrep', 'enrich']),
        # - Significant genes in significant GO terms - p-values - All SNPs
        expand("reports/pbs.{pops}.jointp.sig.GO.sig.genes.{ont}.tex",
            pops=['Batwa.Anda'],
            ont=ONTS),
        # compute_empirical_joint_p:
        expand("results/{root}.pvals.{ont}.{testtype}.results.wGOiterp.txt",
            root=GENE_LIST_ROOTS_PBS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("results/convergent_permuted_pval_{pairs}.{ont}.{testtype}.txt",
            pairs=["Batwa.GBR-Anda.LWK", "Bakiga.GBR-BR1.LWK"],
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # compute_empirical_joint_p_cor:
        expand("results_{cortype}/{root}.pvals.{ont}.{testtype}.results.wGOiterp.txt",
            cortype=["sizeCor", "mafCor"],
            root=GENE_LIST_ROOTS_PBS,
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("results_{cortype}/convergent_permuted_pval_{pairs}.{ont}.{testtype}.txt",
            cortype=["sizeCor", "mafCor"],
            pairs=["Batwa.GBR-Anda.LWK", "Bakiga.GBR-BR1.LWK"],
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # compile_latex:
        expand("reports/{root}.{stat}.{ont}.pdf",
            root=GENE_LIST_ROOTS,
            stat=STATS, ont=ONTS),
        expand("reports/ihs.{pop}.txt.high.{ont}.pdf",
            ont=ONTS, pop=['twa', 'kiga']),
        "reports/joint_SNP_pvalues.pdf",
        "reports/joint_GO_pvalues.pdf",
        # compute_joint_pval_GO:
        expand("results/joint_GO_pvalues.{testtype}.{ending}",
            testtype=['enrich', 'overrep'],
            ending=['txt', 'tex']),
        # download_ancestral_info:
        "data/Allele.bcp",
        "data/SNPAncestralAllele.bcp",
        # compute_global_maf_afr:
        "results/AGRHUM_EASTERN_100x267251.frq",
        # compute_global_maf_asn:
        "results/GreatApe_Indiafinal_Austosome.sm.recode.frq",
        # compute_pop_freq:
        expand("results/AGRHUM_EASTERN_100x267251.{pop}.frq", pop=['RHG', 'AGR']),
        # compute_plot_daf:
        "results/ref.info.txt",
        "reports/DAF.batwa.density.pdf",
        "reports/DAF.bakiga.density.pdf",
        "reports/DAF.PBS.batwa.pdf",
        "reports/DAF.PBS.bakiga.pdf",
        "reports/DAF.batwa.bakiga.pdf",
        # combine_bp_mf_results:
        expand("results{cortype}/pbs.{pop}.pvals.BPMF.{testtype}.results.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pop=['Batwa.GBR', 'Bakiga.GBR', 'Anda.LWK', 'BR1.LWK'],
            testtype=['overrep', 'enrich']),
        expand("results{cortype}/convergent_permuted_pval_{pair}.BPMF.{testtype}.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pair=['Batwa.GBR-Anda.LWK', 'Bakiga.GBR-BR1.LWK'],
            testtype=['overrep', 'enrich']),
        # plot_go_results:
        expand("reports/{root}.pvals.{ont}.GO.overrep.pdf",
            root=["pbs.Batwa.GBR","pbs.Bakiga.GBR","pbs.Anda.LWK","pbs.BR1.LWK",
                  "bayenv_BF_pass"],
            ont=ONTS),
        # plot_go_results_obs_vs_exp:
        expand("reports/pbs.Batwa-Anda.pvals.{ont}.GO.overrep.pdf",
            ont=ONTS_BPMF),
        expand("reports/pbs.Bakiga-BR1.pvals.{ont}.GO.overrep.pdf",
            ont=ONTS_BPMF),
        # plot_go_results_obs_vs_exp_corrected:
        expand("reports{cortype}/pbs.Batwa-Anda.pvals.{ont}.GO.overrep.pdf",
            cortype=["_sizeCor", "_mafCor"], ont=ONTS_BPMF),
        expand("reports{cortype}/pbs.Bakiga-BR1.pvals.{ont}.GO.overrep.pdf",
            cortype=["_sizeCor", "_mafCor"], ont=ONTS_BPMF),
        # plot_go_results_enrich:
        expand("reports/pbs.Batwa-Anda.pvals.{ont}.GO.enrich.pdf",
            ont=ONTS_BPMF),
        # plot_go_results_enrich:
        expand("reports{cortype}/pbs.Batwa-Anda.pvals.{ont}.GO.enrich.pdf",
            cortype=["_sizeCor", "_mafCor"], ont=ONTS_BPMF),
        # plot_pbs_by_snp_along_genome:
        "reports/PBS_by_SNP.pdf",
        "results/gene_SNP_counts.txt",
        # plot_pbs_by_gene_along_genome:
        "reports/PBS_by_gene.pdf",
        "reports/PBS_by_gene_length.pdf",
        "reports/PBS_by_SNP_count.pdf",
        "reports/PBS_stddev_by_SNP_count.pdf",
        "reports/PBS_stderr_by_SNP_count.pdf",
        "reports/PBS_boxplots_by_SNP_count.pdf",
        # plot_convergent_pbs_by_snp_along_genome:
        "reports/convergent_PBS_p_by_SNP.pdf",
        # plot_bf_by_snp_along_genome:
        "reports/BF_by_SNP.pdf",
        # plot_bf_by_gene_along_genome:
        "reports/BF_by_gene.pdf",
        # annotate_pbs_for_schematic:
        "results/GBR.FST_for_PBS.txt.anno.bed",
        "results/LWK.FST_for_PBS.txt.anno.bed",
        # plot_pbs_schematic:
        "reports/pbs_schematic.pdf",
        # plot_map:
        "reports/sampling_map.pdf",
        # parse_apriori_gene_list_results
        "reports/parsed_apriori_gene_list_results.txt",

localrules: all, make_sample_lists, combine_impute, find_pure_inds, randomly_select_1000G_pops

# ----------------------------------------------------------------------------------------
# --- Download and make sequence dictionary for hg19 genome
# ----------------------------------------------------------------------------------------

rule get_hg19_genome:
    output:
        "genomes/hg19/hg19.fa",
        "genomes/hg19/hg19.dict",
        "genomes/hg19/hg19.fa.fai",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_hg19.sh"

# ----------------------------------------------------------------------------------------
# --- Download and simplify refGene
# ----------------------------------------------------------------------------------------

rule get_refGene:
    output:
        "refGene/refGene.sort.simple.gtf",
        "refGene/refGene.sort.simple.justGenes.gtf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_simplify_refGene.sh"

# ----------------------------------------------------------------------------------------
# --- Generate sample lists
# ----------------------------------------------------------------------------------------

rule make_sample_lists:
    input:
        "data/East_AGR_RHG.txt"
    output:
        "data/eRHG_IDs.txt",
        "data/eAGR_IDs.txt"
    shell:
        "sh scripts/make_pop_sample_lists.sh"

# ----------------------------------------------------------------------------------------
# --- LiftOver regions of interest from Perry et al 2014
# ----------------------------------------------------------------------------------------

rule liftover_perry2014:
    input:
        expand("data/Perry_etal_2014/{bed}.bed", bed=PERRY_2014_BEDS),
        "/storage/home/cxb585/liftover_chain_files/hg18ToHg19.over.chain"
    output:
        expand("data/Perry_etal_2014/{bed}.hg19.bed", bed=PERRY_2014_BEDS)
    threads: 1
    params: runtime="1"
    shell:
        "sh scripts/liftover_ROIs.sh"

# ----------------------------------------------------------------------------------------
# --- Make list of genes in Perry et al 2014 phenotype associated regions
# ----------------------------------------------------------------------------------------

rule get_perry2014_genes:
    input:
        "refGene/refGene.sort.simple.justGenes.gtf",
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.bed"
    output:
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.genes.bed",
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.justGenes.txt"
    threads: 1
    params: runtime="1"
    shell:
        "sh scripts/get_Perry2014_genes.sh"

# ----------------------------------------------------------------------------------------
# --- Download MGI phenotype data
# ----------------------------------------------------------------------------------------

rule download_mgi:
    output:
        "data/MGI/HMD_HumanPhenotype.rpt"
    threads: 1
    params: runtime="1"
    shell:
        "sh scripts/download_MGI_database.sh"

# ----------------------------------------------------------------------------------------
# --- Make combined list of all a priori growth genes
# ----------------------------------------------------------------------------------------

rule combine_growth_genes:
    input:
        "data/Wood_etal_2014/OMIM_height_genes.txt",
        "data/MGI/HMD_HumanPhenotype.rpt",
        "data/Lui_etal_2012/growth_plate_expressed.txt",
        "data/GO_growth_genes/GO_0040007_genes.txt",
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.justGenes.txt",
    output:
        "data/all_growth_genes.txt",
    threads: 1
    params: runtime="1"
    shell:
        "sh scripts/make_curated_growth_gene_list.sh"

# ----------------------------------------------------------------------------------------
# --- Download Andamanese data
# ----------------------------------------------------------------------------------------

rule download_andamanese:
    output:
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.vcf.gz",
        "data/Mondal_indiv_list.txt",
        "data/Mondal_indiv_list_BR1.txt",
        "data/Mondal_indiv_list_Anda.txt",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.SNP.locs.vcf"
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/download_andamanese.sh"

# ----------------------------------------------------------------------------------------
# --- Annotate SNPs in VCF
# ----------------------------------------------------------------------------------------

rule annotate:
    input:
        "data/AGRHUM_EASTERN_100x267251.vcf.gz",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz",
    output:
        "results/AGRHUM_EASTERN_100x267251.hg19_multianno.vcf",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.hg19_multianno.vcf",
        "reports/AGRHUM_EASTERN_100x267251_SNP_class_counts.txt",
        "reports/GreatApe_Indiafinal_Austosome.sm.recode_SNP_class_counts.txt",
        "reports/AGRHUM_EASTERN_100x267251_Polyphen_prediction_counts.txt",
        "reports/GreatApe_Indiafinal_Austosome.sm.recode_Polyphen_prediction_counts.txt",
        "results/AGRHUM_EASTERN_100x267251.syn.vcf.gz",
        "results/AGRHUM_EASTERN_100x267251.nonsyn.vcf.gz",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.syn.vcf.gz",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.nonsyn.vcf.gz",
        "results/AGRHUM_EASTERN_100x267251_stopgain.genes.txt",
        "results/GreatApe_Indiafinal_Austosome.sm.recode_stopgain.genes.txt"
    threads: 8
    params: runtime="24"
    shell:
        "sh scripts/annotate_SNPs.sh"

# ----------------------------------------------------------------------------------------
# --- Figure out which individuals to include based on ancestry inferred by
# --- Perry et al 2014
# ----------------------------------------------------------------------------------------

rule find_pure_inds:
    input:
        "data/Perry_etal_2014/Perry_etal_2014-individual_info-height_map.csv",
        "data/East_AGR_RHG.txt"
    output:
        "data/eAGR_IDs.iHS.txt",
        "data/eRHG_IDs.iHS.txt",
        "data/eAGR_IDs.iHS.VCF.txt",
        "data/eRHG_IDs.iHS.VCF.txt",
    shell:
        "module load r/3.4; "
        "Rscript scripts/make_pop_lists_for_iHS.R"

# ----------------------------------------------------------------------------------------
# --- Find putative paralogs based on heterozygosity
# ----------------------------------------------------------------------------------------

rule find_paralogs:
    input:
        "results/AGRHUM_EASTERN_100x267251.hg19_multianno.vcf",
        "data/eAGR_IDs.iHS.VCF.txt",
        "data/eRHG_IDs.iHS.VCF.txt",
    output:
        expand("results/all_het_sites.e{pop}.txt", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe.filtered.bed", pop=['AGR', 'RHG']),
        expand("results/HWE_info.e{pop}.hwe.filtered.justGenes.txt", pop=['AGR', 'RHG']),
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/find_putative_paralogs.sh"

# ----------------------------------------------------------------------------------------
# --- Compute FST on a per-SNP and windowed basis
# ----------------------------------------------------------------------------------------

# Also can be run with PBS script:
# qsub -t 1-3 pbs/compute_Fst.pbs

rule fst:
    input:
        vcf="data/AGRHUM_EASTERN_100x267251.vcf.gz",
        popAGR="data/eAGR_IDs.txt",
        popRHG="data/eRHG_IDs.txt"
    output:
        "results/eAGR_eRHG.weir.fst",
        "results/eAGR_eRHG.windowed.weir.fst"
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/compute_Fst.sh {input.vcf}"

# ----------------------------------------------------------------------------------------
# --- Merge in 1M SNP chip data from Perry et al 2014 PNAS paper
# ----------------------------------------------------------------------------------------

rule merge_1M:
    input:
        "data/AGRHUM_EASTERN_100x267251.vcf.gz",
        "data/Perry_etal_2014-PNAS-1M/Batwa_Kiga.913651pos.230samples.PNAS2014.bed",
    output:
        "results/AGRHUM_EASTERN_100x267251.1M.bed",
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.gen", chr=range(1,23)),
        expand("results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.sample", chr=range(1,23))
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/combine_exome_and_snp_chip.sh"

# ----------------------------------------------------------------------------------------
# --- Phase with IMPUTE2 using 1000 Genomes reference panel (done in hunks in parallel)
# ----------------------------------------------------------------------------------------

# Also can be run with script that submits jobs:
# sh scripts/submit_impute_jobs.sh

# If IMPUTE fails because there are no SNPs, touch output file.

IMPUTE2 = "/storage/home/cxb585/bin/impute_v2.3.2_x86_64_static/impute2"
DIR_1000G = "/gpfs/cyberstar/ghp3/Bergey/1000g_phase3/1000GP_Phase3"

rule run_impute:
    input:
        "data/hg19.genome",
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.gen"
    output:
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_pt{part}"
    threads: 1
    params: runtime="12",
            mem=",mem=20gb"
    shell:
        "chunk_start=$(({wildcards.part}*5000000));"
        "chunk_end=$(($chunk_start+4999999));"
        "{IMPUTE2} -m {DIR_1000G}/genetic_map_chr{wildcards.chr}_combined_b37.txt "
        "-g results/impute/AGRHUM_EASTERN_100x267251.1M.{wildcards.chr}.gen "
        "-int $chunk_start $chunk_end "
        "-Ne 20000 "
        "-o {output} "
        "-phase "
        "-h {DIR_1000G}/1000GP_Phase3_chr{wildcards.chr}.hap.gz "
        "-l {DIR_1000G}/1000GP_Phase3_chr{wildcards.chr}.legend.gz;"
        "sh scripts/touch_missing_IMPUTE_output.sh {wildcards.chr} {wildcards.part}"

# ----------------------------------------------------------------------------------------
# --- Combine IMPUTE output (done in hunks in parallel)
# ----------------------------------------------------------------------------------------

rule combine_impute:
    input:
        lambda wildcards: IMPUTE_OUTPUT_BY_CHR[wildcards.chr]
    output:
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2",
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_haps"
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/combine_impute_results.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Make population-specific haps files for input into iHS
# ----------------------------------------------------------------------------------------

rule make_pop_haps:
    input:
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_haps",
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.sample",
        "data/eAGR_IDs.iHS.txt",
        "data/eRHG_IDs.iHS.txt"
    output:
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_kiga_haps",
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_twa_haps"
    threads: 1
    params: runtime="4"
    shell:
        "sh scripts/make_pop_haps_files.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Compute iHS
# ----------------------------------------------------------------------------------------

# Also can be run with PBS script:
# qsub -t 1-22 pbs/run_selscan.pbs

rule compute_ihs:
    input:
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_twa_haps",
        "results/impute/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.gen.impute2_kiga_haps",
        "/gpfs/cyberstar/ghp3/Bergey/1000g_phase3/1000GP_Phase3/genetic_map_chr{chr}_combined_b37.txt"
    output:
        "results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.twa.selscan.ihs.out",
        "results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.kiga.selscan.ihs.out"
    threads: 1
    params: runtime="20",
            mem=",mem=5gb"
    shell:
        "sh scripts/run_selscan.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Combine and normalize iHS
# ----------------------------------------------------------------------------------------

rule combine_normalize_iHS:
    input:
        expand("results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.twa.selscan.ihs.out", chr=range(1,23)),
        expand("results/selscan/AGRHUM_EASTERN_100x267251.1M.{chr}.1000g.kiga.selscan.ihs.out", chr=range(1,23))
    output:
        "results/ihs.twa.txt",
        "results/ihs.kiga.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/combine_normalize_iHS_output.sh"

# ----------------------------------------------------------------------------------------
# --- Prepare for 1000G SNP calling
# ----------------------------------------------------------------------------------------

rule prep_for_1000G_SNP_calling:
    input:
        "data/AGRHUM_EASTERN_100x267251.vcf.gz"
    output:
        "genomes/hs37d5/hs37d5.fa",
        "genomes/hs37d5/hs37d5.dict",
        "genomes/hs37d5/hs37d5.fa.fai",
        "data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/prepare_for_1000G_SNP_calling.sh"

# ----------------------------------------------------------------------------------------
# --- Download 1000G files to figure out individuals to use
# ----------------------------------------------------------------------------------------

rule download_1000G_indiv_info:
    output:
        "data/1000genomes/integrated_call_samples.20130502.seq.ped"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_1000G_individual_info.sh"

# ----------------------------------------------------------------------------------------
# --- Select random 1000G individuals for outside population in PBS
# ----------------------------------------------------------------------------------------

rule randomly_select_1000G_pops:
    input:
        "data/1000genomes/integrated_call_samples.20130502.seq.ped"
    output:
        "data/1000genomes/subset_GBR.txt",
        "data/1000genomes/subset_LWK.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/randomly_select_1000G_pops.sh"

# ----------------------------------------------------------------------------------------
# --- Download BAMs for 1000G populations
# ----------------------------------------------------------------------------------------

rule download_1000G_BAMs_GBR:
    input:
        expand("data/1000genomes/BAMs/{id}.txt", id=inds_GBR),
    output:
        "data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.GBR.exome.bam",
    threads: 1
    params: runtime="10",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_1000G_BAMs.sh {wildcards.id}"

rule download_1000G_BAMs_LWK:
    input:
        expand("data/1000genomes/BAMs/{id}.txt", id=inds_LWK),
    output:
        "data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.LWK.exome.bam"
    threads: 1
    params: runtime="10",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_1000G_BAMs.sh {wildcards.id}"

# ----------------------------------------------------------------------------------------
# --- Call and filter 1000G alleles at Batwa/Bakiga SNP locations
# ----------------------------------------------------------------------------------------

rule call_1000G_SNPs_GBR:
    input:
        "genomes/hs37d5/hs37d5.fa",
        expand("data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.GBR.exome.bam", id=inds_GBR),
        "data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf"
    output:
        "data/1000genomes/hs37d5_snps/GBR.chr{chr}.pass.snp.vcf"
    threads: 8
    params: runtime="12",
            mem=",mem=8gb"
    shell:
        "sh scripts/call_and_filter_1000G_SNPs.sh {wildcards.chr} GBR"

rule call_1000G_SNPs_LWK:
    input:
        "genomes/hs37d5/hs37d5.fa",
        expand("data/1000genomes/BAMs/{id}.mapped.ILLUMINA.bwa.LWK.exome.bam", id=inds_LWK),
        "data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf"
    output:
        "data/1000genomes/hs37d5_snps/LWK.chr{chr}.pass.snp.vcf"
    threads: 8
    params: runtime="12",
            mem=",mem=8gb"
    shell:
        "sh scripts/call_and_filter_1000G_SNPs.sh {wildcards.chr} LWK"

# ----------------------------------------------------------------------------------------
# --- Combine the 1000G alleles we have newly called
# ----------------------------------------------------------------------------------------

rule combine_1000G_SNPs_GBR:
    input:
        expand("data/1000genomes/hs37d5_snps/GBR.chr{chr}.pass.snp.vcf", chr=range(1,23))
    output:
        "data/1000genomes/hs37d5_snps/hs37d5.GBR.pass.snp.vcf.gz"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/cat_1000_genomes_SNPs.sh GBR"

rule combine_1000G_SNPs_LWK:
    input:
        expand("data/1000genomes/hs37d5_snps/LWK.chr{chr}.pass.snp.vcf", chr=range(1,23))
    output:
        "data/1000genomes/hs37d5_snps/hs37d5.LWK.pass.snp.vcf.gz"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/cat_1000_genomes_SNPs.sh LWK"

# ----------------------------------------------------------------------------------------
# --- Merge 1000G alleles and Batwa/Bakiga SNPs
# ----------------------------------------------------------------------------------------

rule merge_1000G_batwa_bakiga_SNPs_GBR:
    input:
        "data/1000genomes/hs37d5_snps/hs37d5.{pop}.pass.snp.vcf.gz",
        "data/AGRHUM_EASTERN_100x267251.vcf.gz"
    output:
        "data/AGRHUM_EASTERN_100x267251.{pop}.vcf"
    threads: 1
    params: runtime="12",
            mem=",mem=5gb"
    shell:
        "sh scripts/merge_1000G_ingroup_SNPs.sh {wildcards.pop}"

rule merge_1000G_batwa_bakiga_SNPs_LWK:
    input:
        "data/1000genomes/hs37d5_snps/hs37d5.{pop}.pass.snp.vcf.gz",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz"
    output:
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.{pop}.vcf"
    threads: 1
    params: runtime="12",
            mem=",mem=5gb"
    shell:
        "sh scripts/merge_1000G_ingroup_SNPs.sh {wildcards.pop}"

# ----------------------------------------------------------------------------------------
# --- Run ADMIXTURE
# ----------------------------------------------------------------------------------------

rule run_admixture:
    input:
        "data/AGRHUM_EASTERN_100x267251.GBR.vcf",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf"
    output:
        expand("results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf.{k}.Q",
            k=range(2,7)),
        expand("results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.{k}.Q",
            k=range(2,7))
    threads: 8
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/run_admixture.sh"

# ----------------------------------------------------------------------------------------
# --- Plot ADMIXTURE results
# ----------------------------------------------------------------------------------------

rule plot_admixture:
    input:
        expand("results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf.{k}.Q",
            k=range(2,7)),
        expand("results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.{k}.Q",
            k=range(2,7))
    output:
        expand("results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf.{k}.pdf",
            k=range(2,7)),
        expand("results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.{k}.pdf",
            k=range(2,7))
    threads: 8
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_admixture.R results/admixture/AGRHUM_EASTERN_100x267251.GBR.vcf. data/AGRHUM_EASTERN_100x267251.GBR.vcf.fam;"
        "Rscript scripts/plot_admixture.R results/admixture/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf. data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf.fam"

# ----------------------------------------------------------------------------------------
# --- Compute Fst for PBS
# ----------------------------------------------------------------------------------------

rule compute_Fst_for_PBS_GBR:
    input:
        "data/AGRHUM_EASTERN_100x267251.GBR.vcf",
        "data/eAGR_IDs.txt",
        "data/eRHG_IDs.txt",
        "data/1000genomes/subset_GBR.txt"
    output:
        "results/PBS_Bakiga_Batwa.weir.fst",
        "results/PBS_Bakiga_GBR.weir.fst",
        "results/PBS_Batwa_GBR.weir.fst"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_Fst_for_PBS.sh GBR"

rule compute_Fst_for_PBS_LWK:
    input:
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf",
        "data/Mondal_indiv_list_BR1.txt",
        "data/Mondal_indiv_list_Anda.txt",
        "data/1000genomes/subset_LWK.txt"
    output:
        "results/PBS_Anda_BR1.weir.fst",
        "results/PBS_BR1_LWK.weir.fst",
        "results/PBS_Anda_LWK.weir.fst"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_Fst_for_PBS.sh LWK"

# ----------------------------------------------------------------------------------------
# --- Compute PBS
# ----------------------------------------------------------------------------------------

rule compute_PBS_GBR:
    input:
        "results/PBS_Bakiga_Batwa.weir.fst",
        "results/PBS_Bakiga_GBR.weir.fst",
        "results/PBS_Batwa_GBR.weir.fst",
        "refGene/refGene.sort.simple.justGenes.gtf"
    output:
        "results/pbs.Batwa.GBR.bed",
        "results/pbs.Bakiga.GBR.bed",
        "reports/PBS_outliers.GBR.pdf",
        "reports/PBS_Batwa_vs_Fst.pdf",
        "reports/PBS_Bakiga_vs_Fst.pdf",
        "results/GBR.FST_for_PBS.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compute_pbs.R GBR"

rule compute_PBS_LWK:
    input:
        "results/PBS_Anda_BR1.weir.fst",
        "results/PBS_BR1_LWK.weir.fst",
        "results/PBS_Anda_LWK.weir.fst",
        "refGene/refGene.sort.simple.justGenes.gtf"
    output:
        "results/pbs.Anda.LWK.bed",
        "results/pbs.BR1.LWK.bed",
        "reports/PBS_outliers.LWK.pdf",
        "reports/PBS_Anda_vs_Fst.pdf",
        "reports/PBS_BR1_vs_Fst.pdf",
        "results/LWK.FST_for_PBS.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compute_pbs.R LWK"

# ----------------------------------------------------------------------------------------
# --- Compute within-population site diversity (pi)
# ----------------------------------------------------------------------------------------

rule compute_pi:
    input:
        "data/AGRHUM_EASTERN_100x267251.GBR.vcf",
        "data/eAGR_IDs.txt",
        "data/eRHG_IDs.txt",
        "data/1000genomes/subset_GBR.txt",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf",
        "data/Mondal_indiv_list_BR1.txt",
        "data/Mondal_indiv_list_Anda.txt",
        "data/1000genomes/subset_LWK.txt",
    output:
        expand("results/pi_{pop}.sites.pi", pop=['Anda', 'Bakiga', 'Batwa', 'BR1',
            'GBR', 'GBR_all_pops', 'LWK', 'LWK_all_pops_BR1'])
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "sh scripts/compute_pi.sh"

# ========================================================================================
# --- Run Bayenv
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Merge Batwa/Bakiga and Anda/BR1 VCFs and make Bayenv input file
# ----------------------------------------------------------------------------------------

rule merge_all_pops:
    input:
        "data/AGRHUM_EASTERN_100x267251.vcf.gz",
        "data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz"
    output:
        "data/all.pops.merged.vcf.gz",
        "data/all.pops.merged.clean.recode.vcf",
        "data/all.pops.merged.cleaned.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.bayenv.SNPmap",
        "data/all.pops.merged.clean.thin.recode.vcf",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPmap",
        "data/all.pop.list.txt"
    threads: 1
    params: runtime="12",
            mem=",mem=24gb"
    shell:
        "sh scripts/merge_all_pops.sh"

# ----------------------------------------------------------------------------------------
# --- Run Bayenv
# ----------------------------------------------------------------------------------------

rule run_bayenv:
    input:
        "data/all.pops.merged.cleaned.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.bayenv.SNPmap",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPfile",
        "data/all.pops.merged.cleaned.thin.bayenv.SNPmap"
    output:
        "results/bayenv_BF_full.txt"
    threads: 1
    params: runtime="12",
            mem=",mem=24gb"
    shell:
        "sh scripts/run_bayenv.sh"

# ----------------------------------------------------------------------------------------
# --- Filter Bayenv results for missingness
# ----------------------------------------------------------------------------------------

rule filter_bayenv_results:
    input:
        "results/bayenv_BF_full.txt"
    output:
        "results/all.pops.merged.clean.overly.missing.uniq.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=4gb"
    shell:
        "sh scripts/find_overly_missing_bayenv_snps.sh"

# ----------------------------------------------------------------------------------------
# --- Plot Bayenv and write BED
# ----------------------------------------------------------------------------------------

rule plot_bayenv:
    input:
        "results/bayenv_BF_full.txt",
        "results/all.pops.merged.clean.overly.missing.uniq.txt"
    output:
        "reports/bayenv_BFs.pdf",
        "results/bayenv_BF_pass.txt",
        "results/bayenv_BF_pass.bed",
    threads: 1
    params: runtime="1",
            mem=",mem=4gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_bayenv_results.R"

# ----------------------------------------------------------------------------------------
# --- Plot genotype frequencies by population for Bayenv outliers
# ----------------------------------------------------------------------------------------

rule plot_bayenv_geno:
    input:
        "results/bayenv_BF_pass.txt",
        "data/all.pops.merged.clean.recode.vcf",
        "data/all.pop.list.txt"
    output:
        "results/bayenv_outlier_plotting_sentinel.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=4gb"
    shell:
        "sh scripts/plot_bayenv_outlier_genotype_freq.sh"

# ========================================================================================
# --- Prepare to do overrep test
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Generate gene lists for overrep tests AND make annotated BED file along the way
# ----------------------------------------------------------------------------------------

rule make_gene_lists_fst:
    input:
        "results/eAGR_eRHG.weir.fst",
    output:
        expand("results/eAGR_eRHG.weir.fst.{stat}.txt", stat=STATS),
        expand("results/eAGR_eRHG.weir.fst.anno.bed")
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_gene_lists.pl {input}"

rule make_gene_lists_ihs:
    input:
        "results/ihs.{pop}.txt",
    output:
        expand("results/ihs.{{pop}}.txt.{stat}.txt", stat=STATS),
        "results/ihs.{pop}.txt.high.txt",
        "results/ihs.{pop}.txt.anno.bed",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_gene_lists.pl {input}"

rule make_gene_lists_pbs:
    input:
        "results/pbs.{pops}.bed",
    output:
        expand("results/pbs.{{pops}}.{stat}.txt", stat=STATS),
        "results/pbs.{pops}.anno.bed"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_gene_lists.pl {input}"

rule make_gene_lists_bayenv:
    input:
        "results/bayenv_BF_pass.bed",
    output:
        expand("results/bayenv_BF_pass.{stat}.txt", stat=STATS),
        "results/bayenv_BF_pass.anno.bed"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_gene_lists.pl {input}"

# ----------------------------------------------------------------------------------------
# --- Split stats into N and S groups (using annotated files from make_gene_lists step)
# ----------------------------------------------------------------------------------------

rule split_stats_into_N_S:
    input:
        expand("results/AGRHUM_EASTERN_100x267251.{group}.vcf.gz",
            group=['syn', 'nonsyn']),
        expand("results/GreatApe_Indiafinal_Austosome.sm.recode.{group}.vcf.gz",
            group=['syn', 'nonsyn']),
        expand("results/{root}.anno.bed", root=GENE_LIST_ROOTS)
    output:
        expand("results/AGRHUM_EASTERN_100x267251.{group}.bed",
            group=['syn', 'nonsyn']),
        expand("results/GreatApe_Indiafinal_Austosome.sm.recode.{group}.bed",
            group=['syn', 'nonsyn']),
        expand("results/{root}.anno.{group}.bed",
            root=GENE_LIST_ROOTS, group=['syn', 'nonsyn'])
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/split_stats_into_N_and_S.sh"

# ========================================================================================
# --- Find extreme SNPs in single populations and among two populations
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Find extreme SNPs and flag those that overlap GWAS regions or have predicted impact
# ----------------------------------------------------------------------------------------

rule find_outlier_snps:
    input:
        "results/{root}.anno.bed",
        "data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.bed",
        "data/Wood_etal_2014/OMIM_height_genes.txt",
        "data/MGI/HMD_HumanPhenotype.rpt",
        "results/AGRHUM_EASTERN_100x267251.hg19_multianno.vcf"
    output:
        "reports/{root}.anno.paralog.flt.pdf",
        expand("reports/{{root}}.anno.extreme.vals{flt_suf}{all_suf}.txt", flt_suf=['.flt',''], all_suf=['.all','']),
        expand("reports/{{root}}.anno.extreme.vals{flt_suf}.tex", flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.gwas.region.extreme.vals{flt_suf}.{end}", end=['txt','tex'], flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.wood.region.extreme.vals{flt_suf}.txt", flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.MGI.region.extreme.vals{flt_suf}.txt", flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.lui.extreme.vals{flt_suf}.txt", flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.growthgo.extreme.vals{flt_suf}.txt", flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.fisher.overrep{flt_suf}.{end}", end=['txt','tex'], flt_suf=['.flt','']),
        expand("reports/{{root}}.anno.wilcoxshift{flt_suf}.txt", flt_suf=['.flt','']),
    threads: 8
    params: runtime="8",
            mem=",mem=24gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/find_outlier_SNPs.R results/{wildcards.root}.anno.bed TRUE;"
        "Rscript scripts/find_outlier_SNPs.R results/{wildcards.root}.anno.bed FALSE"

# ----------------------------------------------------------------------------------------
# --- Compile LaTeX for outlier SNP tables
# ----------------------------------------------------------------------------------------

rule compile_latex_outlier_snp:
    input:
        expand("reports/{root}.anno.extreme.vals{flt_suf}.tex",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.anno.gwas.region.extreme.vals{flt_suf}.tex",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.anno.deleterious.extreme.vals{flt_suf}.tex",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
        expand("reports/{root}.anno.fisher.overrep{flt_suf}.tex",
            root=GENE_LIST_ROOTS, flt_suf=['.flt','']),
    output:
        "reports/extreme_SNP_results_tables.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/make_tables_SNP_results.sh"

# ----------------------------------------------------------------------------------------
# --- Compute joint p-value for SNP-based statistics
# ----------------------------------------------------------------------------------------

rule compute_joint_pval_snps:
    input:
        expand("reports/{root}.anno.extreme.vals.all.txt",
            root=GENE_LIST_ROOTS_PBS_IHS),
    output:
        "results/joint_SNP_pvalues.txt",
        "reports/joint_SNP_pvalues.tex",
        "reports/joint_SNP_pvalues.deleterious.tex"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compute_joint_pval_SNPs.R"

# ----------------------------------------------------------------------------------------
# --- Annotate outlier convergent SNPs with their rsids
# ----------------------------------------------------------------------------------------

rule get_rsid_joint_pval_snps:
    input:
        "results/joint_SNP_pvalues.txt",
    output:
        "results/joint_SNP_pvalues.rsid.txt",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/add_rsid_joint_pval_SNPs.sh > {output}"

# ----------------------------------------------------------------------------------------
# --- Permute SNPs and genes to get p-value for each gene
# ----------------------------------------------------------------------------------------

rule permute_to_get_pval:
    input:
        "results/{root}.anno.bed"
    output:
        "results/{root}.pvals.txt"
    threads: 20
    params: runtime="24",
            mem=",mem=24gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/permute_snps_to_get_gene_pval.R {input}"

rule permute_to_get_pval_size_corrected:
    input:
        "results/{root}.anno.bed"
    output:
        "results_sizeCor/{root}.pvals.txt"
    threads: 20
    params: runtime="24",
            mem=",mem=24gb"
    shell:
        "module load r/3.4; "
        "mkdir -p results_sizeCor; "
        "Rscript scripts/permute_snps_to_get_gene_pval_matched_by_size.R {input}"

rule permute_to_get_pval_maf_corrected:
    input:
        "results/{root}.anno.bed",
        "results/AGRHUM_EASTERN_100x267251.frq",
        "results/GreatApe_Indiafinal_Austosome.sm.recode.frq"
    output:
        "results_mafCor/{root}.pvals.txt"
    threads: 20
    params: runtime="24",
            mem=",mem=24gb"
    shell:
        "module load r/3.4; "
        "mkdir -p results_mafCor; "
        "Rscript scripts/permute_snps_to_get_gene_pval_matched_by_MAF.R {input}"

# ----------------------------------------------------------------------------------------
# --- Permute genes to get p-value for each gene
# ----------------------------------------------------------------------------------------

rule permute_to_get_joint_pval:
    input:
        batwa="results/pbs.Batwa.GBR.anno.bed",
        anda="results/pbs.Anda.LWK.anno.bed",
    output:
        ba_res="results/pbs.Batwa.Anda.jointp.txt"
    threads: 20
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/permute_snps_to_get_gene_joint_pval.R {input.batwa} {input.anda} {output.ba_res}.tmp;"
        "cut -f 1,4 {output.ba_res}.tmp > {output.ba_res}"

# ----------------------------------------------------------------------------------------
# --- Find outlier genes and flag those that overlap GWAS regions or have predicted impact
# ----------------------------------------------------------------------------------------

rule find_outlier_genes:
    input:
        "results/{root}.pvals.txt",
        "data/Wood_etal_2014/OMIM_height_genes.txt",
        "data/MGI/HMD_HumanPhenotype.rpt",
        expand("results/HWE_info.e{pop}.hwe.filtered.justGenes.txt", pop=['AGR', 'RHG']),
        "data/OPHID_GH1.txt",
        "data/OPHID_IGF1.txt"
    output:
        "reports/{root}.pvals.paralog.flt.pdf",
        expand("reports/{{root}}.pvals.extreme.vals{flt_suf}{all_suf}.txt",
            flt_suf=['.flt',''], all_suf=['.all','']),
        expand("reports/{{root}}.pvals.extreme.vals{flt_suf}.tex",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.gwas.region.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.wood.region.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.MGI.region.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.lui.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.growthgo.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.ophid_igf1.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.ophid_gh1.extreme.vals{flt_suf}.txt",
            flt_suf=['.flt','']),
        expand("reports/{{root}}.pvals.fisher.overrep{flt_suf}.{end}",
            flt_suf=['.flt',''], end=['txt','tex']),
        expand("reports/{{root}}.pvals.wilcoxshift{flt_suf}.txt",
            flt_suf=['.flt','']),
    threads: 8
    params: runtime="1",
            mem=",mem=24gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/find_outlier_genes.R results/{wildcards.root}.pvals.txt FALSE;"
        "Rscript scripts/find_outlier_genes.R results/{wildcards.root}.pvals.txt TRUE"

# ----------------------------------------------------------------------------------------
# --- Compute joint p-value for gene-based statistics
# ----------------------------------------------------------------------------------------

rule compute_joint_pval_genes:
    input:
        expand("reports/{root}.pvals.extreme.vals.all.txt", root=GENE_LIST_ROOTS),
    output:
        "results/joint_gene_pvalues.txt",
        "reports/joint_gene_pvalues.tex",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compute_joint_pval_genes.R"

# ========================================================================================
# --- Do overrep tests on all single pop stats
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Do overrep tests
# ----------------------------------------------------------------------------------------

rule do_overrep:
    input:
        expand("results/{{root}}.{stat}.txt",
            stat=STATS),
        "results/{root}.pvals.txt"
    output:
        # P-values - All SNPs
        expand("reports/{{root}}.pvals.{ont}.{testtype}.tex",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # P-values - All SNPs (simple text version)
        expand("results/{{root}}.pvals.{ont}.{testtype}.results.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports/{{root}}.pvals.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("reports/{{root}}.pvals.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # Significant genes in significant GO terms - p-values - All SNPs
        expand("reports/{{root}}.pvals.sig.GO.sig.genes.{ont}.tex",
            ont=ONTS),
    threads: 20
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/run_overrep_enrich_tests.sh results/{wildcards.root}.pvals.txt"

# Could be combined with above
rule do_overrep_corrected:
    input:
        "results_{whatCor}Cor/{root}.pvals.txt"
    output:
        # P-values - All SNPs
        expand("reports_{{whatCor}}Cor/{{root}}.pvals.{ont}.{testtype}.tex",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # P-values - All SNPs (simple text version)
        expand("results_{{whatCor}}Cor/{{root}}.pvals.{ont}.{testtype}.results.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports_{{whatCor}}Cor/{{root}}.pvals.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("reports_{{whatCor}}Cor/{{root}}.pvals.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # Significant genes in significant GO terms - p-values - All SNPs
        expand("reports_{{whatCor}}Cor/{{root}}.pvals.sig.GO.sig.genes.{ont}.tex",
            ont=ONTS),
    threads: 20
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/run_overrep_enrich_tests.sh results_{wildcards.whatCor}Cor/{wildcards.root}.pvals.txt"

rule do_overrep_jointp:
    input:
        "results/pbs.{pops}.jointp.txt"
    output:
        # P-values - All SNPs
        expand("reports/pbs.{{pops}}.jointp.{ont}.{testtype}.tex",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # Joint p-values - All SNPs (simple text version)
        expand("results/pbs.{{pops}}.jointp.{ont}.{testtype}.results.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports/pbs.{{pops}}.jointp.{ont}.bootstrapped.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        expand("reports/pbs.{{pops}}.jointp.{ont}.bootstrapped-geneGO.GO.p-values.{testtype}.txt",
            ont=ONTS,
            testtype=['overrep', 'enrich']),
        # Significant genes in significant GO terms - p-values - All SNPs
        expand("reports/pbs.{{pops}}.jointp.sig.GO.sig.genes.{ont}.tex",
            ont=ONTS),
    threads: 20
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/run_overrep_enrich_tests.sh results/pbs.{wildcards.pops}.jointp.txt"

# ----------------------------------------------------------------------------------------
# --- Compute empirical joint p-value for two populations (e.g. Batwa and Anda)
# ----------------------------------------------------------------------------------------

rule compute_empirical_joint_p:
    input:
        # P-values - All SNPs (simple text version)
        expand("results/{root}.pvals.{{ont}}.{{testtype}}.results.txt",
            root=GENE_LIST_ROOTS_PBS),
        # P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports/{root}.pvals.{{ont}}.bootstrapped.GO.p-values.{{testtype}}.txt",
            root=GENE_LIST_ROOTS_PBS),
        expand("reports/{root}.pvals.{{ont}}.bootstrapped-geneGO.GO.p-values.{{testtype}}.txt",
            root=GENE_LIST_ROOTS_PBS),
    output:
        expand("results/{root}.pvals.{{ont}}.{{testtype}}.results.wGOiterp.txt",
            root=GENE_LIST_ROOTS_PBS),
        "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.{testtype}.txt",
        "results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.{testtype}.txt"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compare_pvals_to_joint_distributions.R Batwa.GBR Anda.LWK {wildcards.ont} {wildcards.testtype} "";"
        "Rscript scripts/compare_pvals_to_joint_distributions.R Bakiga.GBR BR1.LWK {wildcards.ont} {wildcards.testtype} "";"

rule compute_empirical_joint_p_cor:
    input:
        # P-values - All SNPs (simple text version)
        expand("results_{{cortype}}/{root}.pvals.{{ont}}.{{testtype}}.results.txt",
            root=GENE_LIST_ROOTS_PBS),
        # P-values - All SNPs (bootstrapped GO p-value distributions)
        expand("reports_{{cortype}}/{root}.pvals.{{ont}}.bootstrapped.GO.p-values.{{testtype}}.txt",
            root=GENE_LIST_ROOTS_PBS),
        expand("reports_{{cortype}}/{root}.pvals.{{ont}}.bootstrapped-geneGO.GO.p-values.{{testtype}}.txt",
            root=GENE_LIST_ROOTS_PBS),
    output:
        expand("results_{{cortype}}/{root}.pvals.{{ont}}.{{testtype}}.results.wGOiterp.txt",
            root=GENE_LIST_ROOTS_PBS),
        "results_{cortype}/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.{testtype}.txt",
        "results_{cortype}/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.{testtype}.txt"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/compare_pvals_to_joint_distributions.R Batwa.GBR Anda.LWK {wildcards.ont} {wildcards.testtype} _{wildcards.cortype};"
        "Rscript scripts/compare_pvals_to_joint_distributions.R Bakiga.GBR BR1.LWK {wildcards.ont} {wildcards.testtype} _{wildcards.cortype};"

# ----------------------------------------------------------------------------------------
# --- Compile LaTeX tables
# ----------------------------------------------------------------------------------------

rule compile_latex:
    input:
        expand("reports/{root}.{stat}.{ont}.tex",
            root=GENE_LIST_ROOTS, stat=STATS, ont=ONTS),
        expand("reports/ihs.{pop}.txt.high.{ont}.tex",
            ont=ONTS, pop=['twa', 'kiga']),
        "reports/joint_SNP_pvalues.tex",
        "reports/joint_SNP_pvalues.deleterious.tex",
        expand("results/joint_GO_pvalues.{testtype}.tex",
            testtype=['enrich', 'overrep'])
    output:
        expand("reports/{root}.{stat}.{ont}.pdf",
            root=GENE_LIST_ROOTS, stat=STATS, ont=ONTS),
        expand("reports/ihs.{pop}.txt.high.{ont}.pdf",
            ont=ONTS, pop=['twa', 'kiga']),
        "reports/joint_SNP_pvalues.pdf",
        "reports/joint_SNP_pvalues.deleterious.pdf",
        expand("results/joint_GO_pvalues.{testtype}.pdf",
            testtype=['enrich', 'overrep'])
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/compile_latex.sh"

# ========================================================================================
# --- Compute joint p-value for GO terms
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute joint p-value for PBS and N vs. S stat by GO term
# ----------------------------------------------------------------------------------------

rule compute_joint_pval_GO:
    input:
        expand("results{corType}/{root}.pvals.{ont}.{testtype}.results.txt",
            corType=["", "_sizeCor", "_mafCor"],
            root=["pbs.Batwa.GBR","pbs.Bakiga.GBR","pbs.Anda.LWK","pbs.BR1.LWK",],
            ont=ONTS,
            testtype=['enrich', 'overrep']),
    output:
        expand("results/joint_GO_pvalues.{testtype}.{ending}",
            testtype=['enrich', 'overrep'],
            ending=['txt', 'tex'])
    threads: 1
    params: runtime="4",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/compute_joint_pval_GO_terms.R"

# ========================================================================================
# --- Explore derived allele frequency as a sanity check
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Download ancestral allele info
# ----------------------------------------------------------------------------------------

rule download_ancestral_info:
    output:
        "data/Allele.bcp",
        "data/SNPAncestralAllele.bcp"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_ancestral_allele_info.sh"

# ----------------------------------------------------------------------------------------
# --- Compute "global" allele frequencies
# ----------------------------------------------------------------------------------------

rule compute_global_maf_afr:
    input:
        vcf="data/AGRHUM_EASTERN_100x267251.vcf",
    output:
        "results/AGRHUM_EASTERN_100x267251.frq",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcftools --vcf {input.vcf} --freq --out results/AGRHUM_EASTERN_100x267251"

rule compute_global_maf_asn:
    input:
        vcf="data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf",
    output:
        "results/GreatApe_Indiafinal_Austosome.sm.recode.frq",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcftools --vcf {input.vcf} --freq --out results/GreatApe_Indiafinal_Austosome.sm.recode;"
        "mv results/GreatApe_Indiafinal_Austosome.sm.recode.frq results/GreatApe_Indiafinal_Austosome.sm.recode.error.frq;"
        "awk '{{ if ($3 == 2) print $0 }}' results/GreatApe_Indiafinal_Austosome.sm.recode.error.frq > results/GreatApe_Indiafinal_Austosome.sm.recode.frq"

# ----------------------------------------------------------------------------------------
# --- Compute population-specific allele frequencies
# ----------------------------------------------------------------------------------------

rule compute_pop_freq:
    input:
        vcf="data/AGRHUM_EASTERN_100x267251.vcf",
        pop_list="data/e{pop}_IDs.txt"
    output:
        "results/AGRHUM_EASTERN_100x267251.{pop}.frq",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcftools --vcf {input.vcf} --keep {input.pop_list} --freq --out results/AGRHUM_EASTERN_100x267251.{wildcards.pop}"

# ----------------------------------------------------------------------------------------
# --- Compute and plot DAF
# ----------------------------------------------------------------------------------------

rule compute_plot_daf:
    input:
        "data/Allele.bcp",
        "data/SNPAncestralAllele.bcp",
        expand("results/AGRHUM_EASTERN_100x267251.{pop}.frq", pop=['RHG', 'AGR']),
        "data/AGRHUM_EASTERN_100x267251.vcf",
        expand("results/pbs.Batwa.GBR.anno.{type}.bed", type=['syn', 'nonsyn'])
    output:
        "results/ref.info.txt",
        "reports/DAF.batwa.density.pdf",
        "reports/DAF.bakiga.density.pdf",
        "reports/DAF.PBS.batwa.pdf",
        "reports/DAF.PBS.bakiga.pdf",
        "reports/DAF.batwa.bakiga.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=24gb"
    shell:
        "module load r/3.4; "
        "Rscript scripts/explore_and_plot_DAF.R"

# ========================================================================================
# --- Plot results
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Combine BP and MF results
# ----------------------------------------------------------------------------------------

rule combine_bp_mf_results:
    input:
        expand("results{cortype}/pbs.{pop}.pvals.{ont}.{testtype}.results.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pop=['Batwa.GBR', 'Bakiga.GBR', 'Anda.LWK', 'BR1.LWK'],
            ont=['BP', 'MF'],
            testtype=['overrep', 'enrich']),
        expand("results{cortype}/convergent_permuted_pval_{pair}.{ont}.{testtype}.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pair=['Batwa.GBR-Anda.LWK', 'Bakiga.GBR-BR1.LWK'],
            ont=['BP', 'MF'],
            testtype=['overrep', 'enrich']),
    output:
        expand("results{cortype}/pbs.{pop}.pvals.BPMF.{testtype}.results.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pop=['Batwa.GBR', 'Bakiga.GBR', 'Anda.LWK', 'BR1.LWK'],
            testtype=['overrep', 'enrich']),
        expand("results{cortype}/convergent_permuted_pval_{pair}.BPMF.{testtype}.txt",
            cortype=['', '_sizeCor', '_mafCor'],
            pair=['Batwa.GBR-Anda.LWK', 'Bakiga.GBR-BR1.LWK'],
            testtype=['overrep', 'enrich'])
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/combine_BP_MF_results.sh"

# ----------------------------------------------------------------------------------------
# --- Plot GO results - Observed vs Expected for one stat/population
# ----------------------------------------------------------------------------------------

rule plot_go_results:
    input:
        expand("results/{root}.pvals.{ont}.overrep.results.txt",
            root=["pbs.Batwa.GBR","pbs.Bakiga.GBR","pbs.Anda.LWK","pbs.BR1.LWK",
                  "bayenv_BF_pass"],
            ont=ONTS),
    output:
        expand("reports/{root}.pvals.{ont}.GO.overrep.pdf",
            root=["pbs.Batwa.GBR","pbs.Bakiga.GBR","pbs.Anda.LWK","pbs.BR1.LWK",
                  "bayenv_BF_pass"],
            ont=ONTS),
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "sh scripts/plot_all_GO_results.sh"

# ----------------------------------------------------------------------------------------
# --- Plot GO results observed vs expected ratio
# ----------------------------------------------------------------------------------------

rule plot_go_results_obs_vs_exp:
    input:
        "results/pbs.Batwa.GBR.pvals.{ont}.overrep.results.txt",
        "results/pbs.Anda.LWK.pvals.{ont}.overrep.results.txt",
        "results/pbs.Bakiga.GBR.pvals.{ont}.overrep.results.txt",
        "results/pbs.BR1.LWK.pvals.{ont}.overrep.results.txt",
        "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.overrep.txt",
        "results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.overrep.txt",
    output:
        "reports/pbs.Batwa-Anda.pvals.{ont}.GO.overrep.pdf",
        "reports/pbs.Bakiga-BR1.pvals.{ont}.GO.overrep.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_GO_exp_vs_obs.R results/pbs.Batwa.GBR.pvals.{wildcards.ont}.overrep.results.txt results/pbs.Anda.LWK.pvals.{wildcards.ont}.overrep.results.txt results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{wildcards.ont}.overrep.txt reports/pbs.Batwa-Anda.pvals.{wildcards.ont}.GO.overrep.pdf;"
        "Rscript scripts/plot_GO_exp_vs_obs.R results/pbs.Bakiga.GBR.pvals.{wildcards.ont}.overrep.results.txt results/pbs.BR1.LWK.pvals.{wildcards.ont}.overrep.results.txt results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{wildcards.ont}.overrep.txt reports/pbs.Bakiga-BR1.pvals.{wildcards.ont}.GO.overrep.pdf;"

rule plot_go_results_obs_vs_exp_corrected:
    input:
        "results{cortype}/pbs.Batwa.GBR.pvals.{ont}.overrep.results.txt",
        "results{cortype}/pbs.Anda.LWK.pvals.{ont}.overrep.results.txt",
        "results{cortype}/pbs.Bakiga.GBR.pvals.{ont}.overrep.results.txt",
        "results{cortype}/pbs.BR1.LWK.pvals.{ont}.overrep.results.txt",
        "results{cortype}/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.overrep.txt",
        "results{cortype}/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.overrep.txt",
    output:
        "reports{cortype}/pbs.Batwa-Anda.pvals.{ont}.GO.overrep.pdf",
        "reports{cortype}/pbs.Bakiga-BR1.pvals.{ont}.GO.overrep.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_GO_exp_vs_obs.R results{wildcards.cortype}/pbs.Batwa.GBR.pvals.{wildcards.ont}.overrep.results.txt results{wildcards.cortype}/pbs.Anda.LWK.pvals.{wildcards.ont}.overrep.results.txt results{wildcards.cortype}/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{wildcards.ont}.overrep.txt reports{wildcards.cortype}/pbs.Batwa-Anda.pvals.{wildcards.ont}.GO.overrep.pdf;"
        "Rscript scripts/plot_GO_exp_vs_obs.R results{wildcards.cortype}/pbs.Bakiga.GBR.pvals.{wildcards.ont}.overrep.results.txt results{wildcards.cortype}/pbs.BR1.LWK.pvals.{wildcards.ont}.overrep.results.txt results{wildcards.cortype}/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{wildcards.ont}.overrep.txt reports{wildcards.cortype}/pbs.Bakiga-BR1.pvals.{wildcards.ont}.GO.overrep.pdf;"

# ----------------------------------------------------------------------------------------
# --- Plot GO enrichment results
# ----------------------------------------------------------------------------------------

rule plot_go_results_enrich:
    input:
        "results/pbs.Batwa.GBR.pvals.{ont}.enrich.results.txt",
        "results/pbs.Anda.LWK.pvals.{ont}.enrich.results.txt",
        "results/pbs.Bakiga.GBR.pvals.{ont}.enrich.results.txt",
        "results/pbs.BR1.LWK.pvals.{ont}.enrich.results.txt",
        "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.enrich.txt",
        "results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.enrich.txt",
    output:
        "reports/pbs.Batwa-Anda.pvals.{ont}.GO.enrich.pdf",
        "reports/pbs.Bakiga-BR1.pvals.{ont}.GO.enrich.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_GO_enrich_results.R results/pbs.Batwa.GBR.pvals.{wildcards.ont}.enrich.results.txt results/pbs.Anda.LWK.pvals.{wildcards.ont}.enrich.results.txt results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{wildcards.ont}.enrich.txt reports/pbs.Batwa-Anda.pvals.{wildcards.ont}.GO.enrich.pdf;"
        "Rscript scripts/plot_GO_enrich_results.R results/pbs.Bakiga.GBR.pvals.{wildcards.ont}.enrich.results.txt results/pbs.BR1.LWK.pvals.{wildcards.ont}.enrich.results.txt results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{wildcards.ont}.enrich.txt reports/pbs.Bakiga-BR1.pvals.{wildcards.ont}.GO.enrich.pdf;"

rule plot_go_results_enrich_corrected:
    input:
        "results{cortype}/pbs.Batwa.GBR.pvals.{ont}.enrich.results.txt",
        "results{cortype}/pbs.Anda.LWK.pvals.{ont}.enrich.results.txt",
        "results{cortype}/pbs.Bakiga.GBR.pvals.{ont}.enrich.results.txt",
        "results{cortype}/pbs.BR1.LWK.pvals.{ont}.enrich.results.txt",
        "results{cortype}/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{ont}.enrich.txt",
        "results{cortype}/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{ont}.enrich.txt",
    output:
        "reports{cortype}/pbs.Batwa-Anda.pvals.{ont}.GO.enrich.pdf",
        "reports{cortype}/pbs.Bakiga-BR1.pvals.{ont}.GO.enrich.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_GO_enrich_results.R results{wildcards.cortype}/pbs.Batwa.GBR.pvals.{wildcards.ont}.enrich.results.txt results{wildcards.cortype}/pbs.Anda.LWK.pvals.{wildcards.ont}.enrich.results.txt results{wildcards.cortype}/convergent_permuted_pval_Batwa.GBR-Anda.LWK.{wildcards.ont}.enrich.txt reports{wildcards.cortype}/pbs.Batwa-Anda.pvals.{wildcards.ont}.GO.enrich.pdf;"
        "Rscript scripts/plot_GO_enrich_results.R results{wildcards.cortype}/pbs.Bakiga.GBR.pvals.{wildcards.ont}.enrich.results.txt results{wildcards.cortype}/pbs.BR1.LWK.pvals.{wildcards.ont}.enrich.results.txt results{wildcards.cortype}/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.{wildcards.ont}.enrich.txt reports{wildcards.cortype}/pbs.Bakiga-BR1.pvals.{wildcards.ont}.GO.enrich.pdf;"

# ----------------------------------------------------------------------------------------
# --- Plot PBS by SNP along genome
# ----------------------------------------------------------------------------------------

rule plot_pbs_by_snp_along_genome:
    input:
        expand("results/{root}.anno.bed", root=GENE_LIST_ROOTS_PBS),
    output:
        "reports/PBS_by_SNP.pdf",
        "results/gene_SNP_counts.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_PBS_by_SNP.R"

# ----------------------------------------------------------------------------------------
# --- Plot PBS by gene along genome
# ----------------------------------------------------------------------------------------

rule plot_pbs_by_gene_along_genome:
    input:
        expand("results/{root}.pvals.txt", root=GENE_LIST_ROOTS_PBS),
        "results/gene_SNP_counts.txt"
    output:
        "reports/PBS_by_gene.pdf",
        "reports/PBS_by_gene_length.pdf",
        "reports/PBS_by_SNP_count.pdf",
        "reports/PBS_stddev_by_SNP_count.pdf",
        "reports/PBS_stderr_by_SNP_count.pdf",
        "reports/PBS_boxplots_by_SNP_count.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_PBS_by_gene.R"

# ----------------------------------------------------------------------------------------
# --- Plot convergent PBS p-value by SNP along genome
# ----------------------------------------------------------------------------------------

rule plot_convergent_pbs_by_snp_along_genome:
    input:
        expand("results/joint_SNP_pvalues.txt"),
    output:
        "reports/convergent_PBS_p_by_SNP.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_convergent_PBS_by_SNP.R"

# ----------------------------------------------------------------------------------------
# --- Plot BF by SNP along genome
# ----------------------------------------------------------------------------------------

rule plot_bf_by_snp_along_genome:
    input:
        expand("results/bayenv_BF_pass.anno.bed"),
    output:
        "reports/BF_by_SNP.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_BF_by_SNP.R"

# ----------------------------------------------------------------------------------------
# --- Plot BF by gene along genome
# ----------------------------------------------------------------------------------------

rule plot_bf_by_gene_along_genome:
    input:
        expand("results/bayenv_BF_pass.pvals.txt", root=GENE_LIST_ROOTS_PBS),
    output:
        "reports/BF_by_gene.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_BF_by_gene.R"

# ----------------------------------------------------------------------------------------
# --- Annotate PBS (and Fst and T) file with genes for PBS schematic generation
# ----------------------------------------------------------------------------------------

rule annotate_pbs_for_schematic:
    input:
        "results/GBR.FST_for_PBS.txt",
        "results/LWK.FST_for_PBS.txt"
    output:
        "results/GBR.FST_for_PBS.txt.anno.bed",
        "results/LWK.FST_for_PBS.txt.anno.bed"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "sh scripts/annotate_PBS_file_with_genes_for_schematic.sh"

# ----------------------------------------------------------------------------------------
# --- Plot PBS schematic
# ----------------------------------------------------------------------------------------

rule plot_pbs_schematic:
    input:
        "results/GBR.FST_for_PBS.txt.anno.bed",
        "results/LWK.FST_for_PBS.txt.anno.bed"
    output:
        "reports/pbs_schematic.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_pbs_schematic.R"

# ----------------------------------------------------------------------------------------
# --- Plot map
# ----------------------------------------------------------------------------------------

rule plot_map:
    output:
        "reports/sampling_map.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_map.R"

# ========================================================================================
# --- Gather results and make tables
# ========================================================================================

rule parse_apriori_gene_list_results:
    input:
        expand("reports/{root}.pvals.fisher.overrep.txt",
            root=GENE_LIST_ROOTS),
        expand("reports/{root}.pvals.wilcoxshift.txt",
            root=GENE_LIST_ROOTS),
    output:
        "reports/parsed_apriori_gene_list_results.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=2gb"
    shell:
        "sh scripts/parse_apriori_gene_list_results.sh > {output}"

# ========================================================================================
# --- TO DO:
# ========================================================================================

### Add step to make data/hg19.genome
### Add step to download 1000G data for impute run:
###    {DIR_1000G}/genetic_map_chr*_combined_b37.txt
### Add step to make list of interesting GO terms:
###    data/GO_of_interest.txt

# ========================================================================================
# --- Set default for optional mem parameter
# ========================================================================================

default = ""
name = 'mem'
for r in workflow.rules:
    try:
        getattr(r.params, name)
    except AttributeError:
        r.params.append(default)
        r.params.add_name(name)

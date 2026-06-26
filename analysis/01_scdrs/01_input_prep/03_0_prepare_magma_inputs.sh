# Step 1: download MAGMA software, gene location file, and reference data from
# https://ctg.cncr.nl/software/magma after this step, one should have a folder <MAGMA_DIR>
# with the following files:
# 1) <MAGMA_DIR>/magma 2) <MAGMA_DIR>/g1000_eur.(bed|bim|fam) 3) <MAGMA_DIR>/NCBI37.3.gene.loc
magma_dir=/mnt/isilon/gandal_lab/liaoyd/tools/magma
dat_dir=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/GWAS_result
out_dir=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/GWAS_result/MAGMA

# Step 2: make gene annotation file for MAGMA using the following command, this only needs to be done once for different GWAS summary statistics, and the file will be saved to out/step1.genes.annot
mkdir -p ${out_dir}
# GRCh37
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/g1000_eur.bim \
    --gene-loc ${magma_dir}/NCBI37.3.gene.loc \
    --out ${out_dir}/step1_37

# GRCh38
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/LDSC_GRCH38_plink_files/g1000_eur_hg38.bim \
    --gene-loc ${magma_dir}/NCBI38.gene.loc \
    --out ${out_dir}/step1_38

# Step 3: run MAGMA using the following command, this takes a GWAS file ${trait}.pval,
# which at least has the columns: SNP, P, N, which corresponds to the SNP id
# (matched to the ${magma_dir}/g1000_eur.bim), p-value, sample size. For example,
# <trait>.pval file looks like
#
# CHR     BP      SNP             P           N
# 1       717587  rs144155419     0.453345    279949
# 1       740284  rs61770167      0.921906    282079
# 1       769223  rs60320384      0.059349    281744
#
# After this step, one should obtain a file ${out_dir}/step2/${trait}.gene.out, and the top genes with largest Z-scores can be input to scDRS.

# GRCh37
mkdir -p ${out_dir}/step2
trait="PGC_ASD_2019"
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${dat_dir}/${trait}.pval N=46351 \
    --gene-annot ${out_dir}/step1_37.genes.annot \
    --out ${out_dir}/step2/${trait}

export trait="PGC_ASD_2019"
export out_dir=${out_dir}/
Rscript 02_prepare_magma_gene_scores_grch37.R


# GRCh38
trait="ASD_Matoba_2020"
${magma_dir}/magma \
    --bfile ${magma_dir}/LDSC_GRCH38_plink_files/g1000_eur_hg38 \
    --pval ${dat_dir}/${trait}.pval N=55420 \
    --gene-annot ${out_dir}/step1_38.genes.annot \
    --out ${out_dir}/step2/${trait}


export trait="ASD_Matoba_2020"
export out_dir=${out_dir}
Rscript 02_prepare_magma_gene_scores_grch38.R
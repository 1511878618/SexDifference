#!/bin/bash

hardCallPath=/pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/ukb_cal_allChrs
threads=$1
output=$2

PhenoPath=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff

cat >&1 <<-EOF
Example code:
./step1.sh 30 step1/
EOF

mkdir -p ${output}
# run bt step1

## run step1 for Prot 
## As 2911 is too large, split into several part

# for chunk in {2..6};do  # totally 16; first run 6; then 7..11; finaly 11..16 
#     protChunkOutPutDir=${output}/Prot_part${chunk}
#     mkdir -p ${protChunkOutPutDir}
#     sbatch -J "Prot${chunk}_step1" -c 20 --mem=20G -o step1_Prot_${chunk}.log --wrap """
#     regenie \
#         --step 1 \
#         --threads ${threads} \
#         --bed ${hardCallPath} \
#         --extract /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.snplist \
#         --keep /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.id \
#         --keep /pmaster/xutingfeng/dataset/ukb/phenotype/white.tsv \
#         --qt \
#         --phenoFile ${PhenoPath}/Prot_part${chunk}.regenie \
#         --covarFile /pmaster/xutingfeng/dataset/ukb/phenotype/regenie.cov \
#         --covarColList genotype_array,inferred_sex,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
#         --catCovarList	genotype_array,inferred_sex,assessment_center \
#         --maxCatLevels 30 \
#         --bsize 1000 \
#         --lowmem \
#         --lowmem-prefix /hwmaster/xutingfeng/Prot_part${chunk} \
#         --out ${protChunkOutPutDir}
#         """
# done 
## Step1 for NMR 
for chunk in {2..2};do # totally 2
    nmrChunkOutPutDir=${output}/NMR_part${chunk}
    mkdir -p ${output}/nmrChunkOutPutDir
    sbatch -J "NMR${chunk}_step1" -c ${threads} --mem=20G -o step1_NMR_${chunk}.log --wrap """
    regenie \
        --step 1 \
        --threads ${threads} \
        --bed ${hardCallPath} \
        --extract /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.snplist \
        --keep /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.id \
        --keep /pmaster/xutingfeng/dataset/ukb/phenotype/white.tsv \
        --qt \
        --phenoFile ${PhenoPath}/NMR_part${chunk}.regenie \
        --covarFile /pmaster/xutingfeng/dataset/ukb/phenotype/regenie.cov \
        --covarColList genotype_array,inferred_sex,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
        --catCovarList	genotype_array,inferred_sex,assessment_center \
        --maxCatLevels 30 \
        --bsize 1000 \
        --lowmem \
        --lowmem-prefix /hwmaster/xutingfeng/NMR_part${chunk} \
        --out ${nmrChunkOutPutDir}
        """
done 
## Step1 for RF
# mkdir -p ${output}/bloodAssayAndBasicCharacteristics
# sbatch -J "RF_step1" -c ${threads} --mem=25G -o step1_bloodAssayAndBasicCharacteristics.log --wrap """
# regenie \
#     --step 1 \
#     --threads ${threads} \
#     --bed ${hardCallPath} \
#     --extract /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.snplist \
#     --keep /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.id \
#     --keep /pmaster/xutingfeng/dataset/ukb/phenotype/white.tsv \
#     --qt \
#     --phenoFile ${PhenoPath}/bloodAssayAndBasicCharacteristics.regenie \
#     --covarFile /pmaster/xutingfeng/dataset/ukb/phenotype/regenie.cov \
#     --covarColList genotype_array,inferred_sex,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
#     --catCovarList	genotype_array,inferred_sex,assessment_center \
#     --maxCatLevels 30 \
#     --bsize 1000 \
#     --out ${output}/bloodAssayAndBasicCharacteristics/
#     """

# Step1 for Disease 
mkdir -p ${output}/MainDiseasePheno
sbatch -J "MainDiseasePheno_step1" -c ${threads} --mem=25G -o step1_MainDiseasePheno.log \
    --wrap "
regenie \
    --step 1 \
    --threads ${threads} \
    --bed ${hardCallPath} \
    --extract /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.snplist \
    --keep /pmaster/xutingfeng/dataset/ukb/dataset/snp/geneArray/qc_pass.id \
    --keep /pmaster/xutingfeng/dataset/ukb/phenotype/white.tsv \
    --bt \
    --phenoFile ${PhenoPath}/MainDiseasePheno.regenie \
    --covarFile /pmaster/xutingfeng/dataset/ukb/phenotype/regenie.cov \
    --covarColList genotype_array,inferred_sex,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
    --catCovarList	genotype_array,inferred_sex,assessment_center \
    --maxCatLevels 30 \
    --firth \
    --approx \
    --pThresh 0.01 \
    --bsize 1000 \
    --out ${output}/MainDiseasePheno
    "

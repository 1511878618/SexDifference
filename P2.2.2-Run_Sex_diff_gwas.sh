#!/bin/bash
## group by sex and run regenie
outputPath="./GWASResultV2"
threads=5
memory="15G"

#params
chrDir=/pmaster/xutingfeng/dataset/ukb/dataset/snp/ukb_imputed_v3_qc/qc_pgen_hg19
covarFile=/pmaster/xutingfeng/dataset/ukb/phenotype/regenie.cov

threads=5

# # Run Disease 
# step1=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/step1/MainDiseasePheno/_pred.list
# DiseaseOutputPath=${outputPath}/Disease
# phenoFile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/disease.regenie
# logfile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/log/Burden/disease/

# mkdir -p ${logfile}


# for sex in female male; do
#     keep_files=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/white_train_${sex}.tsv
#     c_outputPath=${DiseaseOutputPath}/${sex}
#     mkdir -p ${c_outputPath}

#     for chr in {1..22}; do  # 22

#         sbatch -J "${sex}_chr${chr}" -c ${threads} --mem=${memory} -o ${logfile}/${sex}_chr${chr}_gwas.log --wrap """
#         regenie --step 2 \
#         --threads=${threads} \
#         --ref-first \
#         --pgen ${chrDir}/ukb_imp_chr${chr}_v3_hg37_qc\
#         --phenoFile ${phenoFile} \
#         --keep ${keep_files} \
#          --bt \
#         --covarFile ${covarFile} \
#         --covarColList genotype_array,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
#         --catCovarList genotype_array,assessment_center \
#         --maxCatLevels 30 \
#         --bsize 1000 \
#         --out ${c_outputPath}/chr${chr} \
#         --minMAC 10 \
#        --firth \
#         --approx \
#         --pThresh 0.01 \
#         --pred ${step1}"""
#     done
# done



# ## Run Prot 

# step1=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/step1/Prot/_pred.list
# ProtOutputPath=${outputPath}/Prot
# phenoFile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/Prot.regenie
# logfile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/log/GWAS/Prot/

# mkdir -p ${logfile}
# for sex in female male; do
#     keep_files=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/white_train_${sex}.tsv
#     c_outputPath=${ProtOutputPath}/${sex}
#     mkdir -p ${c_outputPath}

#     for chr in {2..22}; do  # 22

#     echo sbatch -J "${sex}_chr${chr}" -c ${threads} --mem=${memory} -o ${logfile}/${sex}_chr${chr}_gwas.log --wrap """
#         regenie --step 2 \
#         --threads=${threads} \
#         --ref-first \
#         --pgen ${chrDir}/ukb_imp_chr${chr}_v3_hg37_qc\
#         --phenoFile ${phenoFile} \
#         --keep ${keep_files} \
#          --qt \
#         --covarFile ${covarFile} \
#         --covarColList genotype_array,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
#         --catCovarList genotype_array,assessment_center \
#         --maxCatLevels 30 \
#         --bsize 1000 \
#         --out ${c_outputPath}/chr${chr} \
#         --minMAC 10 \
#         --pred ${step1}"""
#     done
# done


# ## Run NMR 

# step1=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/step1/NMR/_pred.list
# NMROutputPath=${outputPath}/NMR
# phenoFile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/NMR.regenie
# logfile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/log/GWAS/NMR/

# mkdir -p ${logfile}
# for sex in female male; do
#     keep_files=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/white_train_${sex}.tsv
#     c_outputPath=${NMROutputPath}/${sex}
#     mkdir -p ${c_outputPath}

#     for chr in {1..22}; do  # 22

#     sbatch -J "${sex}_chr${chr}" -c ${threads} --mem=${memory} -o ${logfile}/${sex}_chr${chr}_gwas.log --wrap """
#         regenie --step 2 \
#         --threads=${threads} \
#         --ref-first \
#         --pgen ${chrDir}/ukb_imp_chr${chr}_v3_hg37_qc\
#         --phenoFile ${phenoFile} \
#         --keep ${keep_files} \
#          --qt \
#         --covarFile ${covarFile} \
#         --covarColList genotype_array,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
#         --catCovarList genotype_array,assessment_center \
#         --maxCatLevels 30 \
#         --bsize 1000 \
#         --out ${c_outputPath}/chr${chr} \
#         --minMAC 10 \
#         --pred ${step1}"""
#     done
# done

## Run RF 

step1=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/step1/bloodAssayAndBasicCharacteristics/_pred.list
RFOutputPath=${outputPath}/bloodAssayAndBasicCharacteristics
phenoFile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/bloodAssayAndBasicCharacteristics.regenie
logfile=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/log/GWASV2/bloodAssayAndBasicCharacteristics/

mkdir -p ${logfile}
for sex in female male; do
    keep_files=/pmaster/xutingfeng/dataset/ukb/dataset/regenie/sex_diff/white_train_${sex}.tsv
    c_outputPath=${RFOutputPath}/${sex}
    mkdir -p ${c_outputPath}

    for chr in {1..22}; do  # 22

    sbatch -J "${sex}_chr${chr}" -c ${threads} --mem=${memory} -o ${logfile}/${sex}_chr${chr}_gwas.log --wrap """
        regenie --step 2 \
        --threads=${threads} \
        --ref-first \
        --pgen ${chrDir}/ukb_imp_chr${chr}_v3_hg37_qc\
        --phenoFile ${phenoFile} \
        --keep ${keep_files} \
         --qt \
        --covarFile ${covarFile} \
        --covarColList genotype_array,age_visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,assessment_center,age_squared \
        --catCovarList genotype_array,assessment_center \
        --maxCatLevels 30 \
        --bsize 1000 \
        --out ${c_outputPath}/chr${chr} \
        --minMAC 10 \
        --pred ${step1}"""
    done
done







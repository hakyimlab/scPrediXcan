## scPrediXcan: Leveraging Single-Cell Data for Cell-Type–Specific Transcriptome-Wide Association Studies Using Transfer Learning  
#
<p align="center">
  <img height="560" src="Figures/scPrediXcan_workflow.png">
</p>

#
### Description 
Single-cell PrediXcan(scPrediXcan) is framework to perform TWAS at the cell-type level using single-cell data. The workflow of scPrediXcan framework consists of three steps: 1) training a deep learning model named ctPred for epigenomics-to-expression prediction at cell type level, 2) linearizing the ctPred into a SNP-based elastic net model for downstream association tests using GWAS summary statistics, and 3) performing the association tests between genes and trait of interest.

#
### Setup and installation

#
### Usage

#### Step1: Training the ctPred model  

#### Step2: Linearining the ctPred into l-ctPred  
scPrediXcan uses [PrediXcan implementation](https://www.nature.com/articles/ng.3367) to train an elastic-net model for ctPred linearization.
Here are the detailed procedures of step2:

1) Install nextflow into the environment.
```bash
conda install nextflow
```
2) Clone the PredictDb-nextflow repository.
```bash
$ git clone https://github.com/hakyimlab/PredictDb-nextflow.git
```

3) Run the PredictDb nextflow pipeline.
```bash

nextflow run ./main.nf \
--gene_annotation 'Gene_anno.txt' \
--snp_annotation 'snp_annot.txt' \
--genotype 'genotype_file' \
--gene_exp 'ctPred_predicted_gene_expression.csv' \
--outdir results \
--keepIntermediate \
-resume \

```
#
#### Step3: Performing association test between genes and traits 

scPrediXcan uses Summary-PrediXcan(S-PrediXcan) to run the association test. The detailed description of S-PrediXcan are [here](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial). In this step, the input data include: an appropiate **Transcriptome Model Database (i.g., l-ctPred)**, a **GWAS/Meta Analysis summary statistics**, and **SNP covariance matrices**. The l-ctPred database and the SNP covariance matrices are obtained from the last step. Here are the detailed procedures of step3:

1) Clone the S-PrediXcan repository.
```bash
$ git clone https://github.com/hakyimlab/MetaXcan
```

2) Go to the software folder.
```bash
$ cd MetaXcan/software
```

3) Download example [data](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9).

This may take a few minutes depending on your connection: it has to download approximately 200Mb worth of data.
Downloaded data will include an appropiate **Transcriptome Model Database (i.g., l-ctPred)**, a **GWAS/Meta Analysis summary statistics**, and **SNP covariance matrices**.

Extract it with:
```bash
tar -xzvpf sample_data.tar.gz
```

4) Run the High-Level S-PrediXcan Script
```bash
./SPrediXcan.py \
--model_db_path data/DGN-WB_0.5.db \
--covariance data/covariance.DGN-WB_0.5.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/test.csv
```
This should take less than a minute on a 3GHZ computer. For the full specification of command line parameters, you can check the [wiki](https://github.com/hakyimlab/MetaXcan/wiki/MetaXcan's-Command-Line-Reference).


The example command parameters mean:

* `--model_db_path` Path to tissue transriptome model
* `--covariance` Path to file containing covariance information. This covariance should have information related to the tissue transcriptome model.
* `--gwas_folder` Folder containing GWAS summary statistics data.
* `--gwas_file_pattern` This option allows the program to select which files from the input to use based on their name.
...This allows to ignore several support files that might be generated at your GWAS analysis, such as plink logs.
* `--snp_column` Argument with the name of the column containing the RSIDs.
* `--effect_allele_column` Argument with the name of the column containing the effect allele (i.e. the one being regressed on).
* `--non_effect_allele_column` Argument with the name of the column containing the non effect allele.
* `--beta_column` Tells the program the name of a column containing -phenotype beta data for each SNP- in the input GWAS files.
* `--pvalue_column` Tells the program the name of a column containing -PValue for each SNP- in the input GWAS files.
* `--output_file` Path where results will be saved to.

Its output is a CSV file that looks like:

```
gene,gene_name,zscore,effect_size,pvalue,var_g,pred_perf_r2,pred_perf_pval,pred_perf_qval,n_snps_used,n_snps_in_cov,n_snps_in_model
ENSG00000150938,CRIM1,4.190697619877402,0.7381499095142079,2.7809807629839122e-05,0.09833448081630237,0.13320775358,1.97496173512e-30,7.47907447189e-30,37,37,37
...
```
Where each row is a gene's association result:
* `gene`: a gene's id: as listed in the Tissue Transcriptome model.
Ensemble Id for most gene model releases. Can also be a intron's id for splicing model releases.
* `gene_name`: gene name as listed by the Transcriptome Model, typically HUGO for a gene. It can also be an intron's id.
* `zscore`: S-PrediXcan's association result for the gene, typically HUGO for a gene.
* `effect_size`: S-PrediXcan's association effect size for the gene. Can only be computed when `beta` from the GWAS is used.
* `pvalue`: P-value of the aforementioned statistic.
* `pred_perf_r2`: (cross-validated) R2 of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `pred_perf_pval`: pval of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `pred_perf_qval`: qval of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `n_snps_used`: number of snps from GWAS that got used in S-PrediXcan analysis
* `n_snps_in_cov`: number of snps in the covariance matrix
* `n_snps_in_model`: number of snps in the model
* `var_g`: variance of the gene expression, calculated as `W' * G * W`
(where `W` is the vector of SNP weights in a gene's model,
`W'` is its transpose, and `G` is the covariance matrix)

If `--additional_output` is used when running S-PrediXcan, you'll get two additional columns:
* `best_gwas_p`: the highest p-value from GWAS snps used in this model
* `largest_weight`: the largest (absolute value) weight in this model


# 
### Citation

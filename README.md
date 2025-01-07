## scPrediXcan: Leveraging Single-Cell Data for Cell-Type–Specific Transcriptome-Wide Association Studies Using Transfer Learning  
<p align="center">
  <img height="600" src="Figures/scPrediXcan_workflow.png">
</p>

#
### Description 
Single-cell PrediXcan (scPrediXcan) is a framework designed to perform Transcriptome-Wide Association Studies (TWAS) at the cell-type level using single-cell data. This framework utilizes GWAS summary statistics and single-cell RNA-seq data to assess the association between gene expression and disease risk.

The workflow of scPrediXcan framework consists of three steps: 1) training a deep learning model named ctPred for epigenomics-to-expression prediction at cell type level, 2) linearizing the ctPred into a SNP-based elastic net model for downstream association tests using GWAS summary statistics, and 3) performing the association tests between genes and trait of interest. 

The framework is detailed in the [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.11.11.623049v1).

You can download the cell type specific prediction models and covariances [in links from here](https://predictdb.org/post/2024/10/14/scpredixcan-prediction-models/)

#
### Setup and installation
scPrediXcan mainly requires Python, Nextflow and some deep learning packages like Pytorch and Tensorflow. You can install all the packages and softwares needed by creating a new conda environment named **scPrediXcan** using scPrediXcan_env.yml file.
```bash
conda env create -f scPrediXcan_env.yml
conda activate scPrediXcan

```

#
### Usage

#### Step1: Training the ctPred model  
!!Note that if you just want to use the pre-trained models of given cell types for TWAS in certain diseases, please skip step1,2 and directly do step3 with pre-trained l-ctPred models.

ctPred is a multilayer perceptron to predict gene expressions at pseudobulk level from gene epigenomic representations(i.g., Enformer-output epigenomic features). The inputs of this step are a population-average gene expression file and a gene epigenomic features file. They are combined into a single training data csv file since both are relatively small. The output of this step is a pt file storing the model weights.

```bash
python ctPred_train.py --parameters ctPred_train.json --cell_file 'training_data.csv'

```

Here is an example of the training data [file](./Scripts/ctPred/Example_training_data.csv), which is a bit different from the example training data file in the [tutorial](./Scripts/ctPred/Tutorial.ipynb). In the tutorial, the expression data has already been combined to the epigenomics data; here, you only need to provide the expression data and specify the path of epigenomics data path in the json file, and the script will combine them.

or 


```bash
python ctPred_train.py --parameters ctPred_train.json --data_dir 'The_path_of_the_folder' # inside the folder, you have many training data files (.csv) for different cell types

```

For the details, check the code and tutorial [here](https://github.com/hakyimlab/scPrediXcan/tree/master/Scripts/ctPred). The Enformer-predicted epigenomic features of protein-coding genes are shared [here](https://uchicago.app.box.com/folder/289080191619?s=wgyt3wqfjezsqwjhndcceky6nivil8kh).

#
#### Step2: Linearizing the ctPred into l-ctPred  
scPrediXcan uses [PrediXcan implementation](https://www.nature.com/articles/ng.3367) to train an elastic-net model for ctPred linearization. In this step, we utilize the genotype data from 448 Geuvadis individuals along with ctPred-predicted gene expression profiles to fit an elastic-net model for the corresponding cell type. In principle, alternative genotype reference panels can also be employed at this stage. 

This is our [pipeline](https://github.com/hakyimlab/shared_pipelines/tree/main/enformer_pipeline) to generate personalized epigenomic features using Enformer with genotype data (.vcf file). Once users calculate the personalized epigenomics, please use the ctPred models from step 1 to get the in-silico personalized gene expressions. Then, users can generate l-ctPred according to the following instructions.

Here is a nextflow pipeline for l-ctPred generation. The inputs include a genotype file and a ctPred-predicted cell-type-specific gene expression file. The outputs consist of a transcriptome model SQLite database (i.e., l-ctPred) and a SNP covariance matrix file. These output files are intended for use in the final association analysis step.
Here are the detailed procedures of step-2:

1) Clone the PredictDb-nextflow repository.
```bash
git clone https://github.com/hakyimlab/PredictDb-nextflow.git
```

2) Run the PredictDb nextflow pipeline.
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
The detailed descriptions of the pipeline and the used data/output are [here](https://github.com/hakyimlab/PredictDb-nextflow/blob/master/docs/usage.md).

#
#### Step3: Performing association test between genes and traits 

scPrediXcan uses Summary-PrediXcan(S-PrediXcan) to run the association test. The detailed description of S-PrediXcan are [here](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial). In this step, the input data include: a transcriptome model sqlite database (i.g., l-ctPred), a GWAS/Meta Analysis summary statistics, and SNP covariance matrices. The l-ctPred database and the SNP covariance matrices are obtained from the last step. Here are the detailed procedures of step-3:

1) Clone the S-PrediXcan repository and go to the software folder.
```bash
git clone https://github.com/hakyimlab/MetaXcan
cd MetaXcan/software
```

2) Run the High-Level S-PrediXcan Script
```bash
./SPrediXcan.py \
--model_db_path 'l-ctPred_celli.db' \
--covariance 'covariance.txt.gz' \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file 'results/TWAS_result.csv'
```

Pipeline details:

This step should take less than a minute on a 3GHZ computer. For the full specification of command line parameters, you can check the [wiki](https://github.com/hakyimlab/MetaXcan/wiki/MetaXcan's-Command-Line-Reference) and the [turtorial](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial). The output csv file is the TWAS result, and the detailed descriptions of each column are [here](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial)

You can download example data [here](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9). This may take a few minutes depending on your connection: it has to download approximately 200Mb worth of data. Downloaded data will include all the input data needed.

#
### Example files:

Some example input/output files from step1 to step3 are in this [box folder](https://uchicago.box.com/s/hsrb9zenyhpckm6cemxl96if02uy8peu).


# 
### Citation
If you find this code useful, we would appreciate it if you cite the [following publication](https://www.biorxiv.org/content/10.1101/2024.11.11.623049v1)

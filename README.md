## scPrediXcan: Leveraging Single-Cell Data for Cell-Type–Specific Transcriptome-Wide Association Studies Using Transfer Learning
<p align="center">
  <img height="560" src="Figures/scPrediXcan_workflow.png">
</p>

### Description
Single-cell PrediXcan(scPrediXcan) is framework to perform TWAS at the cell-type level using single-cell data. The workflow of scPrediXcan framework consists of three steps: 1) training a deep learning model named ctPred for epigenomics-to-expression prediction at cell type level, 2) linearizing the ctPred into a SNP-based elastic net model for downstream association tests using GWAS summary statistics, and 3) performing the association tests between genes and trait of interest.
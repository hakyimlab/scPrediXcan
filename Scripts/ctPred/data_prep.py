# Import Pandas
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--exp', type=str, required=True, help='Path to expression CSV')
parser.add_argument('--epi', type=str, required=True, help='Path to epigenomics CSV')
parser.add_argument('--output', type=str, default='./training_cell.csv', help='Path to save output')
args = parser.parse_args()

Exp = pd.read_csv(args.exp)
epi = pd.read_csv(args.epi, index_col=0)

# Select the genes in the epigenomic features file and fill those 
# genes not detected in the single-cell data with 0, otherwise it # may affect downstream percentile calculation.
selected_exp = epi.iloc[:, [0]].merge(Exp, on='gene_name', how='left').fillna(0)

# Calculate mean expression
cols_to_average = selected_exp.columns[1:] 
selected_exp['mean_expression'] = selected_exp[cols_to_average].mean(axis=1)
selected_exp = selected_exp.drop(columns=cols_to_average)
# Convert to rank-based percentiles and merge with epigenomics file 
selected_exp['mean_expression'] = selected_exp['mean_expression' ].rank(method='average', pct=True)
# Merge expression percentiles with epigenomics file 
data_training = epi.merge(selected_exp, on='gene_name',how='left') 
# Save the file and release memory 
data_training.to_csv(args.output)
print(f"Training data saved to {args.output}")
# Here you can choose the path you like for the training csv file, and name it with the cell type/state.
del [selected_exp, cols_to_average, Exp, epi]
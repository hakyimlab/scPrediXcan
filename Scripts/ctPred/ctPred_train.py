import torch
from torch.utils import data
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from ctPred_utils import *
import os
import argparse
import json
mpl.rcParams['pdf.fonttype'] = 42

# specify the device that you'll use cpu or gpu
device = torch.device('cpu')



def ctPred_train(exp_matrix_p, params):

    # Load parameters from JSON
    Model_path = params.get("Model_path")
    fig_path = params.get("fig_path")
    train_set = params.get("training_set")
    val_set = params.get("valid_set")
    test_set = params.get("test_set")
    epi_p = exp_matrix_p 

    cell_type = os.path.basename(exp_matrix_p).split('.csv')[0]
    data = pd.read_csv(epi_p, index_col = 0) # assume that the epigenomics and expression files are merged already

    train_epi, train_exp, val_epi, val_epi, val_exp, test_epi, test_exp = data_prepare(data, train_set, val_set, test_set)
    train_data_iter = dataloader(train_epi, train_exp, batch_size=1000)

    # train the model and save it
    ctPred_test = ctPred().to(device)

    saved_path = os.path.join(Model_path, f'{cell_type}_ctPred.pt')
    print(f"Trained model saved to {saved_path}")
    
    ctPred_training(ctPred_test, train_data_iter, train_epi, train_exp, val_epi, val_exp, cell_type, fig_path, epochs=100, model_path=saved_path)

    plot_prediction(ctPred_test, test_epi, test_exp, cell_type, fig_path)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Train Model for a Cell File')
    parser.add_argument('--parameters', type=str, help='Path to the JSON parameter file')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--cell_file', type=str, help='Path to the cell file')
    group.add_argument('--data_dir', type=str, help='Path to the directory containing cell files')
    args = parser.parse_args()

    # Load parameters from JSON file
    if args.parameters:
        with open(args.parameters) as f:
            params = json.load(f)
    else:
        raise ValueError("No JSON parameter file provided.")


    # If --cell_file argument is provided
    if args.cell_file:
        cell_file = args.cell_file
        ctPred_train(cell_file, params)

    # If --data_dir argument is provided
    elif args.data_dir:
        data_dir = args.data_dir
        cell_files = [os.path.join(data_dir, file) for file in os.listdir(data_dir) if file.endswith(".csv")]
        for cell_file in cell_files:
            ctPred_train(cell_file, params)



if __name__ == "__main__":
    main()
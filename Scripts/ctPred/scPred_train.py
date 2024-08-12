import torch
from torch.utils import data
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scPred_utils import *
import os
import argparse
import json
mpl.rcParams['pdf.fonttype'] = 42

# specify the device that you'll use cpu or gpu
device = torch.device('cpu')



def scPred_train(data, params):

    # Load parameters from JSON
    Model_path = params.get("Model_path")
    fig_path = params.get("fig_path")
    train_set = params.get("training_set")
    val_set = params.get("valid_set")
    test_set = params.get("test_set")

    train_epi, train_exp, val_epi, val_epi, val_exp, test_epi, test_exp = data_prepare(data, train_set, val_set, test_set)
    train_data_iter = dataloader(train_epi, train_exp, batch_size=1000)

    # train the model and save it
    scPred_test = scPred().to(device)

    saved_path = os.path.join(Model_path, 'scPred_new.pt')
    scPred_training(scPred_test, train_data_iter, train_epi, train_exp, val_epi, val_exp, fig_path, epochs=100, model_path=saved_path)

    plot_prediction(scPred_test, test_epi, test_exp, fig_path)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Train Model for a Cell File')
    parser.add_argument('--parameters', type=str, help='Path to the JSON parameter file')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--cell_file', type=str, help='Path to the cell file')
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
        gen_data = pd.read_csv(cell_file, index_col=0)

        scPred_train(gen_data, params)




if __name__ == "__main__":
    main()
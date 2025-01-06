import torch
from torch.utils import data
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from scipy.stats import pearsonr
import numpy as np
from sklearn.linear_model import LinearRegression
mpl.rcParams['pdf.fonttype'] = 42

# specify the device that you'll use cpu or gpu
device = torch.device('cuda')


## This is the code to convert the raw expression values to rank-based values (percentile)
# cell_i_matx['mean_expression'] = cell_i_matx['mean_expression'].rank(method = 'average', pct=True)

def input_prep(epi_p, exp_mat_p):
    test_p = exp_mat_p
    test = pd.read_csv(test_p, index_col = 0)
    test = test.drop_duplicates(subset='gene_name')
    test = test.dropna(subset='gene_name')

    epi = pd.read_csv(epi_p, index_col=0)

    test_mat = test.merge(epi, on = 'gene_name', how = 'left')
    exp_col = test_mat['mean_expression'] 
    test_mat = test_mat.drop(columns = ['mean_expression'])
    test_mat['mean_expression'] = exp_col

    return test_mat


def data_prepare(gen_data, train_set, val_set, test_set, is_normalization=True):
    # split the data
    train_data = gen_data[gen_data['chromo'].isin(train_set)]
    val_data = gen_data[gen_data['chromo'].isin(val_set)]
    test_data = gen_data[gen_data['chromo'].isin(test_set)]


    # choose naive B cell as an example to train a model for gene expression prediction
    train_epi = torch.tensor(train_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    train_exp = torch.tensor(train_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    val_epi = torch.tensor(val_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    val_exp = torch.tensor(val_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    test_epi = torch.tensor(test_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    test_exp = torch.tensor(test_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    if is_normalization:
        # Calculate mean and standard deviation on the training data
        train_mean = train_epi.mean(dim=0)
        train_std = train_epi.std(dim=0)

        # Normalize the input data
        train_epi = (train_epi - train_mean) / train_std
        val_epi = (val_epi - train_mean) / train_std
        test_epi = (test_epi - train_mean) / train_std

        # Ensure the data is in torch.float32 and on the correct device
        train_epi = train_epi.to(dtype=torch.float32, device=device)
        val_epi = val_epi.to(dtype=torch.float32, device=device)
        test_epi = test_epi.to(dtype=torch.float32, device=device)

    return train_epi, train_exp, val_epi, val_epi, val_exp, test_epi, test_exp



# the dataloader
def dataloader(features, labels, batch_size, is_train = True):
    dataset = data.TensorDataset(features, labels)
    return data.DataLoader(dataset, batch_size, shuffle = is_train)


# ctPred model
class ctPred(nn.Module):
    def __init__(self, **kwargs):
        super().__init__()
        
        ctPred_defaults = {
            'num_layers' : 4,
            'input_dim' : 5313,
            'hidden_dim' : 64,
            'output_dim' : 1,
            'reg_lambda' : 5e-4,
            'dropout_rate' : 0.05,
            'learning_rate' : 9e-5,
            'random_seed' : 1024
        }

        ctPred_defaults.update(kwargs)

        for key, value in ctPred_defaults.items():
            setattr(self, key, value)


        torch.manual_seed(self.random_seed)


        # model main
        layers = [nn.Linear(self.input_dim, self.hidden_dim), nn.ReLU(), nn.Dropout(self.dropout_rate)]
        hidden_layer = [nn.Linear(self.hidden_dim, self.hidden_dim), nn.ReLU(), nn.Dropout(self.dropout_rate)]
        
        for _ in range(self.num_layers - 1):
            layers.extend(hidden_layer)
        
        layers.append(nn.Linear(self.hidden_dim, self.output_dim))
        
        self.net = nn.Sequential(*layers)

    
    def custom_loss(self, y_true, y_pred):
        return F.mse_loss(y_true.reshape(-1, 1), y_pred.reshape(-1, 1))

    def forward(self, x):
        return self.net(x)
    
    def compile(self):
        self.optimizer = optim.Adam(self.parameters(), lr = self.learning_rate, weight_decay = self.reg_lambda)




def plot_loss_curve(train_losses, val_losses, epochs, cell_type, fig_folder):

    plt.figure(figsize=(8, 6))
    plt.plot(range(1, epochs + 1), train_losses, label='Train_loss')
    plt.plot(range(1, epochs + 1), val_losses, label='Val_loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Loss curve')
    plt.legend()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(os.path.join(fig_folder, f'{cell_type}_loss.pdf'), bbox_inches = 'tight')
    #plt.show()


def save_model(model, path):
    torch.save(model.state_dict(), path)


def ctPred_training(model, train_data_iter, train_epi, train_exp, val_epi, val_exp, cell_type, fig_folder, epochs=80, plot_loss=True, model_path = 'ctPred.pt'):
    model.compile()
    optimizer = model.optimizer

    train_losses = []
    val_losses = []
    best_val_loss = float('inf')
    
    for _ in range(epochs):
        
        model.eval()
        with torch.no_grad():
            val_loss = model.custom_loss(model(val_epi), val_exp).item()
        val_losses.append(val_loss)

        
        model.train()  

        for batch_x, batch_y in train_data_iter:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = model.custom_loss(outputs, batch_y)
            loss.backward()
            optimizer.step()
        
        train_losses.append(model.custom_loss(model(train_epi), train_exp).item())

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            save_model(model, model_path)
        

    if plot_loss:
        plot_loss_curve(train_losses, val_losses, epochs, cell_type, fig_folder)



def plot_prediction(model, test_features, test_labels, cell_type, fig_folder):
    model.eval()

    predictions = model(test_features).detach()
    if test_labels.is_cuda:  # Check if test_labels is on a GPU
        test_labels = test_labels.to(predictions.device)
    
    predictions_np = predictions.cpu().numpy()
    test_labels_np = test_labels.cpu().numpy()

    correlation, _ = pearsonr(test_labels_np, predictions_np)

    regression_model = LinearRegression()
    regression_model.fit(predictions_np, test_labels_np)

    line_y = regression_model.predict(test_labels_np.reshape(-1, 1))

    plt.figure(figsize=(8, 6))
    plt.scatter(test_labels_np, predictions_np, s=5, label='Data Points', color='skyblue')
    plt.plot(test_labels_np, line_y, color='grey', linestyle='--', label='Best Fit Line')

    plt.xlabel('Observed Expressions')
    plt.ylabel('Predicted Expressions')
    plt.title(f'Model Performance on Test Set \n Pearson correlation: {correlation[0]:.2f}')
    plt.legend()

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.savefig(os.path.join(fig_folder, f'{cell_type}_prediction.pdf'), bbox_inches='tight')
    # plt.show()


def load_model(filepath, model_class):

    model = model_class()
    model.load_state_dict(torch.load(filepath))

    return model
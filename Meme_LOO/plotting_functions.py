import numpy as np
import pandas as pd
from math import log,exp,inf,pi
from sklearn.metrics import mean_squared_error
from statistics import mean,variance
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import subprocess



def plot_distr(gene_ID: tuple, true_sites: list, 
               split_motif: bool, best_n: int, the_n_promt: int, val_scores: dict):

    # plot the distribution of motif scores. 
    # val_scores: a dictrionary generate by Meme_LOO notebook and saved as a pickle 
    # true_sites: the experimental binding sites from Ecocyc
    # split_motif: whether the input dict is the result of split_motif scanning or not
    # best_n: plot red dots on the top n scores
    # the_n_promt: for gene presents in multiple promoters, select which promoter to plot.
    
    fig, ax = plt.subplots();  
    
    if split_motif:
        scores = np.array(val_scores[gene_ID]['D'][the_n_promt])
        scores_inverted = np.array(val_scores[gene_ID]['I'][the_n_promt])
        x=np.tile(np.arange(len(scores)), 1)- 150
        
    else: 
        scores = np.array(val_scores[gene_ID][the_n_promt])
        x=np.tile(np.arange(len(scores)), 1)- 150

    scores_a = scores
    ind = np.argpartition(scores_a, -best_n)[-best_n:]
    ax.plot(x, scores_a)
    
    if split_motif:
        ind_inverted = np.argpartition(scores_inverted, -best_n)[-best_n:]
        ax.plot(x, scores_inverted, alpha= 0.8)
        
    for site in true_sites:
        ax.axvline(site, color='r')

    ax.scatter(x=x[ind], y=scores_a[ind], color='r')
    
    if split_motif:
        ax.scatter(x=x[ind_inverted], y=scores_inverted[ind_inverted], color='r')
        ax.legend(['motif score', 'inverted motif score', 'experimental binding'])
    
    else:
        ax.legend(['motif score', 'experimental binding'])

    ax.set_title('Motif score distribution');
    ax.set_xlabel('location')
    ax.set_ylabel('score')

    ax.set_xlabel('location')
    ax.set_ylabel('score')
    
    return fig
    
def plot_reg_distr(gene_ID: str, true_sites: list, best_n: int, the_n_promt: int, distr_df: pd.DataFrame):
    
    fig, ax = plt.subplots();
    score_dist = []
    for col in list(distr_df.columns):
        if ('_motif_' in col) & ~('im_motif' in col):
            score_dist.append(col)
    
    scores = distr_df[distr_df['ID'] == gene_ID][score_dist].iloc[the_n_promt]
    x=np.tile(np.arange(len(scores)), 1)- 150

    scores_a = np.array(scores).reshape(-1)
    ind = np.argpartition(scores_a, -best_n)[-best_n:]

    ax.plot(x, scores_a)
    for site in true_sites:
        ax.axvline(site, color='r')
    ax.scatter(x=x[ind], y=scores_a[ind], color='r')

    ax.legend(['motif score', '-35box', 'experimental binding'])

    ax.set_title('Motif score distribution');
    ax.set_xlabel('location')
    ax.set_ylabel('score')
    
    return fig
    
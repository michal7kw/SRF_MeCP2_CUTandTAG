### Modules
from collections import Counter
import importlib

import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
from scipy import stats
from scipy.sparse import csr_matrix, isspmatrix

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad

from upsetplot import from_contents, UpSet


def print_dict_summary(data, indent=0, max_list_items=4):
    for key, value in data.items():
        print(' ' * indent + str(key) + ':')
        if isinstance(value, dict):
            print_dict_summary(value, indent + 4, max_list_items)
        elif isinstance(value, list):
            print(' ' * (indent + 4) + '[')
            for i, item in enumerate(value):
                if i >= max_list_items:
                    print(' ' * (indent + 8) + '...')
                    break
                if isinstance(item, dict):
                    print_dict_summary(item, indent + 8, max_list_items)
                else:
                    print(' ' * (indent + 8) + str(item))
            print(' ' * (indent + 4) + ']')
        else:
            print(' ' * (indent + 4) + str(value))

def log_to_raw_counts(adata_in, base=2, log1p=True):
    adata_raw = adata_in.copy()
    X = adata_raw.X
    
    # Convert to dense array if sparse
    if sp.issparse(X):
        X = X.toarray()
    
    # Reverse the log transformation
    if log1p:
        X = np.expm1(X * np.log(base))  # equivalent to base**(X) - 1
    else:
        X = base**X
    
    # Round to nearest integer and convert to integers
    X = np.round(X).astype(int)
    
    # Set negative values to 0 (shouldn't occur, but just in case)
    X[X < 0] = 0
    
    # Update adata_raw.X with the raw counts
    adata_raw.X = X if not sp.issparse(adata_raw.X) else sp.csr_matrix(X)
    
    return adata_raw

### 1. densityQCs

def densityQCs(adataObj, hue=None, 
               min_counts=None, max_counts=None, 
               min_genes=None, max_genes=None,
               pct_counts_mito=None, pct_counts_ribo=None):   
    
    '''
    Returns density plots for the following QCs: 
        n_genes_by_counts (log10), total_counts (log10), pct_counts_mito (linear scale), pct_counts_ribo (linear scale).
    
    Parameters:
    - adataObj: the adata containing the QC metrics above specified. PS: at the moment, if the QC are not present, this raises an error. This point could be improved. 
    - hue: xxx. Default is None. 
    - min_counts
    - ** kwargs: I used it to be flexible in passing anything we want to filter on, arguments should be either:
        - valid arguments for sc.pp-filter_cells()
        - column names from adata.obs (and pass as corresponding value the label(s, list or tuple) that you want to filter)
        '''
    
    #Plot them in line so they take up less space
    fig, ax = plt.subplots(1, 4, figsize=(20,5))
    fig.tight_layout(pad=2)   #space between plots
    if hue != None:
        hue_s = adata.obs[hue].astype('string')
    else:
        hue_s = None

    ### Genes ---------------
    n_genes = adataObj.obs['n_genes_by_counts']
    n_genes_log10 = np.log10(n_genes.replace(0, np.nan))  # Replace 0 with NaN before log10
    d1 = sns.kdeplot(n_genes_log10.dropna(), fill=True, color='cornflowerblue', hue=hue_s, ax=ax[0])
    min_x, max_x = d1.get_xlim() 

    #Threshold lines and fill
    if min_genes != None:
        d1.axvline(np.log10(min_genes), 0, 1, c='red')  #set manually for chosen threshold
        d1.axvspan(min_x, np.log10(min_genes), alpha=0.2, color='red')
    if max_genes != None:
        d1.axvline(np.log10(max_genes), c='red')
        d1.axvspan(np.log10(max_genes), max_x, alpha=0.2, color='red')

    ### UMI ---------------
    total_counts = adataObj.obs['total_counts']
    total_counts_log10 = np.log10(total_counts.replace(0, np.nan))  # Replace 0 with NaN before log10
    d2 = sns.kdeplot(total_counts_log10.dropna(), fill=True, color='forestgreen', hue=hue_s, ax=ax[1])
    min_x, max_x = d2.get_xlim() 
        
    if min_counts != None:
        d2.axvline(np.log10(min_counts), 0, 1, c='red')  #set manually for chosen threshold
        d2.axvspan(min_x, np.log10(min_counts), alpha=0.2, color='red')
    if max_counts != None:
        d2.axvline(np.log10(max_counts), c='red')
        d2.axvspan(np.log10(max_counts), max_x, alpha=0.2, color='red')

    ### Mito % ---------------
    mito_data = adataObj.obs['pct_counts_mito']
    if mito_data.var() > 0:  # Check if there's variance in the data
        d3 = sns.kdeplot(mito_data, fill=True, color='coral', hue=hue_s, ax=ax[2])
        min_x, max_x = d3.get_xlim() 

        #Threshold lines and fill
        if pct_counts_mito != None:
            d3.axvline(pct_counts_mito, 0, 1, c='red')  #set manually for chosen threshold
            d3.axvspan(pct_counts_mito, max_x, alpha=0.2, color='red')
    else:
        ax[2].text(0.5, 0.5, 'No variance in mito data', ha='center', va='center')

    ### Ribo % ---------------
    ribo_data = adataObj.obs['pct_counts_ribo']
    if ribo_data.var() > 0:  # Check if there's variance in the data
        d4 = sns.kdeplot(ribo_data, fill=True, color='orchid', hue=hue_s, ax=ax[3])
        min_x, max_x = d4.get_xlim() 
        
        #Threshold lines and fill
        if pct_counts_ribo != None:
            d4.axvline(pct_counts_ribo, 0, 1, c='red')  #set manually for chosen threshold
            d4.axvspan(pct_counts_ribo, max_x, alpha=0.2, color='red')
    else:
        ax[3].text(0.5, 0.5, 'No variance in ribo data', ha='center', va='center')
    
    #Remove additional legends at need
    if hue != None:
        ax[0].get_legend().remove()
        ax[1].get_legend().remove()
        ax[2].get_legend().remove()
        
    # Remove all borders
    sns.despine(bottom = False, left = True)
 
### 2. filterCellBarplot
   
def filterCellBarplot(adataObj, **kwargs):
    
    '''
    Returns a barplot depicting cells present at the start of each spcified filtering step and how many cells get filtered.
    
    Parameters:
    - adataObj: the adata object to filter. The object original is NOT filtered (a copy is made), only the plot is returned.
    - ** kwargs: I used it to be flexible in passing anything we want to filter on, arguments should be either:
        - valid arguments for sc.pp-filter_cells()
        - column names from adata.obs (and pass as corresponding value the label(s, list or tuple) that you want to filter)
        '''
    
    sc.settings.verbosity = 1  #Don't show filtering messages
    
    #INIZIALIZE
    adata_filt = adataObj.copy()
    #Store the number of cells at each filtering step
    n_cells = [adata_filt.n_obs]
    #Compute how many cells are removed at each step
    removed = []
    #Store names for plot labels
    names = []   
    
    #FILTERING
    #I used kwargs and no named arguments to run everything in a for loop, maybe not too elegant
    #for each key-value preform filtering and append cell number and labels
    
    for key, value in kwargs.items(): 
        #Check exact arguments for sc.pp.filter_cells
        if key in ['min_counts', 'max_counts', 'min_genes', 'max_genes']:  
            sc.pp.filter_cells(adata_filt, **{key: value})
            n_cells.append(adata_filt.n_obs)
            names.append(key)    
        
        else:
            assert key in adataObj.obs.columns, 'Please specify valid adata.obs column names as arguments'
        
            #MORE FLEXIBLE: check for any argument containing mito/ribo
            #elif key in ['pct_counts_mito','pct_counts_ribo']:
            if 'pct_counts_mito' in key or 'pct_counts_ribo' in key:  
                 adata_filt = adata_filt[adata_filt.obs[key] < value, :]
                 n_cells.append(adata_filt.n_obs)
                 names.append(key)   
                
            else:
                #When we filter on more than one label in the same .obs column
                if type(value) is list or type(value) is tuple:
                    for el in value:
                        assert el in adataObj.obs[key].unique(), f"Please specify valid adata.obs['{key}'] values as parameter"
                        
                        adata_filt = adata_filt[adata_filt.obs[key] != el, :]
                        n_cells.append(adata_filt.n_obs)
                        names.append(el)   
                else:
                    assert value in adataObj.obs[key].unique(), f"Please specify valid adata.obs['{key}'] value as parameter"
                    
                    adata_filt = adata_filt[adata_filt.obs[key] != value, :]
                    n_cells.append(adata_filt.n_obs)
                    names.append(value) 
        
        
    #Update lists for plot
    names.append('end')
    removed = n_cells[1:] + [adata_filt.n_obs]
    #Use numpy arrays to subtract vectors
    rem = np.array(removed) - np.array(n_cells) #this order to have negative values
    
    #PLOT
    plt.figure(figsize=(20,8))
    plot1 = sns.barplot(x=names, y=n_cells, color='turquoise',  label = "starting cells")
    plot2 = sns.barplot(x=names, y=rem, color='crimson',  hatch='/',  label = "filtered cells").set(xlabel='Filtering Step', ylabel='Number of cells')
    plt.axhline(0, 0, 1, c='black') 
    plt.legend(frameon = False)
    
    # Annotate the bars
    for p in plot1.patches[:len(names)]:
        plot1.annotate(format(p.get_height(), '.0f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, 10), 
                       textcoords = 'offset points') 
        
    for p in plot1.patches[len(names):]:
        plot1.annotate(format(p.get_height(), '.0f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, -10), 
                       textcoords = 'offset points') 
    
    sc.settings.verbosity = 3
    plt.show()    

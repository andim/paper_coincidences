import numpy as np
import pandas as pd
import scipy.sparse

from pyrepseq import *

df = pd.read_csv(snakemake.input[0], sep='\t') 

bins = np.arange(35)
distances_alpha = pdist(df['CDRL3_AA'])
distances_beta = pdist(df['CDRH3_AA'])
distances_pair = distances_alpha + distances_beta

#### save raw distance matrix ####
mask = squareform(distances_pair) < 5
# exclude diagonal
np.fill_diagonal(mask, False)
for chain, distances in [('alpha', distances_alpha),
                         ('beta', distances_beta)]:
    cdist_sparse = scipy.sparse.coo_matrix((squareform(distances)[mask], mask.nonzero()))
    # filter out correct output file
    out = [item for item in snakemake.output if chain in item][0]
    scipy.sparse.save_npz(out,
        cdist_sparse, compressed=True)

#### save distance histograms ####
res_dict = dict(bins=bins[:-1])
for distances, chain in [(distances_alpha, 'alpha'),
                         (distances_beta, 'beta'),
                         (distances_pair, 'pair')]:
    hist, bins_ = np.histogram(distances, bins=bins)
    norm = hist/np.sum(hist)
    res_dict['count_'+chain] = hist
    res_dict['freq_'+chain] = norm
out = pd.DataFrame.from_dict(res_dict)
out.to_csv(snakemake.output[0], index=False)

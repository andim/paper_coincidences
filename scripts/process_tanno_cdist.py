import numpy as np
import pandas as pd
import scipy.sparse

from pyrepseq import *

dfs = [pd.read_csv(filein, sep='\t') for filein in snakemake.input]

bins = np.arange(35)

distances_alpha = cdist(dfs[0]['CDRL3_AA'], dfs[1]['CDRL3_AA'])
distances_beta = cdist(dfs[0]['CDRH3_AA'], dfs[1]['CDRH3_AA'])
distances_pair = distances_alpha + distances_beta

# filter identical alpha or beta nucleotide sequences
masks = []
for chain in ['H', 'L']:
    A = dfs[0]['CDR{chain}3_NT'.format(chain=chain)]
    B = dfs[1]['CDR{chain}3_NT'.format(chain=chain)]
    mask_nt_chain = np.array(A.isin(B))[:, np.newaxis] * np.array(B.isin(A))[np.newaxis, :]
    masks.append(mask_nt_chain)
mask_nt = masks[0] | masks[1]

#### save raw distance matrix ####

# threshold at pair distance below 5 to save storage space
mask = (distances_pair < 5) & (~mask_nt)
for i, (chain, distances) in enumerate([('alpha', distances_alpha),
                         ('beta', distances_beta)]):
    cdist_sparse = scipy.sparse.coo_matrix((distances[mask], mask.nonzero()))
    scipy.sparse.save_npz(snakemake.output[i+1], cdist_sparse, compressed=True)

#### save distance histograms ####
res_dict = dict(bins=bins[:-1])
for distances, chain in [(distances_alpha, 'alpha'),
                         (distances_beta, 'beta'),
                         (distances_pair, 'pair')]:
    hist, bins_ = np.histogram(distances[~mask_nt], bins=bins)
    norm = hist/np.sum(hist)
    res_dict['count_'+chain] = hist
    res_dict['freq_'+chain] = norm
out = pd.DataFrame.from_dict(res_dict)
out.to_csv(snakemake.output[0], index=False)

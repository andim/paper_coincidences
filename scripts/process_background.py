import itertools
import numpy as np
import pandas as pd
import scipy.stats
import pyrepseq as prs

data_directory ='data/'

bins = np.arange(0, 40, 1)

chains = ['alpha', 'beta']
back_norms = {}
for chain in chains:
    df = pd.read_csv(data_directory + 'minervina/{chain}/W_F1_2018_{chain}.txt.gz'.format(chain=chain), sep='\t')
    df.dropna(subset=['aaSeqCDR3'], inplace=True)
    df = df[df['aaSeqCDR3'].apply(prs.isvalidcdr3)]
    df = df.sample(n=30000, replace=False)
    back_norms[chain] = prs.pcDelta(df['aaSeqCDR3'])
    
fa = back_norms['alpha']
fb = back_norms['beta']
fab_ind = np.convolve(fa, fb, mode='full')[:fa.shape[0]]
back_norms['both'] = fab_ind

df = pd.DataFrame(back_norms)
df.to_csv(data_directory + 'background_minervina.csv')

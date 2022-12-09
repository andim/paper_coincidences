import numpy as np
import pandas as pd

import pyrepseq as prs

df = pd.read_csv(snakemake.input[0], sep='\t')
# sort by clone sizes
df = df.sort_values('Clustered', ascending=False)

valid_light = df['CDRL3_AA'].apply(prs.isvalidcdr3)
valid_heavy = df['CDRH3_AA'].apply(prs.isvalidcdr3)
before = len(df)
df = df[valid_light & valid_heavy]
print(before, np.sum(~valid_light), np.sum(~valid_heavy), len(df))

# what prunings will we do?
invariants = True
memory_naive = True
unique_nt = True

if invariants:
    # now remove the invariant cells, defined by the following alpha chain V/J pairs
    invariant = [('TRAV1-2', 'TRAJ33'), 
                 ('TRAV1-2', 'TRAJ12'),
                 ('TRAV1-2', 'TRAJ20'),
                 ('TRAV10', 'TRAJ18')]
    invariant_joined = [s1+'_'+s2 for s1,s2 in invariant]

    mask = (df['VL'] + '_' + df['JL']).isin(invariant_joined)
    df = df[~mask]
    print('invariant filtered', sum(mask))

if memory_naive:
    if 'naive' in snakemake.input[0]:
        dfmem = pd.read_csv(snakemake.input[0].replace('naive', 'memory'), sep='\t')
        A = df['CDRH3_AA'] + '_' + df['CDRL3_AA']
        B = dfmem['CDRH3_AA'] + '_' + dfmem['CDRL3_AA']
        mask = A.isin(B)
        df = df[~mask]
        print('memory overlap filtered', sum(mask))

if unique_nt:
    before = len(df)
    df = df.drop_duplicates('CDRH3_NT', keep='first')
    df = df.drop_duplicates('CDRL3_NT', keep='first')
    print('nt duplicates filtered', before-len(df))


df.to_csv(snakemake.output[0], index=False, sep='\t')

import pandas as pd

df = pd.read_csv(snakemake.input[0]) 
columns = df.columns.tolist()
df['amino_acid'] = df['TCR BioIdentity'].apply(lambda s: s.split('+')[0])
df['v_gene'] = df['TCR BioIdentity'].apply(lambda s: s.split('+')[1])
df['j_gene'] = df['TCR BioIdentity'].apply(lambda s: s.split('+')[2])
reordered_columns = ['amino_acid', 'v_gene', 'j_gene']
reordered_columns.extend(columns)
df = df[reordered_columns]
df.to_csv(snakemake.output[0], index=False)

import itertools

# number of sequences to generate
nseqs = 100000

tannofiles = ['D1', 'D2', 'E1', 'E2', 'F1', 'F2',
              'X', 'Y', 'Z',
              'A1 memory', 'A1 naive',
              'A2 memory', 'A2 naive',
              'B1 memory', 'B1 naive',
              'B2 memory', 'B2 naive',
              'C1 memory', 'C1 naive',
              'C2 memory', 'C2 naive',
              'X replicate',
              'Y replicate',
              'Z replicate',
              'D1 replicate',
              'D2 replicate',
              'E1 replicate',
              'E2 replicate',
              ]

tannofiles_sorted = [
              'A1 memory', 'A1 naive',
              'A2 memory', 'A2 naive',
              'B1 memory', 'B1 naive',
              'B2 memory', 'B2 naive',
              'C1 memory', 'C1 naive',
              'C2 memory', 'C2 naive',
              ]

chains = ['alpha', 'beta']

localrules:
    all,
    download_mira,
    download_dash,
    download_minervina_repertoire,
    download_minervina_dynamic_clones,
    download_minervina2022,
    download_nguyen,
    download_paired_tanno,

rule all:
    input:
        'data/post_sample_alpha.csv',
        expand('data/tanno/{name}.txt', name=tannofiles),
        expand('data/tanno/pruned/pdist_{name}.csv', name=tannofiles),
        'data/mira/processed/peptide-detail-ci.csv',
        'data/mira/processed/peptide-detail-cii.csv',
#        expand('data/tanno/pruned/pdist_full_{chain}_{name}.npz', chain=chains, name=tannofiles),
        ['data/tanno/pruned/cdist_{name1}_{name2}.csv'.format(name1=name1, name2=name2)
            for name1, name2 in itertools.combinations(tannofiles, 2)],
#        ['data/tanno/pruned/cdist_full_{chain}_{name1}_{name2}.npz'.format(chain=chain, name1=name1, name2=name2)
#            for chain in chains
#            for name1, name2 in itertools.combinations(tannofiles, 2)],
        'data/dash_human.csv',
        'data/nguyen_sarscov2.xlsx',
        'data/bacher.csv',
        'data/minervina2022.csv',
        'data/minervina2022.xlsx',
        expand('data/minervina/Table_S{table}.tsv', table=(3,4,5,6)),
        'data/background_minervina.csv',

rule download_mira:
    output:
        'data/mira/peptide-detail-ci.csv',
        'data/mira/peptide-detail-cii.csv',
    shell:
        'wget https://adaptivepublic.blob.core.windows.net/publishedproject-supplements/covid-2020/ImmuneCODE-MIRA-Release002.1.zip -O data/temp.zip\n'
        'unzip -j data/temp.zip -d data/mira/\n'
        'rm data/temp.zip'

rule download_paired_tanno:
    output:
        expand('data/tanno/{name}.txt', name=tannofiles),
    shell:
        'wget DATALINK -O data/temp.zip\n'
        'unzip -j data/temp.zip -d data/tanno/\n'
        'rm data/temp.zip'

rule download_dash:
    output:
        'data/dash_human.csv',
    shell:
        'wget https://raw.githubusercontent.com/kmayerb/tcrdist3/master/dash_human.csv -O {output}'

rule download_nguyen:
    output:
        'data/nguyen_sarscov2.xlsx',
    shell:
        'wget https://github.com/antigenomics/vdjdb-db/files/6683583/SARS-CoV-2.TCR.sequences.xlsx -O {output}'

rule download_bacher:
    output:
        'data/bacher.csv',
    shell:
        'wget https://www.cell.com/cms/10.1016/j.immuni.2020.11.016/attachment/1e7636f0-6036-40df-b7eb-aafd4b794ebe/mmc2.xlsx -O data/bacher.xlsx\n'
        'xlsx2csv data/bacher.xlsx -s 2 > data/bacher.csv\n'
        'rm data/bacher.xlsx'

rule download_minervina_repertoire:
    output:
        'data/minervina/{chain}/W_F1_2018_{chain}.txt.gz',
    wildcard_constraints: chain='(\w+)',
    shell: 
        'mkdir -p data/minervina/{wildcards.chain}\n'
        'wget https://zenodo.org/record/4065547/files/{wildcards.chain}.zip?download=1 -O data/temp.zip\n'
        'unzip -j data/temp.zip -d data/minervina/{wildcards.chain}\n'
        'rm data/temp.zip'

rule download_minervina_dynamic_clones:
    output:
        'data/minervina/Table_S{table}.tsv',
    shell:
        'wget https://raw.githubusercontent.com/pogorely/Minervina_COVID/master/Table_S{wildcards.table}.tsv -O {output}\n'

rule download_minervina2022:
    output:
        'data/minervina2022.csv',
        'data/minervina2022.xlsx'
    shell:
        'wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-022-01184-4/MediaObjects/41590_2022_1184_MOESM7_ESM.xlsx -O data/minervina2022.xlsx\n'
        'xlsx2csv data/minervina2022.xlsx -s 2 > data/minervina2022.csv\n'

rule process_mira:
   input: 
        'data/mira/peptide-detail-c{mhc}.csv',
   wildcard_constraints: mhc='i{1,2}',
   output:
        'data/mira/processed/peptide-detail-c{mhc}.csv',
   script:
        'scripts/process_bioidentity.py'

rule process_tanno_pruning:
    input:
        'data/tanno/{name}.txt',
    output:
        'data/tanno/{pruned}/{name}.txt',
    script:
        'scripts/process_tanno_pruning.py'

rule process_tanno_pdist:
    input:
        'data/{path}/{name}.txt',
    wildcard_constraints:
        name='|'.join(tannofiles),
        path='|'.join(['tanno/pruned', 'tanno']),
    resources:
        mem_mb=32000
    output:
        'data/{path}/pdist_{name}.csv',
        expand('data/{{path}}/pdist_full_{chain}_{{name}}.npz', chain=chains)
    script:
        'scripts/process_tanno_pdist.py'

rule process_tanno_cdist:
    input:
        'data/{path}/{name1}.txt',
        'data/{path}/{name2}.txt',
    wildcard_constraints:
        path='|'.join(['tanno/pruned', 'tanno']),
        name1='|'.join(tannofiles),
        name2='|'.join(tannofiles),
    resources:
        mem_mb=32000
    output:
        'data/{path}/cdist_{name1}_{name2}.csv',
        expand('data/{{path}}/cdist_full_{chain}_{{name1}}_{{name2}}.npz', chain=chains)
    script:
        'scripts/process_tanno_cdist.py'

rule process_background:
    input:
        expand('data/minervina/{chain}/W_F1_2018_{chain}.txt.gz', chain=chains)
    output:
        'data/background_minervina.csv'
    script:
        'scripts/process_background.py'

rule post_sample_alpha:
    output:
        'data/post_sample_alpha.csv',
    shell:
        'sonia-generate --humanTRA --post -o {output} -n {nseqs} -r 100'


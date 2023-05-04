##CoVAbDav_database
import pandas as pd
import numpy as np

df = pd.read_csv(r'CoV-AbDab_database_file.csv')
df = df.loc[df['Ab or Nb']=='Ab'] # select only rows containing Ab data

##new columns that split 'protein' and 'epitope' etc: S,RBD, select only spike specific antibodies
df[['protein','epitope']]=df['Protein + Epitope'].apply(lambda x: pd.Series(str(x).split(";"))) 
df = df.loc[df['protein']=='S'] 

## rename columns to match tcrdist requirements
df.rename(columns={'Name':'clone_id', 'Binds to':'epitope', 'Heavy V Gene':'v_b_gene', 'Heavy J Gene':'j_b_gene', 'Light V Gene':'v_a_gene', 'Light J Gene':'j_a_gene', 'CDRH3':'cdr3_b_aa','CDRL3':'cdr3_a_aa', 'epitope':'protein_epitope'}, inplace=True)

## split gene from species eg. IGHJ6 (Human)
df[['v_b_gene','vb_species']]=df['v_b_gene'].apply(lambda x: pd.Series(str(x).split()))
df[['j_b_gene','jb_species']]=df['j_b_gene'].apply(lambda x: pd.Series(str(x).split()))
df[['v_a_gene','va_species']]=df['v_a_gene'].apply(lambda x: pd.Series(str(x).split()))
df[['j_a_gene','ja_species']]=df['j_a_gene'].apply(lambda x: pd.Series(str(x).split()))

## select only human data
df = df.loc[df['vb_species']=='(Human)']
df = df.loc[df['va_species']=='(Human)']
df = df.loc[df['jb_species']=='(Human)']
df = df.loc[df['jb_species']=='(Human)']

df.drop(['Ab or Nb',"Doesn't Bind to", 'Not Neutralising Vs','Protein + Epitope', 'Origin', 'VHorVHH', 'VL','Structures', 'ABB Homology Model (if no structure)', 'Date Added', 'Last Updated', 'Update Description', 'Notes/Following Up?', 'protein', 'vb_species', 'jb_species', 'va_species', 'ja_species'], inplace=True, axis=1)
df['subject']='CoVAbDab_db'
## We will also assign a ‘count’ number of 1 to each clonotype in the table.
df['count']=1

## Because tcrdist3 infers CDR1, CDR2, and CDR2.5 from the full gene and allele name
## we must append an allele designation (*01 in this case) when it is unresolved in the input data. 
if '*' not in df['v_a_gene']:
    df['v_a_gene'] += '*01'
if '*' not in df['v_b_gene']:
    df['v_b_gene'] += '*01'
if '*' not in df['j_b_gene']:
    df['j_b_gene'] += '*01'
if '*' not in df['j_a_gene']:
    df['j_a_gene'] += '*01'    

for column in df:
    df[column]=df[column].str.replace(' ','') ##drop white space

df.to_csv(r'CoVAbDab_database.csv', index=False)    
## antigen2barcode

def antigen2barcode(file): 
        import csv
        import gzip
        import os
        import scipy.io
        import pandas as pd

        # define MEX directory
        matrix_dir = os.path.join(file, "sample_filtered_feature_bc_matrix/")
        # read in MEX format matrix as table
        mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))

        # list of transcript ids, e.g. 'ENSG00000243485'
        features_path = os.path.join(matrix_dir, "features.tsv.gz")
        feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

        # list of gene names, e.g. 'MIR1302-2HG'
        feature_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

        # list of feature_types, e.g. 'Gene Expression'
        barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
        barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

        
        # transform table to pandas dataframe and label rows and columns
        matrix = pd.DataFrame.sparse.from_spmatrix(mat)
        matrix.columns=barcodes
        matrix.insert(loc=0, column="feature", value=feature_names)

        # Transverse matrix for features along the top
        df=matrix.T
        df.columns=df.iloc[0] #rename columns to features (currently in row 1 -- index0)
        df.drop(labels='feature', axis=0, inplace=True) #drop unneeded row 1

        #get klickmer info
        klickmers=pd.read_csv(os.path.join(file, "klickmers.csv"))                                                              #(r'C:\Users\BFCRa\OneDrive\Documents\p604\klickmers.csv') #use r for windows files
        klickmers_dict=klickmers.set_index('barcode')['antigen'].to_dict()
        dict_list=list(klickmers_dict.keys())

        # Drop all features columns except covid-specific 
        # df.drop(df.columns.difference(['fBCO604','fBCO605','fBCO607','fBCO213','fBCO214']),axis=1, inplace=True)
        df.drop(df.columns.difference(dict_list),axis=1, inplace=True) 
        df.reset_index(names='barcode',inplace=True)
        df.rename(columns=klickmers_dict, inplace=True)

        epitopes = {} # initialize epitopes dictionary

        for barcode in df.index:
                for i in df.columns[1:]:
                        if df.at[barcode,i] >=1:
                                epitope = i
                                key=df.at[barcode,'barcode']
                                if key in epitopes:
                                        epitopes[key]= epitopes[key]+','+epitope
                                else:
                                        epitopes[key] = epitope


        epitope_df=pd.DataFrame.from_dict(epitopes, orient='index')
        epitope_df=epitope_df.reset_index()
        epitope_df.rename(columns={'index':'barcode', 0:'epitope'}, inplace=True)

        df_full = pd.merge(df,epitope_df,how='left', on='barcode') ## combine epitope matrix with epitope list (for full list of barcodes)
        #df_full.fillna('Unk', inplace=True) ## indicate barcodes with missing epitopes (NaN) are unknown (Unk) epitopes
        df_full.to_csv(os.path.join(file,"barcode_epitopes.csv"),index=False)
        
        #merge with fca
        fca_df=pd.read_csv(os.path.join(file, "filtered_contig_annotations.csv"))
        df_final = pd.merge(fca_df, df_full, how='left', on='barcode')
        df_final['epitope'].fillna('Unk', inplace=True)
        df_final.to_csv(os.path.join(file,"filtered_contig_annotations_epitopes.csv"),index=False)

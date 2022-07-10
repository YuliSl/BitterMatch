import numpy as np
import pandas as pd

def load_A(A_file):
    A = pd.read_csv(A_file, sep=',', index_col=0)
    A.columns = A.columns.astype(int)
    return A

def load_X_Rec(X_Rec_file):
    X_Rec = pd.read_csv(X_Rec_file, sep=',')  
    X_Rec.Organism = pd.get_dummies(X_Rec.Organism)
    X_Rec.Organism = pd.get_dummies(X_Rec.Organism)
    return X_Rec

def load_X_Lig(X_Lig_file):
    X_Lig = pd.read_csv(X_Lig_file, sep=',') 
    # remove duplicate rows by cid
    _, ixs = np.unique(X_Lig.cid, return_index=True)  # indices of unique cid
    X_Lig = X_Lig.iloc[ixs]
    X_Lig = X_Lig.reset_index(drop=True)
    # remove SMILES column
    X_Lig.drop(['SMILES' ], inplace=True , axis=1)
    if 'Unnamed: 0' in X_Lig.columns:
        X_Lig.drop('Unnamed: 0', inplace=True , axis=1)

    return X_Lig

def family_features(families_file):
    families_df = pd.read_csv(families_file, sep=',')
    # include 4 most general structural classifications to chemical families for each ligand
    families_df['family'] = families_df.groupby('CompoundID').cumcount(ascending=False)
    families_df = families_df[families_df.family < 4]
    families_df.chemotype = pd.Categorical(families_df.chemotype).codes
    families_df.family = families_df.family.apply(lambda f: 'family%s' % (f+1,))
    families_df = families_df.pivot(index='CompoundID', columns='family', values='chemotype')
    families_df.index.name = None
    families_df.columns.name = None
    families_df = families_df.replace(np.nan, -1)
    families_df = families_df[:].astype(int)
    return families_df

def read_receptor_similarity(data, from_file=True):
    if from_file:
        df = pd.read_csv(data, sep=',', index_col=0)
    else: 
        df = data
    df = pd.read_csv(data, sep=',', index_col=0, header=None, skipfooter=1, engine='python')
    df[:] = np.tril(df) + np.tril(df, -1).T  # complete the upper triangle
    df.columns = df.index
    df.index.name = None
    return df

def read_ligand_similarity(data, from_file=True):
    if from_file:
        df = pd.read_csv(data, sep=',', index_col=0)
    else: 
        df = data
    _, ixs = np.unique(df.index, return_index=True)  # indices of unique rows
    df = df.iloc[ixs, ixs]
    df.columns = df.index
    df.index.name = None
    return df
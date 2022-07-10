import numpy as np
import pandas as pd

def collaborative_similarity(A, similarity_type):
    """
    :param A: a pandas data frame with columns: ligand (integer id), receptr (integer id), association - either 0.0 or 1.0 (np.float64) 
    :param p: train proportion (float)
    :param similarity_type: either 'Ligand' or 'Receptor' (string)
    :return: similarity matrix S 
    """
    if similarity_type=='Ligand':
        S = pd.DataFrame(index=A.index, columns=A.index)
        A_sum = A.values[:, None, :] + A.values[None, :, :]
        axs = 2
    if similarity_type=='Receptor':
        S = pd.DataFrame(index=A.columns, columns=A.columns)
        A_sum = A.values[:, :, None] + A.values[:, None, :]
        axs = 0
    S[:] = (np.sum(A_sum == 2, axis=axs) + np.sum(A_sum == 0, axis=axs)) / (np.sum(A_sum == 2, axis=axs) + np.sum(A_sum == 0, axis=axs) + np.sum(A_sum == 1, axis=axs))
    return S

def sim_metrics(S, A, axis):
    if axis=='row' or axis==0:
        rows = A.index.intersection(S.index)
        cols = A.columns
        A_vals = A.loc[rows, :].values.copy()
        S_vals = S.loc[rows, rows].values.copy()
        transpose = False
    elif axis=='col' or axis==1:
        rows = A.index
        cols = A.columns.intersection(S.index)
        A_vals = A.loc[:, cols].values.copy().T
        S_vals = S.loc[cols, cols].values.copy()
        transpose = True
    else:
        raise ValueError('axis must be either "row" or 0, or "col" or 1.')
    
    np.fill_diagonal(S_vals, 0)

    # weighted similarities
    W1 = S_vals.dot(np.nan_to_num(A_vals)) 
    # if axis = 'row': a matrix A of size L*R in which the value A_lr is the sum of all similarities to l of ligands that associate with r.
    # if axis = 'col': a matrix A of size R*L in which the value A_rl is the sum of all similarities to r of receptors that associate with l.
    
    W0 = S_vals.dot(np.nan_to_num(1 - A_vals)) 
    # if axis = 'row': a matrix A of size L*R in which the value A_lr is the sum of all similarities to l of ligands that do not associate with r.
    # if axis = 'col': a matrix A of size R*L in which the value A_rl is the sum of all similarities to r of receptors that do not associate with l.
    
    # Nearest neighboors
    M01 = np.array([(np.max(S_vals*(line==0), axis=1), #nn
                       np.max(S_vals*(line==1), axis=1)) for line in A_vals.T]).T
    # line is the opposite of axis
    # if axis = 'row': a 3 dimensional tensor M with size R*L*R in which the value M_l0r is the similarity of l to the closest ligand that does not associate with r. 

    # regardless of axis, transform to size L*R (or R*L*R)
    if transpose:
        W1 = W1.T
        W0 = W0.T
        M01 = M01.T
    
    vals = np.vstack([W1.flatten(), W0.flatten(),
                      M01[:, 1, :].flatten(), M01[:, 0, :].flatten()]).T
    res = pd.DataFrame(vals,
                          index=pd.MultiIndex.from_product([rows, cols], names=['ligand', 'receptor']),
                          columns=['W0', 'W1', 'M1', 'M0'])
    return res
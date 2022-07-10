# BitterMatch

This repository accompanies the paper
[BitterMatch: recommendation systems for matching molecules with bitter taste receptors](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00612-9).

Updated versions of the data can be found at [BitterDB](https://bitterdb.agri.huji.ac.il/dbbitter.php).
Please cite both papers when using this repository.

#### Requirements
The file `requirements.txt` lists the package requirements. 

For convenience the notebooks were adapted for running also in Google Colab. 
_______________________________________________________

##### Filling the Gaps Scenario
- `filling_the_gaps-train.ipynb` - On a single train-test split, the notebook trains all BitterMatch models for filling the gaps, and demontrsates the results.
- `filling_the_gaps-eval.ipynb` - Using a pre-trained model the notebook allows to predict activations for unknown values in the association matrix.

##### New Ligands Scenario
- `new_ligands-train.ipynb` - On a single train-test split, the notebook trains the BitterMatch model for new ligands and demonstrates the results.
- `filling_the_gaps-eval.ipynb` - Using a pre-trained model the notebook loads data for ligands that were not used at training (evaluation data) and predicts activations for them.

All notebooks use the file `similarity.py`, that includes the functions to calculate collaborative similarities and extract similarity based features.
The file `preprocessing.py` includes helper functions to load the data from  external formats as detailed in the paper.


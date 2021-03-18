# FDR-validation
A tool to assess the quality of the atomic models derived from the cryoEM maps based on the False Rate Discovery (FDR) analysis


usage: 

fdr_val.py [-h] [-input_pdb INPUT_PDB] [-input_map INPUT_MAP] [-prune]
                  [-valCA]

optional arguments:
-h, --help            show this help message and exit
-input_pdb INPUT_PDB  Input pdb file
-input_map INPUT_MAP  FDR map to validate the model
-prune, --prune       pruning mode
-valCA, --valCA       validate by CA position

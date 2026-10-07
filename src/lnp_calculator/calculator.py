#%%
import json
import pandas as pd
# %%
with open("./components.json", 'r') as f:
    components = json.load(f)

with open("./formulation_vals.json") as f:
    formulation_vals = json.load(f)
# %%

# Constants
MRNA_LENGTH = 4000
TOTAL_MRNA = 5000
NP_RATIO = 6
MASTER_MIX =1.6
MRNA_CONC = 1000
VOL_RATIO = 3
MRNA_CONC_IN_LNP = 112.5


#%%
def p_ratio(
    
)
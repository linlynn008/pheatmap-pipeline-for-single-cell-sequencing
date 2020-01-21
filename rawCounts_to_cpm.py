# cpm: counts per million

import numpy as np
import pandas as pd 


test_data = pd.read_csv("Downloads/seurat_pbmc_data_V2.csv")

#test_data.iloc[:, 1].sum()
#test_data.iloc[21,1]
for col in test_data.columns[1:]:
    test_data[col] = np.log2(test_data[col] * (10**6) / test_data[col].sum() + 1)

test_data.to_csv("Downloads/seurat_cpm.csv")

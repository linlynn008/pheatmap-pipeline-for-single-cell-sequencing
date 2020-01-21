# cpm:counts per million

import numpy as np
import pandas as pd 

#test_data = [[0,3,1],[1,0,7], [0,1,1], [3,2,0]]

test_data = pd.read_csv("C:/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19-1748-51cells/seurat_pbmc_data_V2.csv")
test_data.head()
test_data_np = test_data.to_numpy()
test_data_npn = test_data_np[:,1:]
#test_data_npn[:5,]

l = []
m = len(test_data_npn)
n = len(test_data_npn[0])
for xi in range(0,n):
    sum = 0
    for yi in range(0,m):
        sum = sum + test_data_npn[yi][xi]
    l.append(sum)
print(l)


df = []
for xi in range(0, n):
    for yi in range(0,m):
        cpm = np.log2((test_data_npn[yi][xi] * 10 ** 6 / l[xi]) + 1)       
        df.append(cpm)
#len(df)


d = []
for i in range(0, n):
    s = df[i*m:(i+1)*m]
    d.append(s)
df_np = np.transpose(d)
#len(df_np)
#len(df_np[0])

col_names = test_data.columns.values[1:]
row_names = test_data_np[:, 0]

df_final = pd.DataFrame(df_np, index = row_names, columns = col_names)

df_final.head()



df_final.to_csv(r"C:/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19-1748-51cells/seurat_cpm.csv")





    



  

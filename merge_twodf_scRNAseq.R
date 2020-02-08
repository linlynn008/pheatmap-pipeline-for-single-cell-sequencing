df1 = read.csv("R:/imPACT Dex Experiments/2_Experiments/2019_EXPERIMENT_FILES/Study-18-185_Smarter_Seq/EXP19002657_PACT264_CD62L_phenotype_Smart103019/DeepSequencing/seurat_pbmc_data.csv")

df2 = read.csv("R:/imPACT Dex Experiments/2_Experiments/2019_EXPERIMENT_FILES/Study-18-185_Smarter_Seq/EXP19002657_PACT264_CD62L_phenotype_Smart103019/DeepSequencing/from Ben/seurat_pbmc_data_with_search.csv")
df1_2 = merge(x = df1,y =df2, by = "X", all=TRUE);
col_names_df = colnames(df1_2)
m = length(col_names_df)
col_names_df = col_names_df[2:m]

sample.names = ''
for (s in col_names_df){
  i = strsplit(s, "_")[[1]][1]
  sample.names = c(sample.names, i)
}
sample.names = sample.names[2:m]

sample.names = c('gene_names', sample.names)
colnames(df1_2) = sample.names 
write.csv(df1_2, 'seurat_pbmc_data_v2merged.csv' )

transform_single <- function(inp, outp){
  df0 <- read.delim(inp, header=FALSE, col.names = c('name', 'unstranded_readcounts', 'counts_1st_read', 'counts_2nd_read'))
  df1<-tail(df0,-4)
  df2 <- merge(gn, df1, by.x="query", by.y="name")
  write.csv(df2, file =outp)
}
transform_single(a, b)

input_root <- "c/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19_1502/gene_counts"
output_root <- "c/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19_1502/gene_counts_csv_files"
for (f in dir(input_root)) {
  inp <- paste(input_root, f, sep="/")
  ouf <- paste0(strsplit(f, "_")[[1]][1], ".csv")
  outp <- paste(output_root, ouf, sep="/")
  print(inp)
  print(outp)
  transform_single(inp, outp)
  
}

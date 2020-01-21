# pheatmap for top highly expressed genes
install.packages('dplyr')
install.packages("pheatmap")
install.packages("tidyverse")


library('tidyverse')
library(ggplot2)
library(pheatmap)
library(dplyr)

# set folder to tpm files, and run the script using "ENSEMBL_GeneID_Conversion_Count_for_hg38_TPM_V2", 
# generating seurat_pbmc_data.csv. 

setwd("Z:/imPACT Technology/Sequencing data/EXP19-1748-51cells/tpm_files")

cd <- read_csv('seurat_pbmc_data-V2.csv')


#cd <- cd[!duplicated(cd$x),]
head(cd)

ncells <- dim(cd)[[2]]
View(colnames(cd[,2:ncells]))

wells <- colnames(cd[,2:ncells])

# filter rows if sum of  96 well expression value >10
cd <- cd %>%
  filter( rowSums(cd[,2:ncells]) > 10)

cd <- cd %>%
  filter(rowSums(cd[,2:ncells] != 0) > 5)
length(rowSums(cd[,2:ncells]!=0))

# either log2 or scale method can be used,and I am using scale method. 
mat <- log2(1+as.matrix(cd[,2:ncells]))
mat <- scale(as.matrix(cd[, 2:ncells]))

cells<-colnames(cd)
cells<-cells[2:ncells]
samples = "samples"
for(i in cells){
  x<-strsplit(i, '_')[[1]][1]
  samples <-c(samples,x)
  print(x)
}

# you can copy samples(c(samples,x)) in word document, use find/replace function of doc to make samples names,
# then assign it to colnames of mat or cd_sorted, then you can genrate a vector for types, colors, colnames, which
# can be used for tsne analysis; 
# for(j in samples[2:length(samples)]){
#   y<-strsplit(j, '8.')[[1]][1]
#   print(y)
# }


# types <-c("single_cell.0.5%lysis.Davis-ERCC",       
#           "single_cell.0.1% lysis.Quake",        
#           "single_cell.0.1% lysis.Quake",         
#           "single_cell.0.1% lysis.Davis",        
#           "single_cell.0.5%lysis.Quake",        
#           "single_cell.0.1% lysis.Quake",         
#           "single_cell.0.5%lysis.Quake-ERCC",       
#           "single_cell.0.5%lysis.Quake",        
#           "single_cell.0.1% lysis.Davis",        
#           "single_cell.0.5%lysis.Quake-ERCC",       
#           "single_cell.0.5%lysis.Davis-ERCC",       
#           "single_cell.0.5%lysis.Davis-ERCC",       
#           "single_cell.0.5%lysis.Davis",        
#           "single_cell.0.1% lysis.Quake-ERCC",        
#           "single_cell.0.5%lysis.Quake-ERCC",       
#           "single_cell.0.1% lysis.Davis-ERCC",        
#           "single_cell.0.5%lysis.Quake",        
#           "single_cell.0.5%lysis.Quake",        
#           "single_cell.0.1% lysis.Davis",      
#           "single_cell.0.5%lysis.Davis",        
#           "single_cell.0.5%lysis.Davis",     
#           "single_cell.0.1% lysis.Quake-ERCC",        
#           "single_cell.0.5%lysis.Quake-ERCC",       
#           "single_cell.0.1% lysis.Quake-ERCC",        
#           "single_cell.0.5%lysis.Quake-ERCC",       
#           "single_cell.0.5%lysis.Davis-ERCC",       
#           "single_cell.0.1% lysis.Davis-ERCC",        
#           "single_cell.0.5%lysis.Davis-ERCC",       
#           "single_cell.0.1% lysis.Quake",        
#           "single_cell.0.5%lysis.Quake",        
#           "single_cell.0.1% lysis.Quake",         
#           "single_cell.0.1% lysis.Davis",         
#           "single_cell.0.5%lysis.Davis",        
#           "single_cell.0.1% lysis.Davis",        
#           "single_cell.0.5%lysis.Davis",        
#           "0.1% lysis.20cells.Quake-ERCC",
#           "0.5%lysis.20cells.Quake-ERCC",
#           "0.5%lysis.20cells.Davis-ERCC",
#           "0.1% lysis.20cells.Davis-ERCC",
#           "0.1% lysis.20cells.Quake", 
#           "0.5%lysis.20cells.Quake", 
#           "0.1% lysis.20cells.Davis", 
#           "0.5%lysis.20cells.Davis")


# colors <-c("yellow",       
#           "blue",        
#            "blue",         
#            "green",        
#            "purple",        
#            "blue",         
#            "purple",       
#            "purple",        
#            "green",        
#            "purple",       
#            "yellow",       
#            "yellow",       
#            "yellow",        
#            "blue",        
#            "purple",       
#            "green",        
#            "purple",        
#            "purple",        
#            "green",      
#            "yellow",        
#            "yellow",     
#            "blue",        
#            "purple",       
#            "blue",        
#            "purple",       
#            "yellow",       
#            "green",        
#            "yellow",       
#            "blue",        
#            "purple",        
#            "blue",         
#            "green",         
#            "yellow",        
#            "green",        
#            "yellow",        
#            "red",
#            "orange",
#            "grey",
#            "black",
#            "red", 
#            "orange", 
#            "black", 
#            "grey")
           



# colnames(mat)<-c("A4.0.5%lysis.Davis-ERCC",       
#                 "A5.0.1% lysis.Quake",        
#                  "A6.0.1% lysis.Quake",         
#                  "A7.0.1% lysis.Davis",        
#                  "B5.0.5%lysis.Quake",        
#                  "B6.0.1% lysis.Quake",         
#                  "C1.0.5%lysis.Quake-ERCC",       
#                  "C6.0.5%lysis.Quake",        
#                  "C8.0.1% lysis.Davis",        
#                  "D2.0.5%lysis.Quake-ERCC",       
#                  "D3.0.5%lysis.Davis-ERCC",       
#                  "D4.0.5%lysis.Davis-ERCC",       
#                  "D8.0.5%lysis.Davis",        
#                  "E1.0.1% lysis.Quake-ERCC",        
#                  "E2.0.5%lysis.Quake-ERCC",       
#                  "E4.0.1% lysis.Davis-ERCC",        
#                  "E5.0.5%lysis.Quake",        
#                  "E6.0.5%lysis.Quake",        
#                  "E7.0.1% lysis.Davis",      
#                  "E7.0.5%lysis.Davis",        
#                  "E8.0.5%lysis.Davis",     
#                  "F1.0.1% lysis.Quake-ERCC",        
#                  "F1.0.5%lysis.Quake-ERCC",       
#                  "F2.0.1% lysis.Quake-ERCC",        
#                  "F2.0.5%lysis.Quake-ERCC",       
#                  "F3.0.5%lysis.Davis-ERCC",       
#                  "F4.0.1% lysis.Davis-ERCC",        
#                  "F4.0.5%lysis.Davis-ERCC",       
#                  ".F5.0.1% lysis.Quake",        
#                  "F5.0.5%lysis.Quake",        
#                  "F6.0.1% lysis.Quake",         
#                  "F7.0.1% lysis.Davis",         
#                  "F7.0.5%lysis.Davis",        
#                  "F8.0.1% lysis.Davis",        
#                  "F8.0.5%lysis.Davis",        
#                  "G1.0.1% lysis.20cells.Quake-ERCC",
#                  "G2.0.5%lysis.20cells.Quake-ERCC",
#                  "G3.0.5%lysis.20cells.Davis-ERCC",
#                  "G4.0.1% lysis.20cells.Davis-ERCC",
#                  "G6.0.1% lysis.20cells.Quake", 
#                  "G6.0.5%lysis.20cells.Quake", 
#                  "G7.0.1% lysis.20cells.Davis", 
#                  "G7.0.5%lysis.20cells.Davis"
#                  ) 

# sorting by sum of counts, in descending order. 
cd_sorted <- cd %>%arrange(desc(rowSums(cd[,2:ncells])))

mat_sorted <- log2(1+as.matrix(cd_sorted[,2:ncells]))
mat_sorted <- scale(as.matrix(cd_sorted[,2:ncells]))

row.names(mat_sorted)<- cd_sorted$X1
colnames(mat_sorted)<-c("A4.0.5%lysis.Davis-ERCC",       
                 "A5.0.1% lysis.Quake",        
                 "A6.0.1% lysis.Quake",         
                 "A7.0.1% lysis.Davis",        
                 "B5.0.5%lysis.Quake",        
                 "B6.0.1% lysis.Quake",         
                 "C1.0.5%lysis.Quake-ERCC",       
                 "C6.0.5%lysis.Quake",        
                 "C8.0.1% lysis.Davis",        
                 "D2.0.5%lysis.Quake-ERCC",       
                 "D3.0.5%lysis.Davis-ERCC",       
                 "D4.0.5%lysis.Davis-ERCC",       
                 "D8.0.5%lysis.Davis",        
                 "E1.0.1% lysis.Quake-ERCC",        
                 "E2.0.5%lysis.Quake-ERCC",       
                 "E4.0.1% lysis.Davis-ERCC",        
                 "E5.0.5%lysis.Quake",        
                 "E6.0.5%lysis.Quake",        
                 "E7.0.1% lysis.Davis",      
                 "E7.0.5%lysis.Davis",        
                 "E8.0.5%lysis.Davis",     
                 "F1.0.1% lysis.Quake-ERCC",        
                 "F1.0.5%lysis.Quake-ERCC",       
                 "F2.0.1% lysis.Quake-ERCC",        
                 "F2.0.5%lysis.Quake-ERCC",       
                 "F3.0.5%lysis.Davis-ERCC",       
                 "F4.0.1% lysis.Davis-ERCC",        
                 "F4.0.5%lysis.Davis-ERCC",       
                 ".F5.0.1% lysis.Quake",        
                 "F5.0.5%lysis.Quake",        
                 "F6.0.1% lysis.Quake",         
                 "F7.0.1% lysis.Davis",         
                 "F7.0.5%lysis.Davis",        
                 "F8.0.1% lysis.Davis",        
                 "F8.0.5%lysis.Davis",        
                 "G1.0.1% lysis.20cells.Quake-ERCC",
                 "G2.0.5%lysis.20cells.Quake-ERCC",
                 "G3.0.5%lysis.20cells.Davis-ERCC",
                 "G4.0.1% lysis.20cells.Davis-ERCC",
                 "G6.0.1% lysis.20cells.Quake", 
                 "G6.0.5%lysis.20cells.Quake", 
                 "G7.0.1% lysis.20cells.Davis", 
                 "G7.0.5%lysis.20cells.Davis"
) 

row.names(mat_sorted)<- cd_sorted$gene_name
colnames(mat_sorted)<- colnames(mat)     


# remove mitochondrial genes 
mat.sorted.data.MTindex = grep(pattern = "^MT-", x = rownames(mat_sorted), value = FALSE);

mat.sorted.woMT.data = mat_sorted[-mat.sorted.data.MTindex,];

write.csv(mat.sorted.woMT.data, "mat_woMT.csv")
# take only top 50 hightly expressed genes. 
#mat2 <-mat.sorted.woMT.data[1:50,]

# take off the genes beginning with "RPS and RPL", this can be Ribosome proteins
RPSindex = grep(pattern = "^RPS*", x = rownames(mat.sorted.woMT.data), value = FALSE);
# RPSindex
mat.woRPSindex = mat.sorted.woMT.data[-RPSindex,];
mat3<-mat.woRPSindex[1:50,]
#RPSindex = RPLindex
#RPLindex = grep(pattern = "^RPL*", x = rownames(mat.sorted.woMT.data), value = FALSE); 



# the default clustering method is average, I am using median here.

pheatmap(mat3, treeheigt_row = 0, treeheight_col = 0, clustering_method = "median", fontsize = 8)
#pheatmap(mat2, treeheigt_row = 0, treeheight_col = 0, clustering_method = "average", fontsize = 8)
# pheatmap can also be used on raw counts by cd_sorted
# pheatmap(cd_sorted[,2:ncells], scale = "row", treeheigt_row = 0, treeheight_col = 0, fontsize = 5)


# before scale: instead using scale function on mat, you can also use pheatmap to scale the data.
# edited by Linlin for pheatmap

#
# kmeans_k = 4, cell width= ""
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
# ph=pheatmap(avgreltlog10, color = colorRampPalette(rev(brewer.pal(n = 7, 
# name = "RdYlBu")))(100), 
# kmeans_k = NA, breaks = NA, border_color = "grey60",
# cellwidth = 10, cellheight=10, scale = "none", cluster_rows=TRUE,
# clustering_method = "average", cutree_rows = 4, cutree_cols= 2,)

# you can download Rcolorbrewer to bring more colors and 
# use colorRamPallete(rev(brewer.pal(n=7, name = "paired")))(100)) funcion

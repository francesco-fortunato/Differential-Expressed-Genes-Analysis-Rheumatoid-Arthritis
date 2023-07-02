rm(list = ls())

library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)

setwd("~/Scrivania/Bioinformatics/Exam/PCA")
path_in <- "~/Scrivania/Bioinformatics/Exam/Results/"
################################################

dirRes = "Results/"
if(!dir.exists(dirRes)){
  dir.create(dirRes)
}else { 
  print(paste("the directory ", dirRes,"already exist"))
}

dataset = "RA" # a folder for my case of study
dirDataset = paste0(dirRes,dataset,"/")
if(!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else { 
  print(paste("the directory ", dirDataset,"already exist"))
}

#file score plot is the plot of the score of pc1 against pc2 that are the old variable in the new 
# reference syste,
file_score_plot <- paste0(dirDataset, "score_plot.pdf")

#pareto scree plot is the plot about the expplained variance of each component
# how many component could be extends to retain?
file_pareto_scree_plot <- paste0(dirDataset, "pareto_scree_plot.pdf")
#loading plot is a mess

#contribution of each variable in each component
#for each component which are the gene that are contributing much? we select them
file_contribution_plot <- paste0(dirDataset, "PC_contribution_plot.pdf")

################################################
# 1. Importing data

data <- read.table(paste0(path_in, dataset, "/matrix_DEG.txt"),
                   header = T, sep = "\t", quote = "",
                   check.names = F, row.names = 1)

list_normal <- read.table (paste0(path_in, dataset, "/normal.txt"),
                           header = F, quote = "", check.names = F)$V1

list_case <- read.table (paste0(path_in, dataset, "/case.txt"),
                           header = F, quote = "", check.names = F)$V1

data <- t(data[,c(list_case, list_normal)])

groups <- c(rep("case", length(list_case)), rep("normal", length(list_normal)))
################################################
# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variables (geni)

pca <- prcomp(data, center = T, scale. = T, retx = T) 
################################################
# 3. Compute score and score plot
# (scores = the coordinates of old data (observations) in the new systems, that are the PCs)

# pca$x = t(data)*pca$rotation 
scores <- pca$x 

# alternative for score computation
# scores <- get_pca_ind(pca)$coord
# colnames(scores) <- paste0('PC', seq(1,ncol(scores)))

pdf(file_score_plot,width=5, height=5) 
g <- ggbiplot(pca, 
              obs.scale = 1, 
              var.axes = F, 
              ellipse = T, 
              groups = groups)
print(g)
dev.off() 

# alternative score plot
fviz_pca_ind(pca,
             col.ind = groups, # "cos2"
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE,
             repel = TRUE, # Avoid text overlapping
)

################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order (from the largest to the smallest)
eigenvalue = pca$sdev^2 

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0('PC', seq(1,length(varS)))

pdf(file_pareto_scree_plot, width =5, height=5)

# pareto chart
pareto.chart(varS[1:10])

# scree plot
fviz_eig(pca,addlabels = TRUE)
mean_lambda = 0.7*mean(eigenvalue)
dev.off()

# alternative for computation of eigenvalues and explained and cumulative variance
#df <- get_eigenvalue(pca)
#rownames(df) <- paste0('PC', seq(1,nrow(df)))
################################################
# 5. Compute loadings (coefficient of each PC)
# the matrix of variable loadings (a matrix whose columns contain the eigenvectors)
# matrice di trasformazione dalle vecchie alle nuove coordinate

#loadings <- pca$rotation

# loadings plot
# high cos2 = good representation of the variable on the principal component. 
# (i.e., the variable is positioned close to the circumference of the correlation circle)
# low cos2 = the variable is not perfectly represented by the PCs
# (i.e., the variable is close to the center of the circumference of the correlation circle)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

#fviz_pca_var(pca, col.var = "cos2", #"contrib" # Color by the quality of representation or by contribution of each avriable to the PC1-2
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#             repel = TRUE # Avoid text overlapping
#)
################################################
# Contributions of variables
# The function fviz_contrib() creates a barplot of row/column contributions. 
# A reference dashed line corresponds to the expected value if the contribution were uniform. 
# For a given dimension, any row/column with a contribution above the reference line could be considered 
# as important in contributing to the dimension.

# Contributions of variables to the PCs
contrib_var <- get_pca_var(pca)$contrib 
colnames(contrib_var) <- paste0('PC', seq(1,ncol(contrib_var)))

contrib_var <- contrib_var[order(contrib_var[,"PC1"], decreasing = T),]

plot.new()  # Create a new empty window
dev.off()   # Close the empty window

pdf(file_contribution_plot, width = 10, height = 10)

corrplot(contrib_var[1:10,1:10], is.corr = FALSE,
         tl.col = "black",
         method = "color",
         col = brewer.pal(n = 9, name = "BuPu"),
         addCoef.col = "black",
         tl.cex = 0.8)  # Adjust the text size if needed

# to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
# to PC3
fviz_contrib(pca, choice = "var", axes = 3, top = 10)

# total (PC1 + PC2 + PC3 + PC4)
fviz_contrib(pca, choice = "var", axes = 1:3, top = 15)

dev.off()




#Load 1 to 1 orthologous gene matrix and mask the first three columns.
OrthoMatrix<-read.csv(file = "ReducedOrthoMatrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rOrthoMatrix<-OrthoMatrix[,c(-1,-2,-3)]
strains<-colnames(rOrthoMatrix)

#Generate 50 random integers from 1 to 4718
set.seed(1)
genes<-runif(50,min = 1, max = 4717)

#For each strain, retireve the 50 orthologous genes and write them out in a file names "strain_50genes.txt"
for (i in 1:ncol(rOrthoMatrix)) {
  writeLines(rOrthoMatrix[genes,i],paste0(strains[i],"_50genes.txt"))  
}



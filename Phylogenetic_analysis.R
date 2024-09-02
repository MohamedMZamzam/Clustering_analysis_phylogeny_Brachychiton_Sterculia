#Needed libraries###
library (seqinr)
library (ape)
library (ggtree)
library(ggplot2)
library(Biostrings)
library(DECIPHER)
#Load the fasta format sequences file###
M <- readDNAStringSet ("D://Al-Azhar University/Dr Ahmed Al Tahir/rcbl/rbcl nonredundant accessions trimmed seq - Copy.fasta")

#align sequences###
aligned <- AlignSeqs(M)
writeXStringSet(aligned, "alignedd.fasta")
#read alignment###
DNA <- read.alignment("alignedd.fasta", format = "fasta")
Distance_M <- dist.alignment(DNA, matrix = "similarity")
temp <- as.data.frame(as.matrix(Distance_M))
#construct the tree###
tre <- nj (Distance_M)
tre <- ladderize(tre)
ggtree(tre, layout = "rectangular") + 
  geom_tiplab(size = 4)
write.tree(tre, "tree")


#test reliability by bootsraping####
boot.tree <- function(data, B = 100, tree = "upgma") {
  library(phangorn)
  func <- function(x) upgma(dist(x, method = "euclidean"), method="average")
  if (tree == "nj") {
    func <- function(x) nj(dist(x, method = "euclidean"))
  }
  if (tree == "hclust") {
    func <- function(x) {
      tr = hclust(dist(x, method = "euclidean"))
      tr = as.phylo(tr)
      return(tr)
    }
  }
  tr_real = func(data)
  plot(tr_real)
  library(ape)
  bp <- boot.phylo(tr_real, data, FUN=func, B=B)
  nodelabels(bp, frame = "none", )
  return(bp)
}

data = t(temp) # columns are the branches
boot.tree(data, B=1000, tree = "nj")



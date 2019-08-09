# Function to create List of distance/dissimilarity matrices 
#Methods in decostand "total", "max", "frequency", "normalize", "range", "rank", "rrank", "standardize", "pa",  
#Methods in vegdist "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".

DistList <- function(x, methods=c("totalbray", "totaleuclidean", "totalmanhattan"), relativize="before", namelist="distlist") {
  library(vegan)
  DistList <- list()
  if (sum(duplicated(methods))>0) {
    stop(print("Duplicated entries in method"))
  } 
  
  if (sum(methods=="totalbray")==1) {
    if (relativize=="before") {
      distmatrix<-vegdist(decostand(t(x), method="total"), 
                        method="bray")
      #Create list containing the dist matrices
      distName <- paste( 'totalbray' )
      DistList[[ distName ]]<-distmatrix
    }
  }
  
  if (sum(methods=="totaleuclidean")==1) {
    if (relativize=="before") {
      distmatrix<-vegdist(decostand(t(x), method="total"), 
                          method="euclidean")
      #Create list containing the dist matrices
      distName <- paste( 'totalbray' )
      DistList[[ distName ]]<-distmatrix
    }    
  }
  
  if (sum(methods=="totalmanhattan")==1) {
    if (relativize=="before") {
      distmatrix<-vegdist(decostand(t(x), method="total"), 
                          method="manhattan")
      #Create list containing the dist matrices
      distName <- paste( 'totalbray' )
      DistList[[ distName ]]<-distmatrix
    }    
  }
}

# Function to append lists together





DistList(method=c("euclidean", "manhattan", "euclidean"))
DistList(method=c("euclidean", "manhattan"))



library(microbiome)
library(dplyr)
data(dietswap)

#Create Metadata only include samples from HE study and at the first time point
Metadata <- data.frame(dietswap@sam_data@.Data)
colnames(Metadata) <- dietswap@sam_data@names
Metadata<-filter(Metadata, group == "HE" & timepoint.within.group == 1)
#table(Metadata$nationality)
Metadata$sample <- gsub("-", ".", Metadata$sample) 

#Create testdata 
Testdata <- data.frame(dietswap@otu_table@.Data) #dim(Testdata) #222 samples and 130 features
#Subset according to selected samples
Testdata<-select(Testdata, one_of(Metadata$sample))

#Remove features if they are not present in the subset of samples 
Testdata <- Testdata[rowSums(Testdata)>0,] #37 samples and 115 features

rm(dietswap)

DistList(Testdata)

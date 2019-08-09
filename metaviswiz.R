# Function to 

DistList <- function(x, method) {
  library(vegan)
  if (sum(duplicated(method))>0) {
    stop(print("Duplicated entries in method"))
  } 
  
  if (sum(method=="bray")==1) {
    print("bray")
  }
  
  if (sum(method=="euclidean")==1) {
    print("euclidean")
  }
  
  if (sum(method=="manhattan")==1) {
    print("manhattan")
  }
}

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


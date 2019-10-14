#############################################################################################
#Create testdata Start#
#############################################################################################

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

#############################################################################################
#Create testdata End#
#############################################################################################

# Function to create List of distance/dissimilarity matrices 
#Methods in decostand "total", "max", "frequency", "normalize", "range", "rank", "rrank", "standardize", "pa",  
#Methods in vegdist "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
#The data should be provided with samples as columns and features as rows 

CreateDistList <- function(x, methods=c("totalbray", "totaleuclidean", "totalmanhattan")) {
  library(vegan) #Have to remove
  DistList <- list() #Create empty list
  #Assess multiple of the same entries
  if (sum(duplicated(methods))>0) { #multiple equal entries evaluates as positive integers 
    methods<-unique(methods) #Remove duplicated entries
    warning("Duplicated entries in method removed")
  } 
  
  #Assess if nonexistant methods are specified
  #Remember to update when adding more methods
  possiblemethods<-c("totalbray", "totaleuclidean", "totalmanhattan", "totalhellingerbray", "totalhellingereuclidean", "totalhellingermanhattan") 
  #Remember to update when adding more methods
  if (length(setdiff(methods, possiblemethods))>0) { 
    for (i in 1:length(setdiff(methods, possiblemethods))) {
      warning(paste("Nonexistant methods specified", setdiff(methods, possiblemethods)[i]))
    }
  }   
  
  if (sum(methods=="totalbray")==1) {
    distmatrix<-vegdist(decostand(t(x), method="total"), 
                        method="bray") 
    #Create list containing the dist matrices
    DistList[[ 'totalbray' ]]<-distmatrix
  }
  
  if (sum(methods=="totaleuclidean")==1) {
    distmatrix<-vegdist(decostand(t(x), method="total"), 
                        method="euclidean")
    #Create list containing the dist matrices
    DistList[[ 'totaleuclidean' ]]<-distmatrix
  }
  
  if (sum(methods=="totalmanhattan")==1) {
    distmatrix<-vegdist(decostand(t(x), method="total"), 
                        method="manhattan")
    #Create list containing the dist matrices
    DistList[[ 'totalmanhattan' ]]<-distmatrix
  }
  
  if (sum(methods=="totalhellingerbray")==1) {
    distmatrix<-vegdist(decostand(decostand(t(x), method="total"), 
                                  method="hellinger"), method="bray")
    #Create list containing the dist matrices
    DistList[[ 'totalhellingerbray' ]]<-distmatrix
  }
  
  if (sum(methods=="totalhellingereuclidean")==1) {
    distmatrix<-vegdist(decostand(decostand(t(x), method="total"), 
                                  method="hellinger"), method="euclidean")
    #Create list containing the dist matrices
    DistList[[ 'totalhellingereuclidean' ]]<-distmatrix
  }
  
  if (sum(methods=="totalhellingermanhattan")==1) {
    distmatrix<-vegdist(decostand(decostand(t(x), method="total"), 
                                  method="hellinger"), method="manhattan")
    #Create list containing the dist matrices
    DistList[[ 'totalhellingermanhattan' ]]<-distmatrix
  }
  
  return(DistList)
}


#############################################################################################
#Testing CreateDistList Start#
#############################################################################################

DistList <- CreateDistList(x=Testdata, method=c("dssd", "sdadg", "totalbray", "totalmanhattan", "totalbray"))

#############################################################################################
#Testing CreateDistList End#
#############################################################################################



# Function to append two lists together, to be used when user have created their own list and want to add to an already created DistList   
# Stop if provided is not lists
# Would like to implement more than two lists can be provided
AppendTwoLists <- function(x, y) {
  if ((class(x)!="list") | (class(y)!="list")) {
    stop("Provide as lists, only two lists can be provided at a time")
  }
  
  for (i in 1:length(x)) {
    if (class(x[[i]])!="dist") {
      stop("Provide the first list containing dist objects. If distance or dissimilarity object
      is provided as a matrix or dataframe convert using eg.
           as.dist(yourdistmatrix, diag = FALSE, upper = FALSE")
    }
  }
  
  for (i in 1:length(y)) {
    if (class(y[[i]])!="dist") {
      stop("Provide the second list containing dist objects. If distance or dissimilarity object
      is provided as a matrix or dataframe convert using eg.
           as.dist(yourdistmatrix, diag = FALSE, upper = FALSE")
    }
  }
  
  AppendedLists<-c(x, y)
  
  #Look for duplicate names 
  if (sum(duplicated(names(AppendedLists)))>0) { 
    stop("Duplicated names in appended List")
  } 
  
  for (i in 1:length(AppendedLists)) {
    if ( sum(labels(AppendedLists[[1]]) != labels(AppendedLists[[i]]))>0 ) {
      print(i)
      warning("Labels in distobjects are not the same.
           Incongruence in comparison of labels in the first dist object
           with the dist obejct corresponding to the number printed above")
    }
  }
  
  return(AppendedLists)
}

#############################################################################################
#Testing AppendTwoLists Start#
#############################################################################################

#Not updated accordingly 
#Realised that if c() ... is provided it already concatenates before creates problems when providing eg. distmatrix
#Have problems with assess elements in the dist objects. can use labels()
#User have to use the dist format otherwise modify the provided matrix

List1<-CreateDistList(Testdata, method=c("totalbray"))
List2<-CreateDistList(Testdata, method=c("totaleuclidean", "totalmanhattan"))

distmatrix<-as.matrix(vegdist(decostand(t(Testdata), method="total"), 
                    method="manhattan"))
distmatrixList<-list(distmatrix)
class(distmatrixList[[1]])

#Does it work
List3<-AppendTwoLists(List1, List2)
#Provided one list as a matrix getting error
List4<-AppendTwoLists(List2, distmatrix)
#Provided second list that contained something that was not a distobject
List4<-AppendTwoLists(List2, distmatrixList)
#Provided first list that contained something that was not a distobject
List4<-AppendTwoLists(distmatrixList, List2)

distmatrixdims<-as.dist(distmatrix, diag = FALSE, upper = TRUE)
#distmatrixdims
distmatrixdims2<-vegdist(decostand(t(Testdata), method="total"), 
        method="manhattan")
#distmatrixdims2
#distmatrixdims==distmatrixdims2 #It is the same don't think it matters then on the following analysis

#Changing labels 
vectorlabels<-labels(List3[[1]])
vectorlabels2<-c(vectorlabels[2:length(vectorlabels)], vectorlabels[1])
vectorlabels3<-c(vectorlabels[2:length(vectorlabels)], "Sample.1231322123")
#Basic one sample is missing
Testdata2<-Testdata[2:length(Testdata), 2:length(Testdata)]
newlabels<-list(vegdist(decostand(t(Testdata2), method="total"), 
                   method="manhattan"))
List3<-AppendTwoLists(List2, newlabels)
#Renaming Testdata columns
Testdata2<-Testdata
colnames(Testdata2)<-vectorlabels2
newlabels<-list(vegdist(decostand(t(Testdata2), method="total"), 
                        method="manhattan"))
List3<-AppendTwoLists(List2, newlabels)

##At some point I would like to make a function that can do it for more than one list
#How do I specify in a function that there is a variable number of entries (x, y, z...)
#AppendMultipleLists <- function(x) {
#  len<-length(x)
#  print(len)
#  for (i in 1:len) {
#    if (class(x[i])!="list") {
#      stop("Provide vector of list names c(List1, List2...")
#      print(i)
#    }
#  }
#}
#List3<-AppendMultipleLists(c(List1, List2))

#############################################################################################
#Testing AppendTwoLists End#
#############################################################################################

#PairProtest
#Function wrapper around protest that creates correlations as distances 1-correlation, raw correlations and sum of squares to a list
#Input List of dist objects 
#The output can then be used as input to MetaVisCor and MetaVisWiz

#Different checkpoints inherited from the AppendMultipleLists function
#If changing do it both places
PairProtest <- function(x) {
  library(vegan) #Have to remove
  
  if ((class(x)!="list")) {
    stop("Provide a list of dist objects")
  }
  
  for (i in 1:length(x)) {
    if (class(x[[i]])!="dist") {
      stop("Provide list containing dist objects. If distance or dissimilarity object
      is provided as a matrix or dataframe convert using eg.
           as.dist(yourdistmatrix, diag = FALSE, upper = FALSE")
    }
  }
  
  #Look for duplicate names 
  if (sum(duplicated(names(x)))>0) { 
    stop("Duplicated names in appended List")
  } 
  
  for (i in 1:length(x)) {
    if ( sum(labels(x[[1]]) != labels(x[[i]]))>0 ) {
      print(i)
      warning("Labels in distobjects are not the same.
           Incongruence in comparison of labels in the first dist object
           with the dist obejct corresponding to the number printed above")
    }
  }
  
  #The actual running of Procrustes analysis
  #Create empty df with col and row names according to names in the dist list
  ProCruPairCordist<-setNames(data.frame(matrix(ncol=length(x), nrow=length(x))), c(names(x)))  
  row.names(ProCruPairCordist)<-c(names(ProCruPairCordist))
  ProCruPairCor <- ProCruPairCordist
  ProCruPairSS <- ProCruPairCordist
  #Create empty list to hold results from protest
  ListProCru <- list()
  for (i in 1:length(x)) { 
    for (j in 1:length(x)) {  
      #Make pairwise protest. Comment nested for loop j is iterated first then i
      prot<-protest(capscale(x[[i]]~1), capscale(x[[j]]~1))
      ProCruPairCordist[j,i]<-(1-prot$t0)
      ProCruPairCor[j,i]<-(prot$t0)
      ProCruPairSS[j,i]<-(prot$ss)
      
      #if (method == "cordist") {
      #  ProCruPair[j,i]<-(1-prot$t0)
      #} else if (method == "cor") {
      #  ProCruPair[j,i]<-(prot$t0)
      #} else if (method == "ss") {
      #  ProCruPair[j,i]<-(1-prot$ss)
      #} else {
      #  print("Specify valid method")
      #} 
    }
    #Percentage finished indicator
    print(paste("Percentage run", (i/length(x))*100))
  }
  
  ListProCru[[ "Cordist" ]]<-as.dist(ProCruPairCordist)
  ListProCru[[ "Cor" ]]<-as.dist(ProCruPairCor)
  ListProCru[[ "SS" ]]<-as.dist(ProCruPairSS)
  #The as.dist function default is diag=FALSE, upper=FALSE, auto_convert_data_frames=TRUE.
  
  return(ListProCru)
}

#############################################################################################
#Testing PairProtest Start#
#############################################################################################

DistList <- CreateDistList(x=Testdata, method=c("totaleuclidean", "totalmanhattan", "totalbray", "totalhellingermanhattan"))
##Does the function method work
#testing1<-PairProtest(DistList, method="cordist")
#testing2<-PairProtest(DistList, method="cor")
#testing3<-PairProtest(DistList, method="ss")
testing<-PairProtest(DistList)

#Is the dist matrix binded correctly 
prot<-protest(capscale(DistList$totalbray~1), capscale(DistList$totalhellingermanhattan~1))
prot$t0
prot2<-protest(capscale(DistList$totaleuclidean~1), capscale(DistList$totalmanhattan~1))
prot2$t0

#############################################################################################
#Testing PairProtest End#
#############################################################################################


#MetaVisWiz
#Function creating the meta PCoA plot from the dist object of sum of squares or 1 - the procrustes correlation
#Input list created from PairProtest that contains 1-cor, cor and SS
#x<-Pairpro
MetaVisWiz <- function(x, method="Cordist") {
  library(vegan) #Have to remove
  library(reshape2) #Have to remove includes melt
  library(tidyr) #Have to remove includes unite
  library(gridExtra) #Have to remove include grid.arrange
  library(corrplot)
  
  #Comparisons number limited by input 
  
  #Create empty list to hold plots
  FigureList <- list()
  
  #Create capscale object
  if (method == "Cordist") {
    PCoAObject<-capscale(x$Cordist~1)
  } else if (method == "SS") {
      PCoAObject<-capscale(x$SS~1)
      } else {stop("Specify valid method")}  
  
  ###############################################
  #Make stressplot
  #Extract ordination distances and merge with observed dissimilarity
  #Correlations
  stress<-stressplot(PCoAObject)
  df <- melt(as.matrix(stress)) #Long format ordination distance
  names(df)<-c("rowOrd", "colOrd", "OrdDist")
  #df<-filter(df, OrdDist>0) #Remove comparisons to same method does include the same comparison twice
  
  if (method == "Cordist") {
    df2 <- melt(as.matrix(x$Cordist))
  } else if (method == "SS") {
      df2 <- melt(as.matrix(x$SS))
      }  else {stop("Specify valid method")}  
  names(df2)<-c("rowObs", "colObs", "ObsDism")
  #df2<-filter(df2, ObsDism>0) #Remove comparisons to same method does include the same comparison twice
  #tidyr::unite: Convenience function to paste together multiple columns into one.
  df<-unite(df, mergecol, c(rowOrd, colOrd), remove=FALSE) 
  df2<-unite(df2, mergecol, c(rowObs, colObs), remove=FALSE)
  ggstress<-merge(df, df2, by="mergecol")
  #Plot stress plot
  FigureList$Stress<-ggplot(ggstress) + 
    geom_point(aes(ObsDism, OrdDist), size=1) + 
    #Size depend on number of comparisons, can implement size differences accordingly in if statement
    ggtitle("Stress plot") + 
    labs(x = "Observed dissimilarity", y = "Ordination distance") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12))
  
  ###############################################
  #Create PCoA
  ##Add eig to plot axes. with cmdscale there are negative values not with capscale
  eig <- PCoAObject$CA$eig
  # Calculate the variation explained by PCoA1, 2, 3 and 4
  # and use it to generate axis labels
  eig_1_2 <- eig[1:4] / sum(eig) * 100
  eig_1 <- paste("PCoA1", round(eig_1_2[1], digits = 2), "% variance")
  eig_2 <- paste("PCoA2", round(eig_1_2[2], digits = 2), "% variance")
  #Additional eigen values if user want to modify code to plot additional MDS
  #eig_3 <- paste("PCoA3", round(eig_1_2[3], digits = 2), "% variance")
  #eig_4 <- paste("PCoA4", round(eig_1_2[4], digits = 2), "% variance")
  ##Pull out coordinates for plotting from the ca object
  #Structuring to add to Metadata2
  PCoACA<-PCoAObject$CA 
  #The ca object contains the actual ordination results: 
  #u ((Weighted) orthonormal site scores), 
  #v ((Weighted) orthonormal species scores) all na in mine, 
  #Xbar (The standardized data matrix after previous stages of analysis), 
  #and imaginary.u.eig ???. Info http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/cca.object.html   
  PCoA<-as.data.frame(PCoACA$u)
  #Change colnames. Now add dis and trans info to names 
  colnames(PCoA) <- c("MDS1","MDS2")
  #If user want to modify code to plot additional MDS
  #colnames(PCoA) <- c("MDS1","MDS2", "MDS3","MDS4")
  #Add row names to df
  PCoA$Sample <- row.names(PCoA)
  
  #Create metadata 
  MetadataProC<-data.frame(Sample=rownames(PCoA))
  MetadataProC$Trans <-
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*totallog", MetadataProC$Sample), "TSS_log",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*totalhellinger", MetadataProC$Sample), "TSS_hellinger",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*clr", MetadataProC$Sample), "clr",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*ilr", MetadataProC$Sample), "ilr",        
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*total", MetadataProC$Sample), "TSS",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*max", MetadataProC$Sample), "max",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*freq", MetadataProC$Sample), "freq",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*normalize", MetadataProC$Sample), "norm",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*pa", MetadataProC$Sample), "pa",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*hellinger", MetadataProC$Sample), "hellinger",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*log", MetadataProC$Sample), "log",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*chi.square", MetadataProC$Sample), "chisq",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*CSS", MetadataProC$Sample), "CSS",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*TMM", MetadataProC$Sample), "TMM",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*DESeq", MetadataProC$Sample), "DESeq",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*Rar", MetadataProC$Sample), "Rarefy",
           "Other"))))))))))))))))
  
  MetadataProC$Dist <-
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*manhattan*", MetadataProC$Sample), "manhattan",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*euclidean*", MetadataProC$Sample), "euclidean",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*canberra*", MetadataProC$Sample), "canberra",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*bray*", MetadataProC$Sample), "bray",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*kulczynski*", MetadataProC$Sample), "kulczynski",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*jaccard*", MetadataProC$Sample), "jaccard",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*gower*", MetadataProC$Sample), "gower",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*altGower*", MetadataProC$Sample), "altGower",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*horn*", MetadataProC$Sample), "horn",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*binomial*", MetadataProC$Sample), "binomial",
    ifelse(seq(along=(MetadataProC$Sample)) %in% grep("*clr", MetadataProC$Sample), "euclidean",
           "Other")))))))))))
  
  #Merge according to Sample
  MetadataProC2<-merge(MetadataProC, PCoA, by="Sample")
  
  #Define coloring schemes to create procrustes PCoA/PCA
  ProCCol<-c("chisq" = "#E41A1C", 
             "freq" = "#4DAF4A", 
             "max" = "#FFFF33", 
             "norm" = "#A65628", 
             "clr"="#E41A1C", 
             "ilr"="#9f1214", 
             "CSS"="#ed5e5f", 
             "DESeq"="#377EB8", 
             "TMM"="#265880", 
             "hellinger"="#4DAF4A", 
             "TSS_hellinger"="#357933", 
             "log"="#984EA3", 
             "TSS_log"="#5b2e61", 
             "pa"="#FF7F00", 
             "Rarefy"="#FFFF33", 
             "TSS"="#A65628")
  #Define shape scheme to create procrustes PCoA/PCA
  ProCShape<-c("altGower"=9, 
               "binomial"=5, 
               "bray"=3, 
               "canberra"=6, 
               "euclidean"=0, 
               "gower"=8, 
               "horn"=2, 
               "jaccard"=1, 
               "kulczynski"=7, 
               "manhattan"=4) 
  
  #Plot PCoA
  FigureList$PCoA<-ggplot(MetadataProC2) + 
    #geom_line(aes(x=MDS1, y=MDS2, group=Matching_samples), size=0.1, linetype="dotted") +  
    geom_jitter(aes(MDS1, MDS2, col=Trans, shape=Dist), width=0.00, height=0.00, alpha=0.8, size=3, stroke=1.5) +
    #geom_point(aes(MDS1, MDS2, color = Sample_LPSX, group = Sample, shape = Temperature), size=5) +
    scale_color_manual(values=ProCCol) +  
    scale_shape_manual(values=ProCShape) +
    ggtitle("PCoA") + 
    labs(colour="Preprocessing", shape="beta-diversity", x = eig_1, y = eig_2) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12), legend.position="bottom") #+
  #stat_ellipse(data=filter(MetadataProC2[,1:6], grepl("clr|ilr", Trans)), aes(MDS1, MDS2), level=0.80)
  #scale_y_reverse() #If you want the y scale reversed, to make plots easier to compare
  #scale_x_reverse() #If you want the x scale reversed, to make plots easier to compare
  #ggsave(paste("CapscalePCoAYoungCOMPARE", "Genus", i, OrgFlt, ".pdf", sep=""), height=6, width=12)
  #ggtitle(paste("PCoA ", i, " n orgs > 1% = ", nrow(Tax2))) #Used with filtering
  
  ###############################################
  #Create Scree plot
  screeplot<-data.frame(PCoAObject$CA$eig)
  colnames(screeplot)<-c("eig")
  screeplot$eig <- screeplot$eig[1:length(screeplot$eig)] / sum(screeplot$eig) * 100
  screeplot<-tibble::rownames_to_column(screeplot, "MDS")
  screeplot$MDS <- factor(screeplot$MDS, levels=c(sprintf("MDS%d", 1:length(screeplot$eig))))

  #Plot screeplot
  FigureList$Scree<-ggplot(screeplot, aes(x=MDS, y=eig)) + 
    geom_bar(stat="identity") + 
    labs(x ="MDS", y ="eig (%)") + 
    ggtitle(paste("Scree plot ")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.title=element_text(size=12), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  #ggsave(filename=paste("ScreeplotCapscalePCoA", "Genus", i, OrgFlt, ".pdf", sep=""), height=6, width=12)
  
  ###############################################
  #Create histogram
  #Create histogram / density plot for correlations. Sensitivity analysis 
  #The histogram and density plot is to be used as an assessment of further investigations are needed in order to decide on how to process metagenomics data. 
  #Correls<-melt_dist(matrix(x$Cor)) #Previous solution, but can now remove harrietr
  Correls<-data.frame(dist=as.vector(x$Cor))
  FigureList$DensityCor<-ggplot(Correls, aes(x=dist)) + 
    geom_histogram(aes(y=..density..), binwidth=0.01, colour="black", fill="white") + #..density.. is to have it scale with geom_density
    geom_density(alpha=.25, fill="#FF6666") + 
    labs(x ="Correlation", y ="Density") + 
    ggtitle(paste("Density plot ")) +
    xlim(0,1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=12)) 
  FigureList$DensityCo
  ###############################################
  #Create output figure with all plots
  #Have the plots stored in FigureList
  lay <- rbind(c(1,1,1,1),
               c(1,1,1,1),
               c(2,2,3,4))
  ## Create figure with correlations
  pdf(paste("Metaviswiz", ".pdf", sep=""), width=12.5, height=10)
  grid.arrange(FigureList$PCoA, FigureList$DensityCor, FigureList$Stress, FigureList$Scree, layout_matrix = lay)
  dev.off()
  
  ###############################################
  #Create correlelograms
  #User can modify directly with corrplot 
  pdf(paste("MetaviswizCorgram", ".pdf", sep=""), width=12.5, height=10)
  corrplot(as.matrix(x$Cor), 
           type = "upper", 
           order = "hclust", 
           tl.col = "black", 
           tl.srt = 25,
           addCoef.col = "white",
           diag=FALSE, 
           cl.lim = c(0, 1))
  dev.off()
  
  return(FigureList)
}

#############################################################################################
#Testing MetaVisWiz Start#
#############################################################################################

DistList <- CreateDistList(x=Testdata, method=c("totaleuclidean", "totalmanhattan", "totalbray", "totalhellingermanhattan", "totalhellingerbray",  "totalhellingereuclidean"))
Pairpro<-PairProtest(DistList)
VisWiz<-MetaVisWiz(Pairpro)



#############################################################################################
#Testing MetaVisWiz End#
#############################################################################################


#AllOrdi
#Function to create individual PCA and PCoA
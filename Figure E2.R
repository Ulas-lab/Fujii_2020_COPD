###All .fcs files can be found http://flowrepository.org/ as ID FR-FCM-Z2JL and exported to .csv files with FlowJo software.

##Figure E2A
#Load packages
BiocManager::install("flowCore")
require(flowCore)
install.packages("umap")
require(umap)
install.packages("ggplot2")
require(ggplot2)
install.packages("ggsci")
require(ggsci)
install.packages("wesanderson")
require(wesanderson)
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("JinmiaoChenLab/Rphenograph")
require(Rphenograph)
devtools::install_github('VPetukhov/ggrastr')
require(ggrastr)

#Transform function
nfTransform <- function(dataA, dataB){
  dataA1<-dataA
  dataB1<-dataB
  for(i in as.character(keeptable[,1])) {
    q<-0.05
    m<-4.5
    d <- dataA[,i]
    w <- 0
    t <- max(d)
    nd <- d[d < 0]
    nThres <- quantile(nd, 0.25) - 1.5 * IQR(nd)
    nd <- nd[nd >= nThres]
    if (length(nd)) {
      r <- .Machine$double.eps + quantile(nd, q)
      if (10^m * abs(r) <= t) {
        w <- 0
      }
      else {
        w <- (m - log10(t/abs(r)))/2
        if (is.nan(w) || w > 2) {
          warning(paste0("autoLgcl failed for channel: ",
                         p, "; using default fluor logicle transformation!"))
          w <- 0.1
          t <- 500000
          m <- 4.5
        }
      }
    }
    templgcl <- logicleTransform(w=w, t=t, m=4.5, a=0)
    dataNum <- which(colnames(dataA)==i)
    temp <- apply(dataA[,dataNum,drop=F],2, templgcl)
    dataA1[,dataNum] <- temp
    dataB1[,dataNum] <- temp
    print(paste0(i, " w= ",w," t= ",t))
  }
  return(list(dataA1=dataA1, dataB1=dataB1))
}


##Downsample
ceil = 10000

##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs/Pt141/"
LoaderPATH <- paste(PATH, "fcs", sep="/")
FcsFileNames <- list.files(path = LoaderPATH, pattern = ".fcs")
fs = list()
for(FileNum in 1:length(FcsFileNames)){
  fs[[FileNum]] <- read.FCS(paste0(LoaderPATH,"/",FcsFileNames[FileNum]),transformation =FALSE,ignore.text.offset=T,truncate_max_range=FALSE,emptyValue = F)
}
filenames <- FcsFileNames
NumBC <- length(fs)
FFdata <- NULL
OrigNames <- fs[[1]]@parameters$name
FFt <- exprs(fs[[1]])
if (ceil==0){FFa <- FFt} else {
  down.sample <- sample(c(1:nrow(FFt)), ceil)      
  FFa <- FFt[down.sample,]
}
colnames(FFa) <- fs[[1]]@parameters$desc
empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ")
colnames(FFa)[empties] <- fs[[1]]@parameters$name[empties]
fs[[1]]@parameters$desc <- colnames(FFa)
FFa <- cbind(FFa,rep(1,dim(FFa)[1]))
colnames(FFa)[dim(FFa)[2]] <- "InFile"
FFdata <- rbind(FFdata,FFa)
prepData <- list(FFdata=FFdata, OrigNames=OrigNames, forNewFF=fs[[1]], NumBC=NumBC, FcsFileNames = FcsFileNames)
forNewFF <- prepData$forNewFF
keeptable <- read.csv(paste(PATH, "Pt141.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt141 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs/Pt141/export_Myeloid_BAL_015_CD45 positive cells.csv", header = TRUE)
data.facs.pt141 <- data.facs.pt141[down.sample,]

##Normalize Data
nfTransOut <- nfTransform(data, FFdata)
data1 <- nfTransOut$dataA1
FFdata1 <- nfTransOut$dataB1

##UMAP
umapMat <- NULL
umap_conf <- umap.defaults
umap_conf$random_state <- 123
umap_conf$n_neighbors <- 15
umap_conf$n_components <- 2
umap_conf$min_dist <- 0.2
umap_conf$metric <- "euclidean"

umapMat <- umap(d = data1, config = umap_conf)
umap.pt141 <- umapMat$layout
colnames(umap.pt141) <- c("UMAP1","UMAP2")
plot(umap.pt141[, 1], umap.pt141[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=3)

#Featureplot visualization
library(RColorBrewer)
UMAP_FACS <- cbind(umap.pt141, data.facs.pt141)
color = colorRampPalette(rev(brewer.pal(n = 13, name = "RdYlBu")))(10)
for (i in (colnames(UMAP_FACS))) {
  gg <- ggplot(UMAP_FACS, aes(x=UMAP1 , y=UMAP2,color=get(i)))+
    geom_point_rast(size=0.5,alpha=1)+
    scale_color_gradientn(colours=color)+labs(title=i)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
  
  print(gg)
}

##Figure E2C
#Transform function
nfTransform <- function(dataA, dataB){
  dataA1<-dataA
  dataB1<-dataB
  for(i in as.character(keeptable[,1])) {
    q<-0.05
    m<-4.5
    d <- dataA[,i]
    w <- 0
    t <- max(d)
    nd <- d[d < 0]
    nThres <- quantile(nd, 0.25) - 1.5 * IQR(nd)
    nd <- nd[nd >= nThres]
    if (length(nd)) {
      r <- .Machine$double.eps + quantile(nd, q)
      if (10^m * abs(r) <= t) {
        w <- 0
      }
      else {
        w <- (m - log10(t/abs(r)))/2
        if (is.nan(w) || w > 2) {
          warning(paste0("autoLgcl failed for channel: ",
                         p, "; using default fluor logicle transformation!"))
          w <- 0.1
          t <- 500000
          m <- 4.5
        }
      }
    }
    templgcl <- logicleTransform(w=w, t=t, m=4.5, a=0)
    dataNum <- which(colnames(dataA)==i)
    temp <- apply(dataA[,dataNum,drop=F],2, templgcl)
    dataA1[,dataNum] <- temp
    dataB1[,dataNum] <- temp
    print(paste0(i, " w= ",w," t= ",t))
  }
  return(list(dataA1=dataA1, dataB1=dataB1))
}

##Concatenate the FCS-files
PATH <- "/Users/fujiiwataru/Desktop/My article/Pt141/fcs/"
FcsFileNames <- list.files(path = PATH, pattern = ".fcs")
fs = list()
for(FileNum in 1:length(FcsFileNames)){
  fs[[FileNum]] <- read.FCS(paste0(PATH,"/",FcsFileNames[FileNum]),transformation =FALSE,ignore.text.offset=T,truncate_max_range=FALSE,emptyValue = F)
}
filenames <- FcsFileNames
FFdata <- NULL
for(i in 1:length(filenames)){
  FFa <- exprs(fs[[i]])
  colnames(FFa) <- fs[[i]]@parameters$desc
  empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ")
  colnames(FFa)[empties] <- fs[[i]]@parameters$name[empties]
  fs[[i]]@parameters$desc <- colnames(FFa)
  FFa <- cbind(FFa,rep(1,dim(FFa)[1]))
  colnames(FFa)[dim(FFa)[2]] <- "InFile"
  FFa <- as.data.frame(FFa)
  FFa$gate <- filenames[i]
  FFdata <- as.data.frame(rbind(FFdata,FFa))
}
row.names(FFdata) <- paste(seq(1,nrow(FFdata),1), FFdata$gate, sep="_")
FFdata$gate <- NULL 

##Define gates to keep
keeptable <- read.csv("/Users/fujiiwataru/Desktop/My article/Pt141/BAL_myeloid.csv")
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]

##Normalize Data
nfTransOut <- nfTransform(data, FFdata)
data1 <- nfTransOut$dataA1
FFdata1 <- nfTransOut$dataB1

##UMAP
umapMat <- NULL
umap_conf <- umap.defaults
umap_conf$random_state <- 123
umap_conf$n_neighbors <- 15
umap_conf$n_components <- 2
umap_conf$min_dist <- 0.2
umap_conf$metric <- "euclidean"
umapMat <- umap(d = data1, config = umap_conf)
umap <- umapMat$layout
colnames(umap) <- c("UMAP1","UMAP2")
plot(umap[, 1], umap[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=0.1)


##Plot self-made UMAP
umap <- as.data.frame(umap)
umap$gate <- row.names(umap)
umap$gate <- gsub("\\..*","", gsub(".*\\_", "", umap$gate))
cells.col <- c("alveolar-macrophage" = "#43A132",
               "eosinophils" = "#3079B2",
               "dendritic-cell" = "#AACFE2",
               "mast-cells" = "#F5989A",
               "monocytes" = "#DB0024",
               "neutrophil" = "#F8BD74",
               "lymphocytes" = "#B5DF8D")
gg <- ggplot(umap, aes(UMAP1, UMAP2, color=gate)) + geom_point_rast(size = .3) +
  scale_colour_manual(values = cells.col) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5))
plot(gg)

##Remove legends
gg <- ggplot(umap, aes(UMAP1, UMAP2, color=gate)) + geom_point_rast(size = .3) +
  scale_colour_manual(values = cells.col) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.position = "none")
plot(gg)

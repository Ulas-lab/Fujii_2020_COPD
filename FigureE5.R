###All .fcs files can be found http://flowrepository.org/ as ID FR-FCM-Z2JL and exported to .csv files with FlowJo software.

##Figure E5A
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
ceil = 5000

##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs2/Pt162/"
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
keeptable <- read.csv(paste(PATH, "Pt162.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt162 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs2/Pt162/export_Lymphoid cells_BAL162_015_Lin negative cells.csv", header = TRUE)
data.facs.pt162 <- data.facs.pt162[down.sample,]

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
umap.pt162 <- umapMat$layout
colnames(umap.pt162) <- c("UMAP1","UMAP2")
plot(umap.pt162[, 1], umap.pt162[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=3)

##Phenograph clustering
Rphenograph_out <- Rphenograph(data1, k = 28)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
if(dim(Rphenographmat)[1] != dim(umap.pt162)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt162)[1]-dim(Rphenographmat)[1]), newclustnum))
  losts <- setdiff(1:dim(umap.pt162)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt162 <- as.data.frame(umap.pt162)
if(grepl("Phenograph", names(umap.pt162)[length(names(umap.pt162))])){umap.pt162$Phenograph <- NULL}
umap.pt162$Phenograph <- NULL
umap.pt162 <- cbind(umap.pt162, Rphenographmat)
umap.pt162 <- as.data.frame(umap.pt162)
umap.pt162$Phenograph <- as.character(umap.pt162$Phenograph)
gg <- ggplot(umap.pt162, aes(UMAP1, UMAP2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt162$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
plot(gg)

##Figure E5B
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
ceil = 5000

##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs2/Pt162/"
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
keeptable <- read.csv(paste(PATH, "Pt162.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt162 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs2/Pt162/export_Lymphoid cells_BAL162_015_Lin negative cells.csv", header = TRUE)
data.facs.pt162 <- data.facs.pt162[down.sample,]

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
umap.pt162 <- umapMat$layout
colnames(umap.pt162) <- c("UMAP1","UMAP2")
plot(umap.pt162[, 1], umap.pt162[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=3)

##Phenograph clustering
Rphenograph_out <- Rphenograph(data1, k = 28)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))

#Deal with cells that sometimes go missing with Phenograph
if(dim(Rphenographmat)[1] != dim(umap.pt162)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt162)[1]-dim(Rphenographmat)[1]), " events lost during running of Rphenograph - dunno why. Assigned as new cluster #", newclustnum))
  losts <- setdiff(1:dim(umap.pt162)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt162 <- as.data.frame(umap.pt162)
if(grepl("Phenograph", names(umap.pt162)[length(names(umap.pt162))])){umap.pt162$Phenograph <- NULL}
umap.pt162$Phenograph <- NULL
umap.pt162 <- cbind(umap.pt162, Rphenographmat)
umap.pt162 <- as.data.frame(umap.pt162)
umap.pt162$Phenograph <- as.character(umap.pt162$Phenograph)
gg <- ggplot(umap.pt162, aes(UMAP1, UMAP2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt162$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
plot(gg)

##Combined clusters
combined.clust <- umap.pt162
# combined.clust[combined.clust$Phenograph == 1,"ID"] <- "Airway epithelial cells"
# combined.clust[combined.clust$Phenograph %in% c(2,16,12),"ID"] <- "Lymphocytes"
combined.clust[combined.clust$Phenograph %in% c(3),"ID"] <- "B cells"
combined.clust[combined.clust$Phenograph %in% c(13),"ID"] <- "Double negative T cells"
combined.clust[combined.clust$Phenograph %in% c(4,5,6,7,14),"ID"] <- "CD4 T cells"
combined.clust[combined.clust$Phenograph %in% c(1),"ID"] <- "NK cells"
combined.clust[combined.clust$Phenograph %in% c(9,11,12),"ID"] <- "CD8 T cells"
combined.clust[combined.clust$Phenograph %in% c(10),"ID"] <- "non-defined"
combined.clust[combined.clust$Phenograph %in% c(2),"ID"] <- "ILCs"
combined.clust[combined.clust$Phenograph %in% c(8),"ID"] <- "Double positive T cells"
pal <- brewer.pal(name="Paired", n=8)
gg <- ggplot(combined.clust, aes(UMAP1, UMAP2, color = ID)) + 
  geom_point_rast(size=0.5)+
  # scale_color_manual(values = c(pals::stepped(n=7)))+
  scale_colour_manual(values = pal)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Clusters")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) 

plot(gg)

##Coloring by clusters
clusters <- as.character(unique(umap.pt186$Phenograph))
for(i in 1:length(clusters)){
  tmp.df <- umap.pt186
  tmp.df$color <- "tmp" 
  tmp.df[tmp.df$Phenograph != clusters[i], "color"] <- "white"
  tmp.df[tmp.df$Phenograph == clusters[i], "color"] <- "red"
  tmp.df <- tmp.df[order(tmp.df$color, decreasing = T),]
  print(ggplot(tmp.df, aes(x=UMAP1, y=UMAP2)) + geom_point(colour = ifelse(tmp.df$color == "red", "red", "gray80"), size=0.01) + ggtitle(clusters[i]) +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.x = element_line(color="black", size = 1),
                axis.title.x=element_blank(),
                axis.line.y = element_line(color="black", size = 1)))
}


##Figure E5C
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

###Template sample
##Downsample
ceil = 5000

##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs2/Pt162/"
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
keeptable <- read.csv(paste(PATH, "Pt162.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt162 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs2/Pt162/export_Lymphoid cells_BAL162_015_Lin negative cells.csv", header = TRUE)
data.facs.pt162 <- data.facs.pt162[down.sample,]

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
umap.pt162 <- umapMat$layout
colnames(umap.pt162) <- c("UMAP1","UMAP2")
plot(umap.pt162[, 1], umap.pt162[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=3)

##Phenograph clustering
Rphenograph_out <- Rphenograph(data1, k = 28)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
if(dim(Rphenographmat)[1] != dim(umap.pt162)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt162)[1]-dim(Rphenographmat)[1]), newclustnum))
  losts <- setdiff(1:dim(umap.pt162)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt162 <- as.data.frame(umap.pt162)
if(grepl("Phenograph", names(umap.pt162)[length(names(umap.pt162))])){umap.pt162$Phenograph <- NULL}
umap.pt162$Phenograph <- NULL
umap.pt162 <- cbind(umap.pt162, Rphenographmat)
umap.pt162 <- as.data.frame(umap.pt162)
umap.pt162$Phenograph <- as.character(umap.pt162$Phenograph)
gg <- ggplot(umap.pt162, aes(UMAP1, UMAP2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt162$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
plot(gg)

###Next sample
##Downsample
ceil = 5000

#data1 is template dataset and data2 the projected
lisa.knnClust <- knn(train = data1, test = data2, k = 1, cl = Rphenographmat)

#Visualization
umap.pt177 <- as.data.frame(cbind(umap.pt177, lisa.knnClust))
umap.pt177$umap.pt177.clusters <- as.character(umap.pt177$lisa.knnClust)
names(umap.pt177)[4] <- "phenograph"
ggplot(umap.pt177, aes(V1, V2, color = phenograph)) +
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt177$phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - predicted phenograph")


##Phenograph clustering
Rphenograph_out <- Rphenograph(data2, k = 40)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
if(dim(Rphenographmat)[1] != dim(umap.pt177)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt177)[1]-dim(Rphenographmat)[1]), newclustnum))
  losts <- setdiff(1:dim(umap.pt177)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt177 <- as.data.frame(umap.pt177)
if(grepl("Phenograph", names(umap.pt177)[length(names(umap.pt177))])){umap.pt177$Phenograph <- NULL}
umap.pt177$Phenograph <- NULL
umap.pt177 <- cbind(umap.pt177, Rphenographmat)
umap.pt177 <- as.data.frame(umap.pt177)
umap.pt177$Phenograph <- as.character(umap.pt177$Phenograph)
ggplot(umap.pt177, aes(V1, V2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt177$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")

##Repeat until last sample
#Combine all datasets
umap.pt162$ID <- "162"
umap.pt162$disease <- "COPD"
names(umap.pt162)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt177$ID <- "177"
umap.pt177$disease <- "Control"
names(umap.pt177)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt172$ID <- "172"
umap.pt172$disease <- "Control"
names(umap.pt172)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt167$ID <- "167"
umap.pt167$disease <- "Control"
names(umap.pt167)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt123$ID <- "123"
umap.pt123$disease <- "Control"
names(umap.pt123)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt116$ID <- "116"
umap.pt116$disease <- "Control"
names(umap.pt116)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt111$ID <- "111"
umap.pt111$disease <- "Control"
names(umap.pt111)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt77$ID <- "77"
umap.pt77$disease <- "Control"
names(umap.pt77)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt51$ID <- "51"
umap.pt51$disease <- "Control"
names(umap.pt51)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt115$ID <- "115"
umap.pt115$disease <- "COPD"
names(umap.pt115)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt103$ID <- "103"
umap.pt103$disease <- "COPD"
names(umap.pt103)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt94$ID <- "94"
umap.pt94$disease <- "COPD"
names(umap.pt94)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt83$ID <- "83"
umap.pt83$disease <- "COPD"
names(umap.pt83)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt80$ID <- "80"
umap.pt80$disease <- "COPD"
names(umap.pt80)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt54$ID <- "54"
umap.pt54$disease <- "COPD"
names(umap.pt54)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt48$ID <- "48"
umap.pt48$disease <- "COPD"
names(umap.pt48)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt44$ID <- "44"
umap.pt44$disease <- "COPD"
names(umap.pt44)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.combined <- as.data.frame(rbind(umap.pt162, umap.pt177, umap.pt172, umap.pt167, umap.pt123, umap.pt116, umap.pt111, umap.pt77, umap.pt51, umap.pt115, umap.pt103, umap.pt94, umap.pt83, umap.pt80, umap.pt54, umap.pt48, umap.pt44))

# Create data frame with cell counts
test<-umap.combined[!is.na(umap.combined$phenograph),]
df<-test %>% group_by(phenograph, ID, disease) %>% summarise(count = n())
df

# Figure 2C
gg <- ggplot(df, aes(x=phenograph, y=count, fill=disease)) + 
  geom_bar(stat='identity',colour='black', position = position_dodge())+
  xlab("Cluster")+
  ylab("Number of cells")+
  ggtitle("")+
  scale_fill_manual(values=c('blue2','#E31A1C'),labels=c('Control','COPD'))+
  theme(text=element_text(size=14),
        plot.title=element_text(hjust=0.5,face='bold'),
        axis.title.x=element_text(size=14,face='bold'),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=14,colour='black'),
        axis.text.y=element_text(size=14,colour='black'),
        legend.title=element_text(size=13,face='bold',hjust=0),
        legend.background = element_rect(fill='white',size=2),
        legend.key=element_rect(colour='black',fill='white'),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_rect(fill=NA))

##Figure E5D
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
ceil = 5000

##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs2/Pt162/"
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
keeptable <- read.csv(paste(PATH, "Pt162.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt162 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs2/Pt162/export_Lymphoid cells_BAL162_015_Lin negative cells.csv", header = TRUE)
data.facs.pt162 <- data.facs.pt162[down.sample,]

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
umap.pt162 <- umapMat$layout
colnames(umap.pt162) <- c("UMAP1","UMAP2")
plot(umap.pt162[, 1], umap.pt162[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=3)

##Phenograph clustering
Rphenograph_out <- Rphenograph(data1, k = 28)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))

#Deal with cells that sometimes go missing with Phenograph
if(dim(Rphenographmat)[1] != dim(umap.pt162)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt162)[1]-dim(Rphenographmat)[1]), " events lost during running of Rphenograph - dunno why. Assigned as new cluster #", newclustnum))
  losts <- setdiff(1:dim(umap.pt162)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt162 <- as.data.frame(umap.pt162)
if(grepl("Phenograph", names(umap.pt162)[length(names(umap.pt162))])){umap.pt162$Phenograph <- NULL}
umap.pt162$Phenograph <- NULL
umap.pt162 <- cbind(umap.pt162, Rphenographmat)
umap.pt162 <- as.data.frame(umap.pt162)
umap.pt162$Phenograph <- as.character(umap.pt162$Phenograph)
gg <- ggplot(umap.pt162, aes(UMAP1, UMAP2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt162$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
plot(gg)

#Featureplot visualization
library(RColorBrewer)
UMAP_FACS <- cbind(umap.pt162, data.facs.pt162)
color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(10)
for (i in (colnames(UMAP_FACS))) {
  gg <- ggplot(UMAP_FACS, aes(x=UMAP1 , y=UMAP2,color=get(i)))+
    geom_point_rast(size=0.5,alpha=1)+
    scale_color_gradientn(colours=color)+labs(title=i)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=.5)) 
  
  print(gg)
}

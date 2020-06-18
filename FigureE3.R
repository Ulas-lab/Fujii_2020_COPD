###All .fcs files can be found http://flowrepository.org/ as ID FR-FCM-Z2JL and exported to .csv files with FlowJo software.
##Figure E3A
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
library(flexclust)
library(class)
library(reshape2)
library(ggplot2)
library(dplyr)

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

##Tamplate sample
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

##Phenograph clustering
Rphenograph_out <- Rphenograph(data1, k = 40)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
if(dim(Rphenographmat)[1] != dim(umap.pt141)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt141)[1]-dim(Rphenographmat)[1]), newclustnum))
  losts <- setdiff(1:dim(umap.pt141)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt141 <- as.data.frame(umap.pt141)
if(grepl("Phenograph", names(umap.pt141)[length(names(umap.pt141))])){umap.pt141$Phenograph <- NULL}
umap.pt141$Phenograph <- NULL
umap.pt141 <- cbind(umap.pt141, Rphenographmat)
umap.pt141 <- as.data.frame(umap.pt141)
umap.pt141$Phenograph <- as.character(umap.pt141$Phenograph)
ggplot(umap.pt141, aes(UMAP1, UMAP2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt141$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")

##Next sample
##Folder Preparation
PATH <- "/Users/fujiiwataru/Desktop/umap_facs/Pt130/"
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
keeptable <- read.csv(paste(PATH, "Pt130.csv", sep="/"))
data <- FFdata[,which (colnames(FFdata) %in% as.character(keeptable[,1]))]
data.facs.pt130 <- read.csv("/Users/fujiiwataru/Desktop/umap_facs/Pt130/export_Myeloid cells_BAL 130_015_CD45 positive cells.csv", header = TRUE)
data.facs.pt130 <- data.facs.pt130[down.sample,]

##Normalize Data
nfTransOut <- nfTransform(data, FFdata)
data2 <- nfTransOut$dataA1
FFdata2 <- nfTransOut$dataB1

#data1 is template dataset and data2 the projected
lisa.knnClust <- knn(train = data1, test = data2, k = 1, cl = Rphenographmat)

#Visualization
umap.pt130 <- as.data.frame(cbind(umap.pt130, lisa.knnClust))
umap.pt130$umap.pt130.clusters <- as.character(umap.pt130$lisa.knnClust)
names(umap.pt130)[4] <- "phenograph"
ggplot(umap.pt130, aes(V1, V2, color = phenograph)) +
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt130$phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - predicted phenograph")

##Phenograph clustering
Rphenograph_out <- Rphenograph(data2, k = 50)
Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
if(dim(Rphenographmat)[1] != dim(umap.pt130)[1]){
  rpnames <- as.numeric(row.names(Rpmat))
  newclustnum <- max(Rphenographmat)+1
  print(paste0((dim(umap.pt130)[1]-dim(Rphenographmat)[1]), newclustnum))
  losts <- setdiff(1:dim(umap.pt130)[1],rpnames)
  for(li in 1:length(losts)) {
    Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
}
umap.pt130 <- as.data.frame(umap.pt130)
if(grepl("Phenograph", names(umap.pt130)[length(names(umap.pt130))])){umap.pt130$Phenograph <- NULL}
umap.pt130 <- cbind(umap.pt130, Rphenographmat)
umap.pt130 <- as.data.frame(umap.pt130)
umap.pt130$Phenograph <- as.character(umap.pt130$Phenograph)
umap.pt130$Phenograph <- factor(umap.pt130$Phenograph, levels = c(as.character(seq(1,17, 1))))
ggplot(umap.pt130, aes(V1, V2, color = Phenograph)) + 
  geom_point(size=1)+
  scale_color_manual(values = c(pals::brewer.paired(n=length(unique(umap.pt130$Phenograph)))))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP - Phenograph")

###Repeat until last sample

##Combine all datasets
umap.pt141$ID <- "141"
umap.pt141$disease <- "Control"

umap.pt130$ID <- "130"
umap.pt130$disease <- "COPD"
names(umap.pt130)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt142$ID <- "142"
umap.pt142$disease <- "COPD"
names(umap.pt142)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt177$ID <- "177"
umap.pt177$disease <- "Control"
names(umap.pt177)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt172$ID <- "172"
umap.pt172$disease <- "Control"
names(umap.pt172)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt167$ID <- "167"
umap.pt167$disease <- "Control"
names(umap.pt167)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt146$ID <- "146"
umap.pt146$disease <- "Control"
names(umap.pt146)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt139$ID <- "139"
umap.pt139$disease <- "Control"
names(umap.pt139)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt128$ID <- "123"
umap.pt128$disease <- "COPD"
names(umap.pt128)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt176$ID <- "176"
umap.pt176$disease <- "COPD"
names(umap.pt176)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt173$ID <- "173"
umap.pt173$disease <- "COPD"
names(umap.pt173)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt162$ID <- "162"
umap.pt162$disease <- "COPD"
names(umap.pt162)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt161$ID <- "161"
umap.pt161$disease <- "COPD"
names(umap.pt161)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt160$ID <- "160"
umap.pt160$disease <- "COPD"
names(umap.pt160)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt157$ID <- "157"
umap.pt157$disease <- "COPD"
names(umap.pt157)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt155$ID <- "155"
umap.pt155$disease <- "COPD"
names(umap.pt155)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt152$ID <- "152"
umap.pt152$disease <- "COPD"
names(umap.pt152)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt145$ID <- "145"
umap.pt145$disease <- "COPD"
names(umap.pt145)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt175$ID <- "175"
umap.pt175$disease <- "COPD"
names(umap.pt175)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt164$ID <- "164"
umap.pt164$disease <- "Control"
names(umap.pt164)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt135$ID <- "135"
umap.pt135$disease <- "Control"
names(umap.pt135)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt133$ID <- "133"
umap.pt133$disease <- "Control"
names(umap.pt133)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.pt137$ID <- "137"
umap.pt137$disease <- "COPD"
names(umap.pt137)[c(1,2)] <- c("UMAP1", "UMAP2")

umap.combined <- as.data.frame(rbind(umap.pt141, umap.pt130, umap.pt142, umap.pt177, umap.pt172, umap.pt167, umap.pt146, umap.pt139, umap.pt128, umap.pt176, umap.pt173, umap.pt162, umap.pt161, umap.pt160, umap.pt157, umap.pt155, umap.pt152,umap.pt145, umap.pt175,umap.pt164,umap.pt135, umap.pt133, umap.pt137))

# Create data frame with cell counts
test<-umap.combined[!is.na(umap.combined$phenograph),]
df<-test %>% group_by(phenograph, ID, disease) %>% summarise(count = n())
df

# Add cluster and disease annotation to flow data
test1<-cbind(data1,umap.pt141[,c('phenograph','ID')])
test2<-cbind(data2,umap.pt130[,c('phenograph','ID')])
test3<-cbind(data3,umap.pt142[,c('phenograph','ID')])
test4<-cbind(data4,umap.pt177[,c('phenograph','ID')])
test5<-cbind(data5,umap.pt172[,c('phenograph','ID')])
test6<-cbind(data6,umap.pt167[,c('phenograph','ID')])
test7<-cbind(data7,umap.pt146[,c('phenograph','ID')])
test8<-cbind(data8,umap.pt139[,c('phenograph','ID')])
test9<-cbind(data9,umap.pt128[,c('phenograph','ID')])
test10<-cbind(data10,umap.pt176[,c('phenograph','ID')])
test11<-cbind(data11,umap.pt173[,c('phenograph','ID')])
test12<-cbind(data12,umap.pt162[,c('phenograph','ID')])
test13<-cbind(data13,umap.pt161[,c('phenograph','ID')])
test14<-cbind(data14,umap.pt160[,c('phenograph','ID')])
test15<-cbind(data15,umap.pt157[,c('phenograph','ID')])
test16<-cbind(data16,umap.pt155[,c('phenograph','ID')])
test17<-cbind(data17,umap.pt152[,c('phenograph','ID')])
test18<-cbind(data18,umap.pt145[,c('phenograph','ID')])
test19<-cbind(data19,umap.pt175[,c('phenograph','ID')])
test20<-cbind(data20,umap.pt164[,c('phenograph','ID')])
test21<-cbind(data21,umap.pt135[,c('phenograph','ID')])
test22<-cbind(data22,umap.pt133[,c('phenograph','ID')])
test23<-cbind(data23,umap.pt137[,c('phenograph','ID')])

# Merge all the patients together
all.patients<-rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12,test13,test14,test15,test16,test17,test18,test19,test20,test21,test22,test23)

# Subset each patient for all macrophage clusters
tmp<-list()
for (i in c(2,6,7,10,11)){
  tmp[[i]]<-subset(all.patients, phenograph==i)
}
cluster.2<-as.data.frame(tmp[[2]])
cluster.6<-as.data.frame(tmp[[6]])
cluster.7<-as.data.frame(tmp[[7]])
cluster.10<-as.data.frame(tmp[[10]])
cluster.11<-as.data.frame(tmp[[11]])

# Identify which patients are controls and which COPD in the dataset
controls<-unique(all.patients[all.patients$disease=='Control', ]$ID)
copds<-unique(all.patients[all.patients$disease=='COPD', ]$ID)

## Calculate the mean for control and COPD patients in cluster 2
# CD14 for control patients
vector.con<-vector()
for (i in controls){
  vector.con[i]<-mean(cluster.2[cluster.2$ID==i,]$CD14)
}
vector.con

# CD14 for COPD patients
vector.copd<-vector()
for (i in copds){
  vector.copd[i]<-mean(cluster.2[cluster.2$ID==i,]$CD14)
}
vector.copd
table.2<-data.frame(mean(vector.con),sd(vector.con)/sqrt(length(vector.con)),mean(vector.copd),sd(vector.copd)/sqrt(length(vector.copd)))

# HLA-DR for control patients
vector.con<-vector()
for (i in controls){
  vector.con[i]<-mean(cluster.2[cluster.2$ID==i,]$`HLA.DR`)
}
vector.con

# HLA-DR for COPD patients
vector.copd<-vector()
for (i in copds){
  vector.copd[i]<-mean(cluster.2[cluster.2$ID==i,]$`HLA.DR`)
}
vector.copd
aa<-data.frame(mean(vector.con),sd(vector.con)/sqrt(length(vector.con)),mean(vector.copd),sd(vector.copd)/sqrt(length(vector.copd)))
table.2<-cbind(table.2,aa)

# Autofluorescence for control patients
vector.con<-vector()
for (i in controls){
  vector.con[i]<-mean(cluster.2[cluster.2$ID==i,]$Autofluorescence)
}
vector.con

# Autofluorescence for COPD patients
vector.copd<-vector()
for (i in copds){
  vector.copd[i]<-mean(cluster.2[cluster.2$ID==i,]$Autofluorescence)
}
vector.copd
aa<-data.frame(mean(vector.con),sd(vector.con)/sqrt(length(vector.con)),mean(vector.copd),sd(vector.copd)/sqrt(length(vector.copd)))
table.2<-cbind(table.2,aa)

# FSC for control patients
vector.con<-vector()
for (i in controls){
  vector.con[i]<-mean(cluster.2[cluster.2$ID==i,]$`FSC.A`)
}
vector.con

# FSC for COPD patients
vector.copd<-vector()
for (i in copds){
  vector.copd[i]<-mean(cluster.2[cluster.2$ID==i,]$FSC.A)
}
vector.copd
aa<-data.frame(mean(vector.con),sd(vector.con)/sqrt(length(vector.con)),mean(vector.copd),sd(vector.copd)/sqrt(length(vector.copd)))
table.2<-cbind(table.2,aa)

# SSC for control patients
vector.con<-vector()
for (i in controls){
  vector.con[i]<-mean(cluster.2[cluster.2$ID==i,]$`SSC.A`)
}
vector.con

# SSC for COPD patients
vector.copd<-vector()
for (i in copds){
  vector.copd[i]<-mean(cluster.2[cluster.2$ID==i,]$SSC.A)
}
vector.copd
aa<-data.frame(mean(vector.con),sd(vector.con)/sqrt(length(vector.con)),mean(vector.copd),sd(vector.copd)/sqrt(length(vector.copd)))
table.2<-cbind(table.2,aa)

# Format the table for cluster 2
table.2<-melt(table.2)
table.2$marker<-c(rep('CD14',4),rep('HLA.DR',4),rep('Autofluorescence',4),rep('FSC.A',4),rep('SSC.A',4))
table.2$variable<-rep(c('mean','sem'),10)
table.2$disease<-rep(c('Control','COPD'),each=2)
table.2$cluster<-2

# Add the new information of cluster 2 to the final table
final.table<-rbind(final.table,table.2)
final.table

##Repeat until last cluster
# Prepare final table for plotting
final.table<-dcast(final.table, formula = marker+cluster+disease~variable, fun= mean)
final.table$cluster<-factor(final.table$cluster)

# Plot graphs
gg <- ggplot(final.table,aes(x=cluster, y=mean, fill=disease))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),colour='black',width=0.8,cex=0.8,position=position_dodge(0.9)) +
  geom_bar(stat='identity',position=position_dodge(0.9),colour='black')+
  scale_fill_manual(labels=c('Control','COPD'), values=c('#0000FF','#E3142C'))+
  facet_wrap(~marker,ncol=5)+
  ylab('Fluorescence units')+
  theme(text=element_text(size=14),
        axis.title.x=element_text(),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=14,face='bold'),
        axis.text.y=element_text(size=14,colour='black'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray84',colour='black'),
        legend.title=element_blank(),
        legend.position = 'top',
        legend.key=element_rect(colour='black',fill='white'),
        legend.spacing.x = unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='black',fill=NA))
plot(gg)

##Figure E3B
# Numbers for control patients
vector.con.Number<-vector()
for (i in controls){
  vector.con.Number[i]<-nrow(cluster.2[cluster.2$ID==i,])
}
vector.con.Number

# Numbers for COPD patients
vector.copd.Number<-vector()
for (i in copds){
  vector.copd.Number[i]<-nrow(cluster.2[cluster.2$ID==i,])
}
vector.copd.Number
table.2<-data.frame(mean(vector.con.Number),sd(vector.con.Number),mean(vector.copd.Number),sd(vector.copd.Number))

# Format the table for cluster 2
table.2<-melt(table.2)
table.2$marker<-c(rep('Number',4))
table.2$variable<-rep(c('mean','sd'),2)
table.2$disease<-rep(c('Control','COPD'),each=2)
table.2$cluster<-2

##Repeat until last cluster
# Add the new information of cluster 6 to the final table
final.table2<-rbind(table.2,table.6, table.7, table.10, table.11)
final.table2

# Prepare final table for plotting
final.table2<-dcast(final.table2, formula = marker+cluster+disease~variable, fun= mean)
final.table2$cluster<-factor(final.table2$cluster)


# Plot graphs
gg <- ggplot(final.table2,aes(x=cluster, y=mean, fill=disease))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),colour='black',width=0.8,cex=0.8,position=position_dodge(0.9)) +
  geom_bar(stat='identity',position=position_dodge(0.9),colour='black')+
  scale_fill_manual(labels=c('Control','COPD'), values=c('#0000FF','#E3142C'))+
  facet_wrap(~marker,ncol=5)+
  ylab('Number of cells')+
  theme(text=element_text(size=14),
        axis.title.x=element_text(),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=14,face='bold'),
        axis.text.y=element_text(size=14,colour='black'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray84',colour='black'),
        legend.title=element_blank(),
        legend.position = 'top',
        legend.key=element_rect(colour='black',fill='white'),
        legend.spacing.x = unit(0.600, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='black',fill=NA))
plot(gg)

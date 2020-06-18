makeObject <- function(data = lipidomics_data, 
                       data_annotation = sample_annotation, 
                       sample_column_start = 11, 
                       sample_column_end = NULL){
  # preparing species intensitiy data #################
  species.int <- data[data$IS == 0,]
  # remove all species which have NA values
  species.int <- species.int[!is.na(species.int$IS),]
  rownames(species.int) <- species.int$NAME
  
  # species row annotation
  species.int.rowanno <- species.int[,1:10]
  
  # species intensity data
  if(is.null(sample_column_end)){
    sample_column_end <- ncol(species.int)
  }
  species.int <- species.int[,sample_column_start:sample_column_end]
  
  # prepare sample annotation data ###########################
  
  
  if(!is.null(data_annotation)){
    # check if sample names from annotation file and the intensitiy file match
    if(identical(colnames(species.int), data_annotation$sample) == F){
      samples_to_keep <- intersect(colnames(species.int), data_annotation$sample)
      species.int <- species.int[,samples_to_keep]
      data_annotation <- data_annotation[data_annotation$sample %in% samples_to_keep,]
      print(paste(length(samples_to_keep), "samples are included in the analysis:", paste(head(samples_to_keep), collapse = ","), "...", paste(tail(samples_to_keep), collapse = ",")))
      }
    
    species.int.colanno <- data_annotation
    rownames(species.int.colanno) <- species.int.colanno$sample
  } else {
    species.int.colanno <- data.frame(sample = colnames(species.int), row.names = colnames(species.int))
  }
  
  # standards
  is.int <- data[data$IS == 1,]
  # remove rows that do not contain lipid species
  is.int <- is.int[!is.na(is.int$IS),]
  rownames(is.int) <- is.int$NAME
  # standards row annotation
  is.int.rowanno <- is.int[,1:10]
  
  # standards intensity data
  is.int <- is.int[,11:ncol(is.int)]
  
  # combine to summarized expreriment
  
  # create summarized experiment object
  lipids <- SummarizedExperiment(assays=list(intensity=as.matrix(species.int)),
                                 rowData=species.int.rowanno, 
                                 colData = species.int.colanno)
  standards <- SummarizedExperiment(assays=list(intensity=as.matrix(is.int)),
                                    rowData=is.int.rowanno)
  return(list(lipids = lipids, standards = standards))
}


getPmol <- function(lipid_object = lipids, 
                    standards_object = lip_obj$standards, 
                    filterNA = TRUE)
{
  
  # check if pmol slot is already present in the object
  if("pmol" %in% assayNames(lipid_object)){
    print("Pmol have already been calculated and will be overwritten")
  }
  
  # subset standards to contain the same samples as lipids
  standards_object <- standards_object[,colnames(lipid_object)]
  
  # get lipid classes
  classes <- unique(rowData(lipid_object)$CLASS)
  # check if all classes have a standards
  if (sum(!classes %in% rowData(lipid_object)$CLASS) != 0) {
    print(paste("Warning: no standards are provided for the following lipid classes: ", 
                classes[!classes %in% rowData(lipid_object)$CLASS], 
                ". The respective pmol values cannot be calculated and will be set to NA.", sep = " "))
  }
  result <- data.frame()
  for(i in classes){
    print(i)
    species <- rowData(lipid_object)[rowData(lipid_object)$CLASS == i, ]$NAME
    species <- assay(lipid_object)[species,]
    standard <- rowData(standards_object)[rowData(standards_object)$CLASS == i, ]$NAME
    standard <- as.numeric(t(assay(standards_object)[standard,]))
    if(sum(is.na(standard)) != 0){
      print(paste("Warning: no standard available for samples from lipid class", i))
    }
    # divide each sample by its standard
    if(is.vector(species)){
      pmol <- species/standard
    } else {
      pmol <- sweep(as.matrix(species), 2, standard, "/")
    }
    
    # multiply by pmol of standard
    ispmol <- rowData(lipid_object)[rowData(lipid_object)$CLASS == i, ]$ISpmol[1]
    pmol <- pmol * ispmol
    if(is.vector(pmol)){
      name <- rowData(lipid_object)[rowData(lipid_object)$CLASS == i, ]$NAME
      result <- rbind(result, pmol)
      rownames(result)[nrow(result)] <- name
    } else {
      result <- rbind(result, pmol)
    }
  }
  assay(lipid_object, 2) <- result
  assayNames(lipid_object)[2] <- "pmol"
  if(filterNA == T){
    keep <- complete.cases(assay(lipid_object,2))
    print(paste("Filtering out", sum(!keep), "species with no standard. Filtered species: ", paste(rownames(lipids[!keep,]), collapse = " ")))
    lipid_object <- lipid_object[keep,]
  }
  return(lipid_object)
}


getMolpermil <- function(lipid_object = lipids, 
                         stable_classes = c("PCHO","PHSE","SPMY", "PETA")){
  pmol <- as.data.frame(assay(lipid_object, "pmol"))
  pmol_sum <- colSums(pmol, na.rm = T)
  molpermil <- sweep(as.matrix(pmol), 2, pmol_sum/1000, "/")
  assay(lipid_object, 3) <- molpermil
  assayNames(lipid_object)[3] <- "molpermil"
  
  for(i in stable_classes){
    if(!i %in% rowData(lipid_object)$CLASS){
      warning(paste("lipid class ", i, " not found in data", sep = ""))
    }
  }
  selected_species <- rowData(lipid_object)[rowData(lipid_object)$CLASS %in% stable_classes,]$NAME
  pmol_sum <- colSums(pmol[rownames(pmol) %in% selected_species,], na.rm = T)
  molpermil <- sweep(as.matrix(pmol), 2, pmol_sum/1000, "/")
  assay(lipid_object, 4) <- molpermil
  assayNames(lipid_object)[4] <- "molpermil_stable"
  return(lipid_object)
}


getSaturation <- function(lipid_object = lipids){
  
  # check if pmol slot is already present in the object
  if("saturation" %in% assayNames(lipid_object)){
    print("Saturation already has been calculated and will be overwritten")
  }
  # get lipid classes
  classes <- unique(rowData(lipid_object)$CLASS)
  
  result <- data.frame()
  for(i in classes){
    print(paste("normalizing class", i))
    species <- rowData(lipid_object)[rowData(lipid_object)$CLASS == i, ]$NAME
    pmol <- assay(lipid_object, 3)[species,]
    
    double_bonds <- rowData(lipid_object)[species, ]$db
    
    # divide each sample by its double bonds
    saturation <- sweep(as.matrix(pmol), 1, double_bonds, "*")
    saturation <- sweep(as.matrix(saturation), 2, colSums(pmol), "/")
    result <- rbind(result, saturation)
  }
  
  assay(lipid_object, 4) <- saturation
  assayNames(lipid_object)[4] <- "saturation"
  return(lipid_object)
}

getSummary <- function(){
  
  summary <- as.data.frame(assay(lipids))
  summary$species <- rownames(summary)
  summary <- gather(summary, key = sample, value = intensity, -species)
  summary <- merge(summary, sample_annotation, by = "sample")
  summary <- merge(as.data.frame(rowData(lipids)), summary, by.x = "NAME", by.y = "species")
  
  # add identifier 
  summary$ID <- paste(summary$sample, summary$NAME, sep = "_")
  if("pmol" %in% assayNames(lipids)){
    # add pmol
    pmol <- assay(lipids, "pmol")
    pmol$species <- rownames(pmol)
    pmol <- gather(pmol, key = sample, value = pmol, -species)
    pmol$ID <- paste(pmol$sample, pmol$species, sep = "_")
    summary <- merge(summary, pmol[,c("ID", "pmol")], by = "ID")
  }
  
  # add molpermil
  if("molpermil" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil")], by = "ID")
  }
  
  # add molpermil_log
  if("molpermil_log" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil_log"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil_log, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil_log")], by = "ID")
  }
  
  if("molpermil_log_batchrm" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil_log_batchrm"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil_log_batchrm, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil_log_batchrm")], by = "ID")
  }
  
  # add molpermil_batchrm
  if("molpermil_batchrm" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil_batchrm"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil_batchrm, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil_batchrm")], by = "ID")
  }
  
  
  # add molpermil_batchrm_norm
  if("molpermil_batchrm_norm" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil_batchrm_norm"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil_batchrm_norm, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil_batchrm_norm")], by = "ID")
  }
  
  
  
  # add molpermil_stable
  if("molpermil_stable" %in% assayNames(lipids)){
    molpermil <- as.data.frame(assay(lipids, "molpermil_stable"))
    molpermil$species <- rownames(molpermil)
    molpermil <- gather(molpermil, key = sample, value = molpermil_stable, -species)
    molpermil$ID <- paste(molpermil$sample, molpermil$species, sep = "_")
    summary <- merge(summary, molpermil[,c("ID", "molpermil_stable")], by = "ID")
  }
  return(summary)
  
}


scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}
# Network visualization ----------------------------------------------------------------------------------

visualize_network <- function(network, 
                              color.by,
                              select.cluster = NULL,
                              plot.subnetwork = NULL,
                              gene.label = NULL,
                              class.label = NULL,
                              use.layout = "layout_with_fr",
                              save.pdf=T) {
  
  
  # prepare data
  network$TF_info <- factor(x = network$TF_info,
                            levels = c("TF","Chromatin_remodeller","Co_factor","RNBP","noTF"))
  
  network$cluster_color <- as.character(network$cluster_color)
  
  cluster.col <- network$cluster_color
  names(cluster.col) <- network$cluster_color
  
  
  
  #basic plot - full network or only subnetwork
  if (!is.null(plot.subnetwork)) {
    
    plot_viz <- ggplot(data = network[network$cluster_color == plot.subnetwork, ],
                       aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "lightgrey") + 
      theme_blank()
    
  } else {
    
    plot_viz <- ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_edges(color = "lightgrey") + 
      theme_blank()
    
  }
  
  
  
  
  # basic coloring
  if (color.by=="basic"){
    plot_viz <- plot_viz + 
      geom_nodes(fill="#f0f0f0",color="grey",shape=21,size = 2)
  }
  

  
  
  # cluster coloring
  if (color.by=="cluster") {
    if (!is.null(select.cluster)) {
      
      plot_viz <- plot_viz +
        geom_nodes(data = network[network$cluster_color != select.cluster, ],
                   fill="#f0f0f0",color="grey",shape=21,size = 2) +
        geom_nodes(data = network[network$cluster_color %in% select.cluster, ],
                   aes(fill = cluster_color),color="grey",shape=21,size = 3) +
        ggtitle("Network colored by clusters") + 
        scale_fill_manual(values = cluster.col,
                          name = "",
                          breaks = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_color"]))),
                          labels = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_color"]))))
      
    } else {
      
      plot_viz <- plot_viz +
        geom_nodes(aes(fill = cluster_color),color="grey",shape=21,size = 2) +
        scale_fill_manual(values = cluster.col,
                          name = "",
                          breaks = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_color"]))),
                          labels = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_color"]))))
    }
    
  }
  
  
  
  # GFC coloring
  if (grepl("GFC", color.by)){
    
    plot_viz <- plot_viz + 
      geom_nodes(aes_string(fill = network[,color.by]),color="grey",shape=21,size = 2) +
      scale_fill_gradientn(colours = colorRampPalette(c("#1b15c2", "#ffffff", "#ffbb00"))(11), 
                           limits = c(-range_GFC,range_GFC)) +
      ggtitle(color.by) + 
      labs(fill = "GFC")
    
  }
  
  
  
  # TF coloring
  if (color.by=="TF"){
    
    plot_viz <- plot_viz +
      geom_nodes(fill = "#f0f0f0",color = "grey",shape=21,size = 2) +
      geom_nodes(data=network[network$TF_info != "noTF", ],
                 aes(x = x, y = y, fill = TF_info), shape = 21, color = "grey", size = 3) +
      scale_fill_manual(values = c("Chromatin_remodeller" = "#0073C2FF", "Co_factor" = "#EFC000FF",
                                   "RNBP" = "#79AF97FF", "TF" = "#CD534CFF")) + 
      ggtitle("Network colored by TF")
  }
  

  
  # GSEA coloring
  if (grepl("^HALLMARK_|^GO_|^KEGG_|^REACTOME_",x = color.by)){
    
    column_to_choose <- stringr::str_split(color.by, pattern = "_",simplify = T)[,1]
    
    plot_viz <- plot_viz +
      geom_nodes(fill = "#f0f0f0",color = "grey",shape=21,size = 2) +
      geom_nodes(data=network[grepl(color.by, x = network[[column_to_choose]]), ],
                 aes(x = x, y = y),
                 shape=21, color = "grey", fill = "#e31a1c", size = 3) +
      ggtitle(color.by)
    
  }
  
####  
  


  
  
  # Add label
  if (!is.null(gene.label)){

    if(color.by=="cluster") {

      plot_viz <- plot_viz +
       
        geom_nodetext_repel(data = network[network$Gene %in% gene.label,],
                            aes(x = x, y = y, label = Gene),  size = 5)+ #nudge_x = 50,
        geom_nodes(data = network[network$Gene %in% gene.label, ],
                   aes(x = x, y = y, fill = cluster_color),
                   shape=21,color="black", size = 3, stroke=1, show.legend = F) +
        scale_fill_manual(values = cluster.col)#,
                          # name = "",
                          # breaks = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_name"]))),
                          # labels = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_name"]))))


    } else if(grepl("GFC", color.by)) {

      plot_viz <- plot_viz +
       
        geom_nodetext_repel(data = network[network$Gene %in% gene.label,],
                            aes(x = x, y = y, label = Gene) ,nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Gene %in% gene.label,],
                   aes(x = x, y = y, fill=network[network$Gene %in% gene.label, color.by]),
                   shape=21,color="black", size = 3,stroke=1, show.legend = F) +
        scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), 
                             limits = c(-range_GFC,range_GFC))
      
      
    } else if(color.by=="TF") {
      
      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Gene %in% gene.label,],
                            aes(x = x, y = y, label = Gene) ,nudge_x = 50, size = 5) + 
        geom_nodes(data= network[network$Gene %in% gene.label,],
                   aes(x = x, y = y,fill = TF_info), 
                   shape = 21, color = "black", size = 3, stroke=1, show.legend = F) +
        scale_fill_manual(values = c("Chromatin_remodeller" = "#0073C2FF", "Co_factor" = "#EFC000FF",
                                     "RNBP" = "#79AF97FF", "TF" = "#CD534CFF", "noTF" = "#f0f0f0"))
      
      
    } else if(grepl("^HALLMARK_|^GO_|^KEGG_|^REACTOME_",x = color.by)) {
      
      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Gene %in% gene.label, ],
                            aes(x = x, y = y, label = Gene) ,nudge_x = 50, size = 5) +
        geom_nodes(data=network[network$Gene %in% gene.label, ],
                   aes(x = x, y = y),
                   shape=21, color = "black", fill = "#e31a1c", size = 3) 
      
      
      
    } else {

      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Gene %in% gene.label,],
                            aes(x = x, y = y, label = Gene) ,nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Gene %in% gene.label,],
                   aes(x = x, y = y),shape=21,fill=NA,color="black",
                   size = 4,stroke=1)

    }


  }
  
  
  ## add class column
  
  network$Class <- "nope"
  network$Class[grepl(pattern = "PC", x = network$Gene)] <- "PCHO"
  network$Class[grepl(pattern = "CE", x = network$Gene)] <- "CHES"
  network$Class[grepl(pattern = "Cer", x = network$Gene)] <- "CERA"
  network$Class[grepl(pattern = "CL", x = network$Gene)] <- "CARD"
  network$Class[grepl(pattern = "DAG", x = network$Gene)] <- "DACG"
  network$Class[grepl(pattern = "HexCer", x = network$Gene)] <- "HEXC"
  network$Class[grepl(pattern = "DiHexCer", x = network$Gene)] <- "LACC"
  network$Class[grepl(pattern = "LPC", x = network$Gene)] <- "LPCH"
  network$Class[grepl(pattern = "PC-O", x = network$Gene)] <- "PC-O"
  network$Class[grepl(pattern = "MAG", x = network$Gene)] <- "MACG"
  network$Class[grepl(pattern = "LPC-O", x = network$Gene)] <- "LPCO"
  network$Class[grepl(pattern = "PE", x = network$Gene)] <- "PETA"
  network$Class[grepl(pattern = "PE-O", x = network$Gene)] <- "PE-O"
  network$Class[grepl(pattern = "PI", x = network$Gene)] <- "PHIN"
  network$Class[grepl(pattern = "PS", x = network$Gene)] <- "PHSE"
  network$Class[grepl(pattern = "SM", x = network$Gene)] <- "SPMY"
  network$Class[grepl(pattern = "TAG", x = network$Gene)] <- "TACG"
  

  ## add class label
  
  if (!is.null(class.label)){
    
    if(color.by=="cluster") {
      
      plot_viz <- plot_viz +
        
        # geom_nodetext_repel(data = network[network$Class %in% class.label,],
        #                     aes(x = x, y = y, label = Class), nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Class %in% class.label, ],
                   aes(x = x, y = y, fill = cluster_color),
                   shape=21,color="black", size = 3, stroke=1, show.legend = F)+
        # scale_fill_manual(values = cluster.col,
        #                   name = "",
        #                   breaks = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_name"]))),
        #                   labels = gtools::mixedsort(as.character(unique(network[network$cluster_color != "white", "cluster_name"]))))+
        ggtitle(paste("Network colored by Cluster and highlithed for Class", paste(noquote(class.label),collapse=', ')))

      
    } else if(grepl("GFC", color.by)) {
      
      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Class %in% class.label,],
                            aes(x = x, y = y, label = Class), nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Class %in% class.label, ],
                   aes(x = x, y = y, fill=network[network$Gene %in% gene.label, color.by]),
                   shape=21,color="black", size = 3,stroke=1, show.legend = F) +
        scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), 
                             limits = c(-range_GFC,range_GFC))
      
      
    } else if(color.by=="TF") {
      
      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Class %in% class.label,],
                            aes(x = x, y = y, label = Class), nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Class %in% class.label, ],
                   aes(x = x, y = y,fill = TF_info), 
                   shape = 21, color = "black", size = 3, stroke=1, show.legend = F) +
        scale_fill_manual(values = c("Chromatin_remodeller" = "#0073C2FF", "Co_factor" = "#EFC000FF",
                                     "RNBP" = "#79AF97FF", "TF" = "#CD534CFF", "noTF" = "#f0f0f0"))
      
      
    } else if(grepl("^HALLMARK_|^GO_|^KEGG_|^REACTOME_",x = color.by)) {
      
      plot_viz <- plot_viz +
        
        geom_nodetext_repel(data = network[network$Class %in% class.label,],
                            aes(x = x, y = y, label = Class), nudge_x = 50, size = 5)+
        geom_nodes(data = network[network$Class %in% class.label, ],
                   aes(x = x, y = y),
                   shape=21, color = "black", fill = "#e31a1c", size = 3) 
      
      
   
      
    } else {
      
      plot_viz <- plot_viz +
        
        
        geom_nodes(data = network[network$Class %in% class.label, ],
                   aes(x = x, y = y,fill=Class),shape=21,color="black",
                   size = 4,stroke=1)
      
    }
    
    
  }
  
  
  
  #Save pdf
  
  if (save.pdf){
    
    if(length(gene.label) < 6) {
      vec <- c("Network", use.layout, select.cluster, plot.subnetwork, color.by, gene.label)
      vec <- vec[vec!=""]
     
    } else {
      vec <- c("Network", use.layout, select.cluster, plot.subnetwork, color.by, "multiple_genes_labeled")
      vec <- vec[vec!=""]
    }
    
    ggsave(filename = paste0(paste(vec,collapse = "_"),".pdf"), path = paste0(working_directory,cutoff_wd),
           plot = plot_viz, width = 6, height = 5, units = "in", device = cairo_pdf)
    
    return(plot_viz)

  } else {

    return(plot_viz)

   }
  

}





#GFC plots ------------------------------------------------------------------



GFC_colored_network <- function(network,
                                select.cluster,
                                plot.subnetwork,
                                gene.label,
                                use.layout,
                                save.pdf,
                                save.single.pdf){
  
  names(network) <- gsub(pattern = "_2", replacement = " GOLD2", x = names(network))
  names(network) <- gsub(pattern = "_3_4", replacement = " GOLD3/4", x = names(network))
  names(network) <- gsub(pattern = "_ctrl", replacement = " Control", x = names(network))
  
  GFC_names <- colnames(network[grepl("GFC",colnames(network))])
  
  
  plot_list <- lapply(GFC_names, function(x){
    plot_GFC <- visualize_network(network = network,
                                  color.by = x,
                                  select.cluster = select.cluster,
                                  use.layout = use.layout,
                                  plot.subnetwork = plot.subnetwork,
                                  gene.label = gene.label,
                                  save.pdf=save.pdf)
  })
  
  
  if(save.single.pdf) {
    plots <- gridExtra::marrangeGrob(grobs = plot_list, ncol = 1, nrow = 1)  #extra step prevents empty page in pdf file
    
    cairo_pdf(filename = paste0(working_directory, cutoff_wd, "/", 
                                paste("Network_colored_by_GFC", use.layout, sep = "_"), ".pdf"), 
              onefile = T)
    print(plots)
    dev.off()
    dev.off()
  } 
  
  return(plot_list)
  
}



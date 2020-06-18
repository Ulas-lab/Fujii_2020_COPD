heatmap_clusters <- function(data) {
                              
  # Data preparation
  
   mat_heatmap= data %>% 
    dplyr::filter(cluster_included=="yes")%>% 
    dplyr::select(grp_means) %>% 
    tidyr::separate(col=grp_means, sep=",",
                    convert = T,
                    into = 1:(stri_count_regex(data$grp_means[1], pattern = ",")+1) %>% 
                      as.character()) %>% as.matrix()
  
  rownames(mat_heatmap) =data %>% 
    dplyr::filter(cluster_included=="yes") %>% 
    dplyr::mutate(labs=color) %>% dplyr::pull(labs) 
  
  colnames(mat_heatmap) =  data %>% 
    dplyr::filter(cluster_included=="yes") %>% dplyr::select(conditions) %>% purrr::map(1) %>%
    stri_split_regex(pattern = "#") %>% unlist()
  
  
  # # Annotation
  # 
  # annotation = as.data.frame(data %>% dplyr::filter(cluster_included=="yes"))
  # annotation = annotation[ ,c("clusters", "color")]
  # rownames(annotation) <- rownames(mat_heatmap)
  # 
  # annotation_names <- annotation[ , "clusters", drop = F]
  # colnames(annotation_names) <- "ClusterColor"
  # 
  # annotation_colors <- annotation[ , "color"]
  # names(annotation_colors) <- annotation_names$ClusterColor
  # annotation_colors <- list(ClusterColor = annotation_colors)
  # 
  # 
  # 
  # # Plot heatmap 
  # 
  #  hm <- pheatmap::pheatmap(mat = mat_heatmap,
  #                          color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
  #                          cellwidth = 30,
  #                          cellheight = 30,
  #                          cluster_rows =T,
  #                          cluster_cols = T,
  #                          clustering_distance_rows = "euclidean",
  #                          clustering_distance_cols = "euclidean",
  #                          clustering_method = "complete",
  #                          fontsize = 10,
  #                          treeheight_row = 25,
  #                          treeheight_col = 25,
  #                          legend_labels = "GFC",
  #                          annotation_row = annotation_names, 
  #                          annotation_legend = FALSE,
  #                          annotation_colors = annotation_colors,
  #                          angle_col = 90,
  #                          #cutree_rows = 2,
  #                          #cutree_cols = 2,
  #                          filename = "Heatmap_clusters.pdf")

  
  colnames(mat_heatmap) <- c("GOLD2", "GOLD3/4", "Control")
  
  
  #Annotation
  
  ann_row <- as.data.frame(data[data$color != "white", c("color", "gene_no")])
  colnames(ann_row) <- c("Cluster", "Gene_count")
  rownames(ann_row) <- ann_row$Cluster
  
  gene_counts <- ann_row$Gene_count
  names(gene_counts) <- ann_row$Cluster
  
  clusters <- ann_row$Cluster
  names(clusters) <- ann_row$Cluster
  
  rowAnno <- ComplexHeatmap::rowAnnotation("Species counts" = anno_barplot(x = gene_counts, width = unit(3, "cm")),
                                           "Clusters" = anno_simple(x = clusters, col = clusters,
                                                                  simple_anno_size = unit(0.5, "cm") ),
                                           width = unit(4, "cm"), annotation_name_side = "top", 
                                           gap = unit(2, "mm"), annotation_name_rot = 0,
                                           annotation_name_gp = gpar(fontsize = 8))
  
  breakList <- seq(-2, 2, by = .1)
  
  #cairo_pdf(filename, width = 9, height = 8, pointsize = 9)
  
  pdf(file = paste0(working_directory,cutoff_wd, "/Heatmap_clusters.pdf"), width = 5, height = 4, 
      pointsize = 7)
  
  hm_cluster <- ComplexHeatmap::Heatmap(matrix = mat_heatmap,
                         # col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(breakList)),
                          col = colorRampPalette(c("#1b15c2", "#ffffff", "#ffbb00"))(length(breakList)),
               # Seurat colors   col = colorRampPalette(PurpleAndYellow())(length(breakList)), 
                         clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean",
                          clustering_method_rows = "complete",
                          clustering_method_columns = "complete",
                          right_annotation = rowAnno,
                          column_names_rot = 360,
                          column_names_centered = T,
                          row_names_gp = gpar(fontsize = 10),
                          column_names_gp = gpar(fontsize = 10),
                          rect_gp = grid::gpar(col = "black"), 
                          heatmap_legend_param = list(title = "", legend_height = unit(3, "cm")),
               border = T, 
               column_order = c("Control", "GOLD2", "GOLD3/4"))

  print(hm_cluster)
  dev.off()
 
  
  print(hm_cluster)
  #gb = grid.grabExpr(draw(heatmap))
  # pdf(paste0(name_for_pdf,".pdf") , width = 8 , height = 8)
  # grid.arrange(gb,ncol=1)
  #   Sys.sleep(5)
}

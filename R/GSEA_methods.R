library(msigdbr)
library(fgsea)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(stringr)

############################################################
#Choose/make your gene sets for gene set enrichment analysis
############################################################

#' generate_gene_set_object
#'
#' This function generates a list of gene sets from msigdbr, which are formatted for GSEA.
#'
#' @param species Choose species = 'human' or 'mouse'. Default to human.
#' @param gene_sets Choose your gene sets from the list of: "Hallmark", "GO:BP", "Reactome", "GO:MF", "TFT-val", "TFT-pred", or "Cell types". Default = "Hallmark" and "Reactome"
#' @return A list of gene sets by name, with their associated genes in each list element.
#' @import msigdbr
#' @export

generate_gene_set_object <- function(species = "human", gene_sets = c("Hallmark", "Reactome")){
  if(species == "human"){
    sp <- "Homo sapiens"
  }else{
    sp <- "mouse"
  }
  
  gene_set_df <- data.frame("Selection" = c("Hallmark", "GO:BP", "Reactome", "GO:MF", "TFT-val", "TFT-pred", "Cell types"),
                          "Cat" = c("H", "C5", "C2", "C5", "C3", "C3", "C8"),
                          "subCat" = c("", "GO:BP", "CP:REACTOME", "GO:MF", "TFT:TFT_Legacy", "TFT:GTRD", ""))
  
  if(any(!(gene_sets %in% gene_set_df$Selection))){
    stop("Sorry, but you have at least one gene set mismatch. Please select from: ", paste0(gene_set_df$Selection, sep = " ,"))
  }
  
  gene_set_df <- gene_set_df[match(gene_sets, gene_set_df$Selection),]
  
  
  ##START here
  gene_set_list <- list()
  
  if (nrow(gene_set_df)>1){
    for (i in 1:nrow(gene_set_df)){
      if (gene_set_df$subCat[i] == ""){
        gene_set_list[[i]] <- msigdbr(species = sp, category = gene_set_df$Cat[i])
      }else{
        gene_set_list[[i]] <- msigdbr(species = sp, category = gene_set_df$Cat[i], subcategory = gene_set_df$subCat[i])
      }
    }
  }else{
    if (gene_set_df$subCat == ""){
      gene_set_df <- msigdbr(species = sp, category = gene_set_df$Cat)
    }else{
      gene_set_df <- msigdbr(species = sp, category = gene_set_df$Cat, subcategory = gene_set_df$subCat)
    }
  }
  
  #Condense if larger than 1 gene set
  if (length(gene_sets) > 1){
    gene_set_df <- do.call(rbind, gene_set_list)
  }
  
  msigdbr_list = split(x = gene_set_df$gene_symbol, f = gene_set_df$gs_name)
  
  return(msigdbr_list)
}




#' scores_fun
#'
#' This function generates a named, ordered vector of scores from your differential expression result. There are multiple ways to calcualte this and this function has two.
#'
#' @param difex_res A data frame of a differential expression result. The difex_res must have gene names as the rownames.
#' @param method Choose your method of calculating the score 'Effect.size.only' will use things like logFC or attention differential and 'p.value.adjusted' will combine effect size and p.value.
#' @param effect.size.col The name of your effect size column.
#' @param sig.col The name of the column you want to use for significance test.
#' @return An ordered vector of gene scores with gene names.
#' @export

############################################################
#Generate a scores object for GSEA
############################################################
#Make a score list
scores_fun <- function(difex_res, method, effect.size.col, sig.col){
  methods <- c("Effect.size.only", "p.value.adjusted")
  
  if(!(method %in% methods)){
    stop("Please choose from the following methods: 'Effect.size.only' or 'p.value.adjusted'")
  }
  
  if(method == "Effect.size.only"){
    scores <- as.vector(difex_res[,effect.size.col])
  }
  
  if(method == "p.value.adjusted"){
    #Correct for 0 values
    p.vals <- difex_res[,sig.col]
    p.vals[p.vals == 0] <- min(p.vals[p.vals > 0])
    
    #Generate the score
    scores <- as.vector(difex_res[,effect.size.col]*(-log(p.vals)))
  }
  
  #Change the names of the elements
  names(scores) <- rownames(difex_res)
  
  #Remove NA values and any duplicates
  scores <- scores[!is.na(names(scores))]
  scores <- scores[!duplicated(names(scores))]
  
  #Reorder the scores
  scores <- sort(scores, decreasing = TRUE)
  
  ######Quick tests
  if(length(scores) < 1000){
    warning("Warning!!!! You have very few genes in this result. Reconsider the test and the differential expression method.")
  }
  
  if(all(scores > 0)){
    warning("Warning!!!! You have an all positive score set, make sure you adjust the GSEA accordingly!")
  }
  if(all(scores < 0)){
    warning("Warning!!!! You have an all negative score set, make sure you adjust the GSEA accordingly!")
  }
  if (any(scores > 0) & any(scores < 0)){
    if(!(all(scores[1:(0.1*floor(length(scores)))] > 0)) | !(all(scores[floor(0.9*length(scores)):length(scores)] < 0))){
      warning("Warning!!!! Your scores are both positive and negative but > 90% of them skew one direction. Please adjust your scoring.")
    }
  }
  ######
  #There's room to build in more safety checks here later
  
  #######
  
  return(scores)
}




############################################################
#Run the GSEA
############################################################

#' gsea_run
#'
#' This function runs gene set enrichment analysis on a gene scores vector generated from differential expression.
#'
#' @param gene_set_list A list of gene sets where the elements are named and each element contains a vector of gene names.
#' @param scores_list A vector of gene scores, ordered from high to low. The vector elements must be named.
#' @param minimum_GS_size The minimum number of overlapping genes a gene set must have to be tested. Default = 30.
#' @param maximum_GS_size The maximum number of overlapping genes a gene set can have to be tested. Default = 250.
#' @return A data frame of gene set enrichment results. This result is unfiltered for significance.
#' @import fgsea
#' @export

gsea_run <- function(gene_set_list, scores_list, minimum_GS_size = 30, maximum_GS_size = 250){

  fgseaRes <- fgsea(pathways = gene_set_list, 
                    stats = scores, 
                    maxSize=maximum_GS_size, minSize = minimum_GS_size,
                    scoreType = "std", #If using attention scores, make sure to change to "pos" if all attention scores are non-negative
                    eps = 0)
  
  fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing = TRUE),]
  
  return(fgseaRes)
}



############################################################
#GSEA Plot
############################################################
#' gsea_plot
#'
#' This function generates a standard gsea plot where the area under the curve is calculated.
#'
#' @param GSEA_result A data frame of results from gene set enrichment analysis. This must have been generated from fgsea.
#' @param gene_set_list A list of gene sets used for the GSEA which generated the GSEA_result.
#' @param scores The scores vector used for GSEA.
#' @param row_number The row number in the GSEA_result data frame, which you would like to plot.
#' @return A standard GSEA plot of the normalized enrichment of a single gene set and the scores that generated the result.
#' @import msigdbr
#' @import cowplot
#' @import ggplot2
#' @export

gse_plot <- function(GSEA_result, gene_set_list, scores, row_number){
  Gene_set_name <- GSEA_result$pathway[row_number]
  
  enrichment_plot <- plotEnrichment(gene_set_list[[Gene_set_name]], scores, gseaParam = 1, ticksSize = 0.2) + 
    labs(title = Gene_set_name, subtitle = paste0("NES = ", round(GSEA_result$NES[row_number], digits = 2), ",  Adjusted pval = ", formatC(GSEA_result$padj[row_number], format = "e", digits = 2))) +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 16, color = "black"))
  
  score_plot <- data.frame("Score" = "", 
                           "Genes" = factor(names(scores), levels = names(scores)), 
                           "Rank" = scores) %>%
    ggplot(aes(x = Genes, y = Rank)) +
    geom_bar(stat = "identity") +
    labs(y = "Score") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.title.y=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none",
          plot.margin=margin(t=-0.3,unit="cm"))
  
  
  plot_grid(enrichment_plot, score_plot, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(3,1), vjust = 0, greedy = TRUE)
  
}


############################################################
#Generate a barplot
############################################################
#' gsea_barplot
#'
#' This function generates a barplot of specified size to display multiple gene set enrichment pathways and their significance.
#'
#' @param gsea_res A data frame of results from gene set enrichment analysis. This must have been generated from fgsea.
#' @param number_up The number of +NES gene sets you want to display.
#' @param number_down The number of -NES gene sets you want to display.
#' @return A ggplot object with bars representing the normalized enrichemnt score and colored by the -log(p.val). The number of +NES and -NES is selected in order from highest abs(NES) to lowest.
#' @import stringr
#' @import ggplot2
#' @export
gsea_barplot <- function(gsea_res, number_up = 10, number_down = 10){
  
  # Subset the data for top up and down pathways
  gsea_res <- gsea_res[c(1:number_up, (nrow(gsea_res)-number_down+1):nrow(gsea_res)), ]
  
  #Sub the underscore
  gsea_res$pathway <- sapply(gsea_res$pathway, function(x) gsub("_", " ", x))
  
  # Convert pathway column to a factor to control the order of the bars
  gsea_res$pathway <- factor(gsea_res$pathway, levels = rev(gsea_res$pathway))
  
  # Plot with continuous fill based on -log(padj) using the YlOrRd palette
  ggplot(gsea_res) + 
    aes(x = NES, y = pathway, fill = -log(padj)) +  # Use -log(padj) for fill
    geom_bar(stat = "identity", width = 0.7) + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    geom_vline(xintercept = 0, linetype="solid", color = "black", linewidth=0.25) +
    labs(x = "NES", y = "") + 
    theme_bw() +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "-log(padj)", limits = c(0, max(-log(gsea_res$padj), na.rm = TRUE)))  # Use YlOrRd color scale
}








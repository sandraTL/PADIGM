

#' function that output a barplot graph related to one specific gene with all
#' the shortest distances from that gene to all metabolites
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#' for param gene : is a gene in data ex: hsa:8801
#' @param pathwayId, data(gene, metabolites), gene
#' @keywords  Graph, barplot, shortestDistance
#' @export
#' @examples barplotFunction(hsa01100, associatedGeneMetaDF,
#'  completeMetaboliteDF, gene)

barplotFunction <- function(pathwayId, associatedGeneMetaDF,
                            completeMetaboliteDF, gene){
    ############################################################
    #' Serious need to refactor this function in multiple ones
    #' It is way to long
    ############################################################

    geneCommonName <- getCommonNames(c(gene), "gene")


    # get all shortest oaths from data entry
    shortestsPathsDF <- getDistanceAll(pathwayId,
                     associatedGeneMetaDF[associatedGeneMetaDF$gene == gene,],
                     completeMetaboliteDF  );

    associatedMetabo <- getAssociatedMetaboByGene(associatedGeneMetaDF,gene)

    aM <- data.frame(associatedMetabo);

    shortestsPathsDF <- t(shortestsPathsDF);
    shortestsPathsDF <- data.frame(shortestsPathsDF);
    gene1 <- gsub(":", ".", gene);

    # add metabolite row
    shortestsPathsDF[ "metabolites" ] <- rownames(shortestsPathsDF);
    # print(shortestsPathsDF)
    # get a subset of shortestsPathsDF contaning only geneOf interest, gene
    # and metaboltie column
    spDF.sub2 <- subset(shortestsPathsDF ,select = c(gene1, "metabolites"))
    spDF.sub2[is.na(spDF.sub2)] <- Inf;

    print(spDF.sub2)

     infVal <- which(spDF.sub2[,1] == Inf)
     temp <- spDF.sub2
     temp[infVal,] <- -1
     maxVal <- max(temp[,1])
#
#     ## get frequency of every value until the maxVal found + Inf val
     frequenceDistDF <- data.frame(table(factor(spDF.sub2[,1],
                                            levels=c(0:maxVal,Inf))))
     print(frequenceDistDF)
     colnames(frequenceDistDF) <- c("Var1", "Freq")


   # frequenceDistDF <- table(spDF.sub2[,gene1]);
   # frequenceDistDF <- data.frame(frequenceDistDF)

    #' initiation of values
    test = FALSE;
    results <- data.frame();
    #' creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(spDF.sub2)){
        test = FALSE;
        for(row2 in 1:nrow(aM)){

            if(spDF.sub2[row1,"metabolites"] == aM[row2,]){
                test<- TRUE;

                break;
            }else test<- FALSE;
        }

        results <- rbind(results,test);
        colnames(results) <- c("Associations")
        return <- results;
    }
    shortestsPathsDF.plot<-spDF.sub2;
    shortestsPathsDF.plot<-cbind(shortestsPathsDF.plot, Associations = results);

    #' initiation of values
    test = FALSE;
    results1 <- data.frame();

    #' creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(frequenceDistDF)){
        test = FALSE;
        for(row2 in 1:nrow(shortestsPathsDF.plot)){

            if(frequenceDistDF[row1,"Var1"] == shortestsPathsDF.plot[row2,gene1]){
                if(shortestsPathsDF.plot[row2,"Associations"] == TRUE){
                    test<- TRUE;

                    break;
                }
            }else test<- FALSE;
        }

        results1 <- rbind(results1,test);
        colnames(results1) <- c("Associations")
        return <- results1;
    }


    # Add a column for the coloring of the bar associated with gene to subgraph
    shortestsPathsDF.plot<-frequenceDistDF;
    shortestsPathsDF.plot<-cbind(shortestsPathsDF.plot, Associations = results1);

    # create barplot
    print(shortestsPathsDF.plot)
    plot <- ggplot2::ggplot(shortestsPathsDF.plot, ggplot2::aes(
        x = Var1,
        y = Freq,
        fill = Associations
    ),environment = environment())

    plot <- (plot + ggplot2::geom_bar(stat="identity", colour="black",
                                      size=0.5, position="identity",width=1)
             + ggplot2::theme_bw()
             + ggplot2::theme(panel.border = ggplot2::element_blank(),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank(),
                              text = ggplot2::element_text(size=12, family="Arial"),
                              axis.line = ggplot2::element_line(colour = "black"))
             + ggplot2::xlab("Distance from Gene")
             + ggplot2::ylab("Metabolite count")
             + ggplot2::ggtitle(geneCommonName)
            + ggplot2::scale_y_sqrt()
        #   +  ggplot2::scale_y_continuous(trans=log_trans(), breaks=c(0,10,20,50),expand = c(0,0))
#              + ggplot2::scale_y_continuous( breaks =trans_breaks("log2",
#                                             function(x) 2^x)
#                                             ,expand = c(0,0))

 #           + ggplot2::scale_y_continuous(expand = c(0,0), trans = 'log2')

             + ggplot2::scale_fill_manual(values = c("FALSE" ="grey",
                                                     "TRUE" = "red3"))
             + ggplot2::guides(fill=FALSE)
    );
    filename = paste0(geneCommonName,".png")
    ggplot2::ggsave(filename, dpi=300)
    # print plot it could change to save the graph image somewhere
    print(plot);

}

getAssociatedMetaboByGene <- function(data, gene){

    selectedRows <- data[data$gene == gene,];
    associatedMetabo <- selectedRows[,2];

    return <- associatedMetabo;

}


##############################   IGRAPH S4   ###################################
#' object igraph instanciated as S4 classe to be able to use igraph objects in
#' other S4 classes (Graphe)
#'
setClass(Class = "igraph")


##########################   GRAPH ELEMENTS S4   ###############################

setClass("GraphElements" ,
         representation = representation(
             pathwayId = "character",
             nodeDF = "list",
             edgeDF = "list"
         ),  contains="list")



setGeneric("createIGraph", function(object)
{
    standardGeneric("createIGraph")
}
)

#' Creation of object igraph with data from XML file of pathway
#'
#' @param GraphElements
#' @keywords graphElements
#' @examples createIGraph(GraphElements)

setMethod("createIGraph", "GraphElements", function(object) {

    g <- igraph::graph.data.frame(object@edgeDF, directed=FALSE,
                                  vertices=object@nodeDF);
    return <- g;

})



##########################   GRAPHELEMENTS S4   ################################

#'S4 classe Graph, contains all information important to create graphe argument
#'and to calcul distances for it
setClass(
    Class= "Graph",
    representation = representation(
        graph = "igraph"
    ),
    contains=c("igraph", "GraphElements")
)



setGeneric("allShortestPaths", function(object, data)
{
    standardGeneric("allShortestPaths")
}
    )

#' fonction qui crée le grpah pour les calculs de distances, ne peut
#' pas etre une methode de la classe Graphe cause des erreurs de profondeurs
#' de calculs
#' @importFrom igraph shortest_paths
#' @param data frame des reactions extraites du KGML en listes
#' @keywords  kegg
#' @examples aLLShortestPaths(igraph, data)


setMethod("allShortestPaths","Graph", function(object, data){

        # Name column names by metabolite ids
         goalsMenu <- paste(data[,2], 1:length(data[,1]), sep="")
         output <- as.data.frame(matrix(rep(0, length(goalsMenu)), nrow=1))

         names(output) <- c(goalsMenu);


         # calcul al distances
         pl <-  apply(data, 1, function(x){

             dfTemp <- data.frame();

             dfTemp <- (igraph::shortest_paths(object@graph , x['mg'],
             as.vector(unlist(data[,c('m')]))))

        return <- dfTemp;
      })


       pl1 <- lapply(pl[], function(x){

           pl2 <- lapply(x$vpath, function(x){

               dfTemp1 <- data.frame();
               #' if length of path is 0 -> couldnt reach a path
               #' if length of path is >0 i have to do length -1
               #' for example path i am lookinf for path from 1523 to 1523
               #' path is 1523 and the length is 1. But the real length is 0.
               if(length(x[]) != 0){
               dfTemp1 <- union(dfTemp1, (length(x[])-1));
                }else if(length(x[]) == 0)
                    dfTemp1 <- union(dfTemp1, NA);
            return <- dfTemp1;
           })

            # create distance vector from 1 geneMetabolite to all metabolites
            output <- rbind(output, as.vector(unlist(pl2)));

        # delete the row of zeros
        output <- output[-c(1), ]
        return <- output;

      })
      # combine all vectors of distances
      df <- do.call(rbind.data.frame, pl1)

     return <- df;
})




    #' Fonction that calculates every shortest distances between each gene and
    #' all metabolites.
    #'
    #' @param data frame des reactions extraites du KGML en listes
    #' @keywords  kegg
    #' @export
    #' @examples getAllShortestDistances(igraph, data)
getAllShortestDistances <- function(pathwayId, data){

    # if doesnt exist yet
    if(isFileInDirectory(pathwayId) == FALSE){
    getPathwayKGML(pathwayId);
    }
#     download.file("rest.kegg.jp/get/hsa01100/kgml",
#                   paste(datahere,"hsa01100.xml",sep=""), method = "curl")
    #' create data frame for the vertices

    finalDF <- data.frame();

    graphe <-  createGraphFromPathway(pathwayId, data);

    finalDF <- getFinalDFSHortestDistance(graphe, data);

  #  barplotFunction(graphe,finalDF,"hsa:6999", data);
  #  print(finalDF);
    finalDF[is.na(finalDF)] <- Inf;
    return <- finalDF;


}

createGraphFromPathway <- function(pathwayId, data){

    # if doesnt exist yet
    # getPathwayKGML("hsa01100");
    #' create data frame for the vertices
    #'
    nodeDF <- getListNodeFromKGML(pathwayId);

    #' create data frame for edges
    edgeDF <- finalReactionEdgeDF(pathwayId);

    graphEl <- new("GraphElements", nodeDF= nodeDF,
                   edgeDF= edgeDF, pathwayId = pathwayId);

    igraphe <- createIGraph(graphEl);

    graphe <- new("Graph", graph = igraphe, graphEl);


    return <- graphe;


}










setGeneric("fromDFEntryToIGraphIdDF", function(object, data,
                                                         indiceMetabolite) {
    standardGeneric("fromDFEntryToIGraphIdDF");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' if indiceMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indiceMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  igraph, node, kgmlId
#' @examples fromDFEntryToIGraphIdDF(g, data, indiceMetabolite)

setMethod("fromDFEntryToIGraphIdDF", "Graph", function(object, data,
                                                    indiceMetabolite){
    #########################################################################
    ##### Ajouter une condition sur indiceMetabolite can only be 1 or 2 #####
    #########################################################################
    mgmDF <- data.frame();
    f <- apply(data,1, function(x){

        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);
        m2 <- getCompoundNodeKgmlId(object@graph, x[2], object@nodeDF);

        if( m1[1] != "na"){

            if(m2[1] != "na"){
                 mgmDF<- c(mg = m1[indiceMetabolite], m = m2);
                 return<- mgmDF;

            }else{mgmDF <- c(NA,NA)}  # si soit metabolite ou gene n'est
                                          # pas dans le graphe alors on ne peut
                                          # pas calculer la distance des gènes
                                          # associé
        }else{mgmDF <- c(NA,NA)}

    })
    f <- data.frame(f); ## typage
    f <- t(f)
    f <- data.frame(f); ## typage

    f[3]  <- unlist(data[1]); # add gene column
    f[4]  <- unlist(data[2]); # add metabolite column
    f <- f[rowSums(is.na(f)) != 2,] ## delete rows with 2  NA or more
    colnames(f) <- c("mg", "m","genes","metabolites") #set colnames
    return <- f;

})


setGeneric("getFinalDFSHortestDistance", function(object, data)
{
    standardGeneric("getFinalDFSHortestDistance");
}
)

#' A igraph function
#'
#'
#' Calls fromDFEntryToIGraphIdDF() twice to get :
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' where indiceMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' #' where indiceMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' This is because genes are edges and we calculate vertex to vertex
#' the distance between genes and metabolites.
#' From each gene we want to take the shortest distance between the 2 vertex
#' attached, so we calcul both and in an other function we take the smallest
#' and that will be the output.
#'
#' This function keeps, between the 2 metabolites of each gene vertexes,
#' the distance that is the smallest between the 2 to the same metabolite.
#'
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  igraph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, data, indiceMetabolite)

setMethod("getFinalDFSHortestDistance", "Graph", function(object, data){

    finalDF <- data.frame();

     # indiceMetabolite = 1 to get the first metabolite attached to gene
     m1g <- fromDFEntryToIGraphIdDF(object, data, 1);
     # indiceMetabolite = 2 to get the second metabolite attached to gene
     m2g <- fromDFEntryToIGraphIdDF(object, data, 2);

     #' get all shortest paths for both ends of gene to all metabolites
     r1 <- allShortestPaths(object, m1g);
     r2 <- allShortestPaths(object, m2g);

     #' take shotest paths between calculation of both end of a gene to the
     #' same metabolite
     for (row in 1:nrow(r1)) {

       r <- mergeVectorsLowerValues(r1[row,], r2[row,]);
       rf <- t(data.frame(r));
       finalDF <- rbind(finalDF,rf);

      return <- finalDF;
    }

     ## take out duplicated genes
     geneVect <- m1g[3][!duplicated(m1g[3]), ]
     geneVect <- as.vector(unlist(geneVect));

     # identification of colnames and rows names with metabolites and genes
     # frome entry list giving by user
    colnames(finalDF) <-  as.vector(unlist(m1g[4]));
    finalDF <- finalDF[!duplicated(finalDF), !duplicated(colnames(finalDF))]
    row.names(finalDF) <- geneVect;


   return <- finalDF;
})

#'function merge 2 numeric vectors and return 1 vector with the smalest
#'values for each position of the vectors
#'
#' @param 2 numeric vectors
#' @keywords  kegg
#' @examples d1 <- c(1,2,1,2,3,4), d2 <- c(2,3,1,1,5,2)
#'  mergeVectorsLowerValues(d1,d2) returns -> 1,2,1,1,3,2

mergeVectorsLowerValues <- function(A,B) {

    dataCombineDF <- rbind(A,B);
    dataCombineDF <- t(dataCombineDF);

        apply(dataCombineDF,1,function(x) {
            result <- data.frame();
            if(is.na(x[1]) && is.na(x[2])){
                result <- cbind(NA);
            }
            else if (x[1] <= x[2]) {
                 result <- cbind(x[1]);
            } else
                result <- cbind(x[2]);

            return <- result;
    })
}


#' Fonction to barplot all the shortest paths between 1 gene of interest to all
#' metabolites in the data value in entry
#'
#' @param keggpathwayId, data.frame() colunm c(gene, metabolites) representing association
#' between them for each row, gene of interest
#' @keywords  igraphe, barplot, shortestDistance
#' @export
#' @examples barplotFunction(hsa01100, data, gene)

barplotFunction <- function(pathwayId, data, gene){


    shortestsPathsDF <- getAllShortestDistances(pathwayId, data);
    #graphe <-  createGraphFromPathway(pathwayId, data);

    associatedMetabo <- getAssociatedMetaboByGene(data,gene)
    aM <- data.frame(associatedMetabo);

    shortestsPathsDF <- t(shortestsPathsDF);
    shortestsPathsDF <- data.frame(shortestsPathsDF);
    gene1 <- gsub(":", ".", gene);

    # add metabolite row
    shortestsPathsDF[ "metabolites" ] <- rownames(shortestsPathsDF);

    # get a subset of shortestsPathsDF contaning only geneOf interest, gene
    # and metaboltie column
    spDF.sub2 <- subset(shortestsPathsDF ,select = c(gene1, "metabolites"))
    spDF.sub2[is.na(spDF.sub2)] <- Inf;


    frequenceDistDF <- table(spDF.sub2[,gene1]);
    frequenceDistDF <- data.frame(frequenceDistDF)


   # initiation of values
    test = FALSE;
    results <- data.frame();
    # creation of vector to fill bar colors automatically
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


    #frequenceDistDF : var1=distances Freq=frequenceOfdistance
    #shortestPathsDF.plot : dataframe with gene - metablite - true or false associations

    # initiation of values
    test = FALSE;
    results1 <- data.frame();
    # creation of vector to fill bar colors automatically
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





    plot <- ggplot2::ggplot(shortestsPathsDF.plot, ggplot2::aes(
        x = Var1,
        y = Freq,
        fill = Associations
      #  color = Associations
        ),environment = environment())

    plot <- (plot + ggplot2::geom_bar(stat="identity"
                                      #,position = "dodge"
                                      )
      + ggplot2::theme_bw()
      + ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"))
      + ggplot2::xlab("Distance from Gene")
      + ggplot2::ylab("Metabolite count")
      + ggplot2::ggtitle(gene)
    #  + ggplot2::scale_color_manual(values = c("FALSE" ="grey", "TRUE" = "red"))
      + ggplot2::scale_fill_manual(values = c("FALSE" ="grey", "TRUE" = "red"))
      + ggplot2::guides(fill=FALSE)
     );







   # print plot it could change to save the graph image somewhere
   print(plot);

}

#' HeatMap function, is using the result data.frame from getAllShortestDistances
#' to create it. It uses the associations data form the users.
#'
#' @param keggpathwayId, data.frame() colunm c(gene, metabolites)
#' representing association between them for each row, gene of interest
#' @keywords  igraphe, heatmap, shortestDistance
#' @export
#' @examples heatmapFunction(hsa01100, data)

heatmapFunction <- function(pathwayId, data){

    AllSP <- getAllShortestDistances(pathwayId, data);

    dat <- data.frame( Row = rep(row.names(b15), each= length(colnames(b15))),
                       Col = rep(colnames(b15), times= length(row.names(b15))),
                       Y = c(t(b15)),
                       Associations = rep(c(TRUE,FALSE,FALSE,
                                    FALSE,FALSE,FALSE,FALSE),times=9));



    # create frame to color edges of associated genes and metabolites
    frames = dat[dat$Associations, c("Row", "Col")]

    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)


    p2 = ggplot2::ggplot(data=dat) +
        ggplot2::geom_raster(ggplot2::aes(x=Row, y=Col, fill=Y)) +
        ggplot2::theme(
               panel.border = ggplot2::element_blank(),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank()) +

        ggplot2::scale_fill_gradient2(low="yellow",mid="orange", high="darkred", na.value="lightgrey", name="") +
        ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
            ggplot2::aes(xmin=Row - 0.5, xmax=Row + 0.5, ymin=Col - 0.5, ymax=Col + 0.5)) +
        ggplot2::geom_text(fill = dat$Y, label = round(dat$Y, 1), ggplot2::aes(x = Row, y = Col)) +
        ggplot2::labs(title="Heatmap")


print(p2);


}


getAssociatedMetaboByGene <- function(data, gene){

    selectedRows <- data[data$gene == gene,];
    associatedMetabo <- selectedRows[,"metabolites"];

    return <- associatedMetabo;

}







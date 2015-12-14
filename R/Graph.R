
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
#' of interest
#' @param GraphElements
#' @keywords graphElements
#' @examples createIGraph(GraphElements)

setMethod("createIGraph", "GraphElements", function(object) {

   # print("createIGraph");
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

setGeneric("allShortestPaths", function(object, associatedGeneMetaDF,
                                        completeMetaboliteDF)
{
    standardGeneric("allShortestPaths")
}
    )

#' function that uses the object Graph and calculs distances for every pair of
#' gene - metabolite from data.
#'
#' Note that for 1 gene is attached 1 enzyme which is an edge in our graph.
#' In order to calcul a distance from vertice to vertice this function is
#' called twice in function 'getFinalDFSHortestDistance'.
#'
#' First call is to calcul all the distance for the one vertice of the gene
#' to the all metabolites, the second call is to calcul the distance
#' from the other vertice to all metabolites and then
#' choose the smallest distance between the 2.
#'
#' The selection of the first vertice and the second vertice is done in a other
#' function 'getIdGeneInGraph'
#'
#' param data where
#'  mg = id (from igraph in object Graph) of 1 metabolite related to the gene
#'  m = id (igraph in object Graph)
#'  genes = hsa:... id from KEGG
#'  metabolites = C.... id from KEGG
#'
#' @importFrom igraph shortest_paths
#' @param object Graph, data(mg, m, genes, metabolites)
#' @keywords  kegg
#' @examples allShortestPaths(Graph, data)


setMethod("allShortestPaths","Graph", function(object, associatedGeneMetaDF,
                                               completeMetaboliteDF){
     # print("allShortestPaths")


        #'Handling dataFrame
        #'Name column names by metabolite ids
         goalsMenu <- paste(completeMetaboliteDF[,2], sep="")
         output <- as.data.frame(matrix(rep(0, length(goalsMenu)), nrow=1))
         names(output) <- c(goalsMenu);


         #'calcul al distances
         pl <-  apply( associatedGeneMetaDF, 1, function(x){

             dfTemp <- data.frame();

             dfTemp <- (igraph::shortest_paths(object@graph , x['mg'],
             as.vector(unlist(completeMetaboliteDF[,1]))))

        return <- dfTemp;
      })
       #' choosing smallest distance between the two metabolites of gene and
       #' all metabolites
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

    df<-  removeDuplicatedColumnDF(df)

     return <- df;
})

setGeneric("associatedShortestPaths", function(object, data)
{
    standardGeneric("associatedShortestPaths")
}
)

#' function that uses the object Graph and calculs distances for every pair of
#' gene - metabolite from data.
#'
#' Note that for 1 gene is attached 1 enzyme which is an edge in our graph.
#' In order to calcul a distance from vertice to vertice this function is
#' called twice in function 'getFinalDFSHortestDistance'.
#'
#' First call is to calcul all the distance for the one vertice of the gene
#' to the all metabolites, the second call is to calcul the distance
#' from the other vertice to all metabolites and then
#' choose the smallest distance between the 2.
#'
#' The selection of the first vertice and the second vertice is done in a other
#' function 'getIdGeneInGraph'
#'
#' param data where
#'  mg = id (from igraph in object Graph) of 1 metabolite related to the gene
#'  m = id (igraph in object Graph)
#'  genes = hsa:... id from KEGG
#'  metabolites = C.... id from KEGG
#'
#' @importFrom igraph shortest_paths
#' @param object Graph, data(mg, m, genes, metabolites)
#' @keywords  kegg
#' @examples associatedShortestPaths(Graph, data)

setMethod("associatedShortestPaths","Graph", function(object, data){

    #'calcul al distances
    pl <-  apply(data,1, function(x){

        dfTemp <- data.frame();
        if(is.na(x['mg']) || is.na(x['m'])){

            dfTemp <- NA;
        }else{
        dfTemp <- (igraph::shortest_paths(object@graph , x['mg'],x['m']))
        }


        return <- dfTemp;
    })

    #' choosing smallest distance between the two metabolites of gene and
    #' all metabolites
    outputFinal <- data.frame();
    output <- NULL;
    pl1 <- lapply(pl, function(x){

       if(!is.na(x)){
       lengthPath <- length(x$vpath[[1]])

        output <- lengthPath;
            #' if length of path is 0 -> couldnt reach a path
            #' if length of path is >0 i have to do length -1
            #' for example path i am lookinf for path from 1523 to 1523
            #' path is 1523 and the length is 1. But the real length is 0.
            if(lengthPath  != 0){
                lengthPath <-  (lengthPath -1);
            }else if(lengthPath  == 0)
                lengthPath  <-  NA;

        output <- lengthPath;
        }else output <- NA
       return <- output;

    })

    outputFinal <- rbind(outputFinal, pl1)
    outputFinal <- t(outputFinal)
    colnames(outputFinal) <- c("lengthShortestPath")

    return <- outputFinal;
})



#' Fonction that calculates distance between each gene-metabolite pairs.
#'
#' The igrpah created to simulate the KEGG pathways as metabolties as nodes
#' and genes (related gene enzymes and reaction) as egdes.
#'
#' The shortest distance is taking from calculation from both vertices related
#' to a gene to the metabolites of interest.
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#'
#' @param data(gene, metabolites )
#' @keywords KEGG
#' @export
#' @examples getDistanceAsso(pathwayId, data)
getDistanceAsso <- function(pathwayId, data, ordered = FALSE){

    #if the xml file was already dowmloaded
    if(isFileInDirectory(pathwayId) == FALSE){
        getPathwayKGML(pathwayId);
    }

    finalDF <- data.frame();

    #graph creation
    graphe <-  createGraphFromPathway(pathwayId);

    #modify function calculate distance directly for association
    finalDF <- getFinalAssoDfSd(graphe, data);

    # Change Na in finalDF to Inf value
    finalDF <- changeDFassosToRigthDistances(finalDF);

    # order result by increasing distances
    finalDF[is.na(finalDF)] <- NaN;

#temporary
   # finalDF <- finalDF[!duplicated(finalDF),]
  #  finalDF <- removeDuplicatedColumnDF(finalDF)
    if(ordered == TRUE){
    finalDF <- finalDF[ order(finalDF[,5]), ]
    }
  #  finalDF <- removeRowsWithConditions(finalDF)
  #  print(finalDF)
    #finalDF <- finalDF[,order('lengthShortestPath')]

    return <- finalDF;

}

changeDFassosToRigthDistances <- function(associatedShortestPathsDF){

    for(row in 1:nrow(associatedShortestPathsDF)){

        # test if both gene and metabolite are in map
        # but no distance is found
        if(associatedShortestPathsDF[row,2] == TRUE &&
           associatedShortestPathsDF[row,4] == TRUE &&
           is.na(associatedShortestPathsDF[row,5])){

             associatedShortestPathsDF[row,5] <- Inf;
        }
    }
 return <- associatedShortestPathsDF;

}

createGraphFromPathway <- function(pathwayId){
    #print("CreateGraphFromPathway")
    #' create df for vertices
    nodeDF <- getListNodeFromKGML(pathwayId);

    #' create df edges
    edgeDF <- finalReactionEdgeDF(pathwayId);

    #' create graphEl objects
    graphEl <- new("GraphElements", nodeDF= nodeDF,
                   edgeDF= edgeDF, pathwayId = pathwayId);


    #' create igraph with graphEl object elements
    igraphe <- createIGraph(graphEl);

    #' create Graph object
    graphe <- new("Graph", graph = igraphe, graphEl);

      return <- graphe;
}


#
#
#     f <- data.frame(f)
#     print(length(f))
#     f[,2] <- unlist(completeMetaboDF[1]);
#     f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
#
#     colnames(f) <- c("m","metabolites") #set colnames
#
#     return <- f;
#
# })


setGeneric("getIdGeneInGraph", function(object, associatedGeneMetaDF,
                                                         indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdGeneInGraph", "Graph", function(object,
                                        associatedGeneMetaDF, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################
   # print("getIdGeneInGraph")
    f <- apply(associatedGeneMetaDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

     })
    f <- t(f)# handling df
    f <- data.frame(f)

    f[,3] <- unlist(associatedGeneMetaDF[1]);

    f <- f[, c(indexMetabolite,3)]

    f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
    colnames(f) <- c("mg","genes") #set colnames
    return <- f;

})

setGeneric("fromAssosDFEntryToIGraphIdDF", function(object, data,
                                               indexMetabolite) {
    standardGeneric("fromAssosDFEntryToIGraphIdDF");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples fromAssosDFEntryToIGraphIdDF(g, data, indexMetabolite)

setMethod("fromAssosDFEntryToIGraphIdDF", "Graph", function(object, data,
                                                       indexMetabolite){
    #########################################################################
    #####  Add condition to insure indexMetabolite can only be 1 or 2   #####
    #########################################################################

    mgmDF <- data.frame();
    f <- apply(data,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

        # ' get metabolites id from Graph of data
        m2 <- getCompoundNodeKgmlId(object@graph, x[2], object@nodeDF);

        temp <- x;

        if(length(m2) > 1){

            f1 <- lapply(m2, function(x){

                mgmDF<- c(mg = m1[indexMetabolite], m = x, temp[1], temp[2]);
                return<- mgmDF;

         })
        }else {

            f1<- c(mg = m1[indexMetabolite], m = m2,  temp[1], temp[2]);

        }

    if(typeof(f1) == 'list'){
        f1 <- do.call(rbind, f1)
    }else{
        f1 <- data.frame(f1);
        f1 <- t(f1)
    }

        return <- f1
    })

    f <- do.call(rbind, f)
    f <- data.frame(f)

 #print(f)
    return <- f;

})




setGeneric("getFinalAssoDfSd", function(object, data)
{
    standardGeneric("getFinalAssoDfSd");
}
)

#' A igraph function
#'
#'
#' Calls getIdGeneInGraph() twice to get :
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' where indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' #' where indexMetabolite = 2
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
#' @examples getCompoundNodeKgmlId(g, data, indexMetabolite)

setMethod("getFinalAssoDfSd", "Graph", function(object, data){

    finalDF <- data.frame();
    tempDF <- data.frame();
    # indexMetabolite = 1 to get the first metabolite attached to gene
    m1g <- fromAssosDFEntryToIGraphIdDF(object, data, 1);
    # indexMetabolite = 2 to get the second metabolite attached to gene
    m2g <- fromAssosDFEntryToIGraphIdDF(object, data, 2);

        #' get all shortest paths for both ends of gene to all metabolites
    r1 <- associatedShortestPaths(object, m1g);
    r2 <- associatedShortestPaths(object, m2g);

    #' Each gene is related to 2 metabolites, choose the shortest distance
    #' between both gene-metabolite to metabolite
    for (row in 1:nrow(r1)) {

        r <- mergeVectorsLowerValues(r1[row,], r2[row,]);
        rf <- t(data.frame(r));
        tempDF <- rbind(tempDF,rf);

        return <- tempDF;
    }

    colnames(tempDF) <- c("distance")
  #create final DF

    m1g <- apply(m1g, 2, function(x) gsub("^$|^ $", NA, x))
    m1g <- m1g[,as.vector(c(3,1,4,2))]
    m1g <- data.frame(m1g)


    m1g[,2] <- !(is.na(m1g[,2]));
    m1g[,4] <- !(is.na(m1g[,4]));
    m1g <- cbind(m1g, as.vector(tempDF));
  #  print(m1g)

   return <- m1g;
})

#'function merge 2 numeric vectors and return a vector with the smalest
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

    # get all shortest oaths from data entry
    shortestsPathsDF <- getDistanceAll(pathwayId, associatedGeneMetaDF,
                                       completeMetaboliteDF  );

    associatedMetabo <- getAssociatedMetaboByGene(associatedGeneMetaDF,gene)
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
#print(shortestsPathsDF.plot)
   # create barplot
    plot <- ggplot2::ggplot(shortestsPathsDF.plot, ggplot2::aes(
        x = Var1,
        y = Freq,
        fill = Associations
      ),environment = environment())

    plot <- (plot + ggplot2::geom_bar(stat="identity")
      + ggplot2::theme_bw()
      + ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"))
      + ggplot2::xlab("Distance from Gene")
      + ggplot2::ylab("Metabolite count")
      + ggplot2::ggtitle(gene)
      + ggplot2::scale_fill_manual(values = c("FALSE" ="grey", "TRUE" = "red"))
      + ggplot2::guides(fill=FALSE)
     );


   # print plot it could change to save the graph image somewhere
   print(plot);

}

#' function that output a heatmap graph showing all results
#'  from getDistanceAll function
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#' @param pathwayId, data(gene, metabolites)
#' @keywords  Graph, heatmap, shortestDistance
#' @export
#' @examples heatmapFunction(hsa01100, data)

heatmapFunction <- function(pathwayId, data){
    data1 <- data.frame(c(data[2]))

    AllSP <- getDistanceAll(pathwayId, data, data1);


    dat <- data.frame( Row = rep(row.names(AllSP), each= length(colnames(AllSP))),
                       Col = rep(colnames(AllSP), times= length(row.names(AllSP))),
                       Y = c(t(AllSP)),
              Associations = rep(c(TRUE,FALSE,FALSE,FALSE,FALSE
                                   ,FALSE,FALSE,FALSE,FALSE,FALSE),times=21)

                       );
#     Associations = rep(c(TRUE,FALSE,FALSE,
#                          FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,
#                          FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,
#                          FALSE,FALSE,FALSE,TRUE),times=10)

    # create frame to color edges of associated genes and metabolites
    frames = dat[dat$Associations, c("Row","Col")]

    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)

    # create heatmap graph
    p2 = ggplot2::ggplot(data=dat) +
        ggplot2::geom_raster(ggplot2::aes(x=Row, y=Col, fill=Y)) +
        ggplot2::theme(
               panel.border = ggplot2::element_blank(),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank()) +
    #      ggplot2::scale_fill_gradient2(low = "#0000FF", mid = "#FFFFFF", high ="#FF0000",
    #                                    na.value="lightgrey", name="")+

        ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
        ggplot2::aes(xmin=Row - 0.5, xmax=Row + 0.5, ymin=Col - 0.5, ymax=Col + 0.5)) +
        ggplot2::geom_text(fill = dat$Y, label = dat$Y, ggplot2::aes(x = Row, y = Col)) +
        ggplot2::labs(title="Heatmap")


        print(p2);

}



getAssociatedMetaboByGene <- function(data, gene){

    selectedRows <- data[data$gene == gene,];
    associatedMetabo <- selectedRows[,"metabolites"];

    return <- associatedMetabo;

}

#' Fonction that calculates every shortest distances between each gene and
#' all metabolites in KEGG pathway of choice
#'
#' The igrpah created to simulate the KEGG pathways as metabolties as nodes
#' and genes (related gene enzymes and reaction) as egdes.
#'
#' The shortest distance is taking from calculation from both vertices related
#' to a gene to the metabolites of interest.
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#'
#' @param data(gene, metabolites )
#' @keywords KEGG
#' @export
#' @examples getDistanceAll(pathwayId, data)
getDistanceAll <- function(pathwayId, associatedGeneMetaDF,
                           completeMetaboliteDF){

 #   print("getDistanceALL")
    finalDF <- data.frame();
    #if the xml file was already dowmloaded
    # look when it was downloaded if it has been to long redownload
    if(isFileInDirectory(pathwayId) == FALSE){
        getPathwayKGML(pathwayId);
    }

    if(!is.data.frame(associatedGeneMetaDF) || length(associatedGeneMetaDF[1,])< 2 ||
       length(associatedGeneMetaDF[1,])> 3){
        e <- simpleError("dataframe dimension is wrong, please enter you data
                         frame with KEGG id of genes (ex : hsa:00001) in first
                         column and associated KEGG id metabolites (ex: C00001)
                         in second column")
        tryCatch(stop(e), finally = print("please try again"))

    }else{
        if(is.null(getKGMLRootNode(pathwayId))){
            print("path you entered do not exist, enter valid hsa number without :");
        }else{


            #graph creation
            graphe <-  createGraphFromPathway(pathwayId);
            finalDF <- getFinalDFSHortestDistance(graphe, associatedGeneMetaDF,
                                                  completeMetaboliteDF );

            # Change Na in finalDF to Inf value
            finalDF[is.na(finalDF)] <- Inf;

            return <- finalDF;
        }
    }

}

setGeneric("getFinalDFSHortestDistance", function(object,  associatedGeneMetaDF,
                                                  completeMetaboliteDF)
{
    standardGeneric("getFinalDFSHortestDistance");
}
)

#' A igraph function
#'
#'
#' Calls getIdGeneInGraph() twice to get :
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' where indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' #' where indexMetabolite = 2
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
#' @examples getCompoundNodeKgmlId(g, data, indexMetabolite)

setMethod("getFinalDFSHortestDistance", "Graph", function(object,
                                 associatedGeneMetaDF, completeMetaboliteDF){
   # print("getFinalDFShortestDistance")
    finalDF <- data.frame();

    # indexMetabolite = 1 to get the first metabolite attached to gene
    idM1g <- getIdGeneInGraph(object,associatedGeneMetaDF, 1);
    # indexMetabolite = 2 to get the second metabolite attached to gene
    idM2g <- getIdGeneInGraph(object,associatedGeneMetaDF, 2);

    idM <- getIdMetabolitesInGraph(object, completeMetaboliteDF)
    idM <- na.omit(idM)



    #' get all shortest paths for both ends of gene to all metabolites
    r1 <- allShortestPaths(object, idM1g  , idM);
    r2 <- allShortestPaths(object, idM2g  , idM);

    #' Each gene is related to 2 metabolites, choose the shortest distance
    #' between both gene-metabolite to metabolite
    for (row in 1:nrow(r1)) {
        r <- mergeVectorsLowerValues(r1[row,], r2[row,]);
        rf <- t(data.frame(r));
        finalDF <- rbind(finalDF,rf);
        return <- finalDF;
    }
   # colnames(finalDF) <-  as.vector(unlist(idM[2]));

    finalDF <- cbind(finalDF,"genes" = as.vector(idM1g[2]) )


    #'  testing the length of finalDF, if bigger then 1 remove duplicate
    #'  if not table of length 1 -> no need to remove duplicates
    if(length(finalDF) > 1){
        finalDF <- finalDF[!duplicated(finalDF), !duplicated(colnames(finalDF))]
    }
    #'identification of rownames with genes from entry data

    row.names(finalDF) <- as.vector(unlist(finalDF['genes']));

    # take last column with gene names
    numColumns <- dim(finalDF)[2] -1;
    finalDF<- subset(finalDF[,1:numColumns]);

    return <- finalDF;
})


setGeneric("getIdGeneInGraph", function(object, associatedGeneMetaDF,
                                        indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdGeneInGraph", "Graph", function(object,
                                                associatedGeneMetaDF, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################
   # print("getIdGeneInGraph")
    f <- apply(associatedGeneMetaDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

    })
    f <- t(f)# handling df
    f <- data.frame(f)

    f[,3] <- unlist(associatedGeneMetaDF[1]);

    f <- f[, c(indexMetabolite,3)]

    f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
    colnames(f) <- c("mg","genes") #set colnames
    return <- f;

})


setGeneric("getIdMetabolitesInGraph", function(object,completeMetaboDF) {
    standardGeneric("getIdMetabolitesInGraph");
}
)

#' A igraph function
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdMetabolitesInGraph", "Graph", function(object,
                                                       completeMetaboDF){

#print("getIdMetabolitesInGraph")
    f <- apply(completeMetaboDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getCompoundNodeKgmlId(object@graph, x[1], object@nodeDF);

        temp <- x;

        if(length(m1) > 1){

            f1 <- lapply(m1, function(x){

                mgmDF<- c(m=x, temp[1]);
                return<- mgmDF;

            })
        }else {

            f1<- c(m=  m1, temp[1]);

        }

        if(typeof(f1) == 'list'){
            f1 <- do.call(rbind, f1)
        }else{
            f1 <- data.frame(f1);
            f1 <- t(f1)
        }

        return <- f1;
    })

    f <- do.call(rbind, f)

    f <- data.frame(f)

    return <- f;
})




#' A igraph function
#'
#' function returning the nodes (sub, prod ) of the egdes representing the gene
#'  of interest
#'
#' @param graphe, gene of interest, dataframe of reaction (edges) in the graph
#' @keywords  graph, edge
#' @examples getHeadTailKgmlIdOfEdge(g, hsa, reactionDF)

getHeadTailKgmlIdOfEdge <- function(g, hsaGene,  reactionDF){

    # return lines of DF that contains hsa gene
    x <-grep(hsaGene, reactionDF$ko)


     # get node object of the same id
    #########################################################################
    #########################################################################
    ###faire ajustement lorsqu'il y a plus d'un noeud ayant ce metabolite####
    #########################################################################
    #########################################################################
    nodesVector1 <- data.frame();
    if(length(x) > 0){
    #print(x)
    nodesVector <- as.vector(igraph::get.edges(g, igraph::E(g)[x[1]]));

    # from data frame id toi kgmlId(used in graph)
    sub <- igraph::V(g)[nodesVector[1]];
    prod <- igraph::V(g)[nodesVector[2]];

    sub <- sub[[1]]$name;
    prod <- prod[[1]]$name;

    nodesVector1 <- c(as.character(sub), as.character(prod));
    }else{
        nodesVector1 <-c("na","na")
    }

    return <- nodesVector1;

}


#' A igraph function
#'
#' function returning the kgmlId of the compound of interest
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  graph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)


getCompoundNodeKgmlId <- function(g, compoundKeggId, nodeDF){
### ajouter le traitement de plusieurs match ... plusieurs fois le compounds
### dans la même map!!
    x <-grep(compoundKeggId, nodeDF$keggId)

    #########################################################################
    #########################################################################
    ###faire ajustement lorsqu'il y a plus d'un noeud ayant ce metabolite####
    #########################################################################
    #########################################################################

    if(length(x) > 0 ){
       # print(x)
        x <- igraph::V(g)[x[1]];
        x <- x[[1]]$name
    } else{
        x <- "na";
    }
    return <- x;

}

#' A igraph function
#'
#' fucntion return distance between 2 compound of the graph, in our case
#' the first (c1) compound is one of the edge of our gene of interest and the second
#' (c2) is the target metabolite
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  graph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)

getShortestDistanceBetween2Compounds <- function(g, c1, c2){


    ## fonctionne seulement lorsqu'il y a une seul distance de calculé
    # retourne la première distance. Si on veut calculer plusieurs distance
    # il faudra ajusté...

    distance <- igraph::shortest.paths(g,
        igraph::V(g)[as.character(c1)], igraph::V(g)[as.character(c2)]);

    return <- distance[1];



}
# ######### mettre comme function de la class GRAPH #################
#
# #' A igraph function
# #'
# #' Modification of dataframe giving in entry by user, from
# #' c("gene", "metabolite") to 2 dataFrame
# #' if indiceMetabolite = 1
# #' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
# #' if indiceMetabolite = 2
# #' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
# #'
# #' @param graphe, kegg id of compound of interest, dataframe of compounds
# #' (nodes)
# #' @keywords  igraph, node, kgmlId
# #' @examples fromDFEntryToIGraphIdDF(g, data, indiceMetabolite)
#
# fromDFEntryToIGraphIdDF<- function(g, data, indiceMetabolite, nodeDF, edgeDF){
#
#     #########################################################################
#     ##### Ajouter une condition sur indiceMetabolite can only be 1 or 2 #####
#     #########################################################################
#
#     mgmDF <- data.frame();
#     f <- apply(data,1, function(x){
#
#
#         m1 <- getHeadTailKgmlIdOfEdge(g, x[1], edgeDF);
#         m2 <- getCompoundNodeKgmlId(g, x[2], nodeDF);
#         mgmDF<- c(m1g = m1[indiceMetabolite], m = m2);
#
#         return<- mgmDF;
#
#     })
#
#
#     return <- f;
#
# }
#
# ######### mettre comme function de la class GRAPH #################
#
#
# #' A igraph function
# #'
# #'
# #' Calls fromDFEntryToIGraphIdDF() twice to get :
# #' Modification of dataframe giving in entry by user, from
# #' c("gene", "metabolite") to 2 dataFrame
# #' where indiceMetabolite = 1
# #' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
# #' #' where indiceMetabolite = 2
# #' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
# #'
# #' This is because genes are edges and we calculate vertex to vertex
# #' the distance between genes and metabolites.
# #' From each gene we want to take the shortest distance between the 2 vertex
# #' attached, so we calcul both and in an other function we take the smallest
# #' and that will be the output.
# #'
# #' This function keeps, between the 2 metabolites of each gene vertexes,
# #' the distance that is the smallest between the 2 to the same metabolite.
# #'
# #' @param graphe, kegg id of compound of interest, dataframe of compounds
# #' (nodes)
# #' @keywords  igraph, node, kgmlId
# #' @examples getCompoundNodeKgmlId(g, data, indiceMetabolite)
#
# getFinalDFSHortestDistance<- function(g, data, nodeDF, edgeDF){
#
#
#     # indiceMetabolite = 1 to get the first metabolite attached to gene
#     m1g <- fromDFEntryTo2DFiGraphID(igraphe, realData, 1, nodeDF, edgeDF)
#     # indiceMetabolite = 2 to get the second metabolite attached to gene
#     m2g <- fromDFEntryTo2DFiGraphID(igraphe, realData, 2, nodeDF, edgeDF)
#
#
#
#
#
#     return <- f;
#
# }


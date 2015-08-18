#' A graph function
#'
#' function returning the nodes (sub, prod ) of the egdes representing the gene of
#' interest
#' @param graphe, gene of interest, dataframe of reaction (edges) in the graph
#' @keywords  graph, edge
#' @export
#' @examples getHeadTailKgmlIdOfEdge(g, hsa, reactionDF)

getHeadTailKgmlIdOfEdge <- function(g, hsaGene,  reactionDF){

    # return lines of DF that contains hsa gene
    x <- pmatch(hsaGene,
                reactionDF$ko, nomatch = NA_integer_, duplicates.ok = FALSE);

    # get node object of the same id
    nodesVector <- as.vector(igraph::get.edges(g, igraph::E(g)[x]));

    # from data frame id toi kgmlId(used in graph)
    sub <- igraph::V(g)[nodesVector[1]];
    prod <- igraph::V(g)[nodesVector[2]];

    sub <- sub[[1]]$name;
    prod <- prod[[1]]$name;

    nodesVector1 <- c(as.character(sub), as.character(prod));

    return <- nodesVector1;

}


#' A graph function
#'
#' function returning the kgmlId of the compound of interest
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  graph, node, kgmlId
#' @export
#' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)


getCompoundNodeKgmlId <- function(g, compoundKeggId, nodeDF){
### ajouter le traitement de plusieurs match ... plusieurs fois le compounds
### dans la même map!!
    x <- pmatch(compoundKeggId,
                nodeDF$keggId, nomatch = NA_integer_, duplicates.ok = FALSE);
    x <- igraph::V(g)[nx];
    x <- x[[1]]$name

    return <- x;

}

#' A graph function
#'
#' fucntion return distance between 2 compound of the graph, in our case
#' the first (c1) compound is one of the edge of our gene of interest and the second
#' (c2) is the target metabolite
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  graph, node, kgmlId
#' @export
#' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)

getShortestDistanceBetween2Compounds <- function(g, c1, c2){


    ## fonctionne seulement lorsqu'il y a une seul distance de calculé
    # retourne la première distance. Si on veut calculer plusieurs distance
    # il faudra ajusté...

    distance <- igraph::shortest.paths(g,
        igraph::V(g)[as.character(c1)], igraph::V(g)[as.character(c2)]);

    return <- distance[1];



}


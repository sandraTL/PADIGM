
#'creation of graphe with iGraph
#'with the list of compounds (nodes) in the KGML
#'and the list of reations related to a gene
#'

createGraph <- function(pathwayId){

    nodeDF <- getListNodeFromKGML(pathwayId);

    reactionDF <- mergeReactionEdgeDF(pathwayId);

    g <- igraph::graph.data.frame(reactionDF, directed=FALSE, vertices=nodeDF);

     return <- g;
}




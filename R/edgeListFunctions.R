

#'Fonction qui sépare les listes de réactions : id R1... R2... R3... ko
#'en plusieurs lignes avec les même id et le même ko.
#' -> id R1 ko
#' -> id R2 ko
#' -> id R3 ko
#'
#' @param data frame des reactions extraites du KGML en listes
#' @keywords  kegg
#' @export
#' @examples completeEdgeList(edgeDataFrame)

unlistEdgeReactionCol <- function(edgeDataFrame){

    s <- strsplit(as.vector(edgeDataFrame$reactions), " ");
    edgeDataFrame <- data.frame(
                       kgmlId = rep(edgeDataFrame$kgmlId, sapply(s, length)),
                       reactions = unlist(s),
                       type = rep(edgeDataFrame$type, sapply(s, length)),
                       ko = rep(edgeDataFrame$ko, sapply(s, length)))

    return <- edgeDataFrame;


}


correctReactionString <- function(edgeDF){

    edgeDF <- lapply(edgeDF,
                           function (x) gsub("rn:","",x))

    return <- edgeDF;

}

#'Since it is not possible to have the related gene when parsing the
#'reaction nodes in the KGML. I had to parse de reaction nodes and the
#'gene entry nodes separatly and then combinde the 2 with the information
#'needed.
#'
#'the results gives a data.frame(substrateId, productId, subtrateName,
#'productName, reactionId, reactionName, ko(gene))
#'
#' @param id of the explored pathway
#' @keywords reaction data.frame
#' @export
#' @examples mergeReactionEdgeDF(pathwayId)

mergeReactionEdgeDF <- function(pathwayId){

    edgeDF <- getListEdgeFromGeneKGML(pathwayId);
    reactionDF <- getListReactionFromKGML(pathwayId);
    reactionDF <- as.data.frame(lapply(reactionDF,
                                       function(X) unname(unlist(X))));

    mergeDF <-  merge(reactionDF, edgeDF, by="reactionId");
    mergeDF <- mergeDF[ -c(8,9)];
    mergeDF <- mergeDF[ c(2,3,4,5,1,6,7,8)];

    return <- mergeDF;

}



### plus tard ortholog
# completeSubProdOrthologReactions <- function(edgeDF){
#
#     #print(edgeDF[2]);
#
#     lapply(edgeDF$reactions,
#            function (x) if(!is.na(x)) getMainReactions(getReactionInfo(x)));
#
#
#
# }


## plus tard ortholog gene
# getFinalEdgeList <- function(pathwayId){
#
#     # get list of edge reaction in KGML type gene and ortholog
#     finalEdgeDF <- rBindOrthologAndGeneDataFrame(pathwayId);
#
#     # separate list of reactions in the same row in multiple rows
#     finalEdgeDF <- unlistEdgeReactionCol(finalEdgeDF);
#
#     # correct string of reaction to utilise KEGGREST search
#     finalEdgeDF <- correctReactionString(finalEdgeDF);
#
#     return <- finalEdgeDF;
#
# }

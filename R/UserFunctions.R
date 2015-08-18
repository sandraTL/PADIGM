
shortestDistanceBtwGeneAndMetabolite <- function(pathwayId, file){

     importFile <- scan("/Users/sandra/Desktop/Metabolomique/
         InbornErrorMetabolism/testGeneMetabolite.txt", skip = 1)
     nodeDF <- getListNodeFromKGML(pathwayId);

     reactionDF <- mergeReactionEdgeDF(pathwayId);

     graphe <- createGraph(pathwayId);


}


importedFileTraitment <- function(file){






}

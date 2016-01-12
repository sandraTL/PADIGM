

#' Permutation test to evauluate if gene-metabolite associations pairs are significantly closer
#' than randomly selected gene-metabolite pairs
#'
#' @param data, data.frame() colunm c(gene, metabolite) representing association between them for each row
#' geneMeasured, list of evaluated genes in study (associated and not associated)
#' metaboliteMeasured, list of measured metabolites in study (associated and not associated)
#' permutation, number of permutation to execute
#' output, "medians" median value of permutated sets - "pvalue" - "histogram"
#'
#' @keywords permutation
#' @export
#' @examples permutationFunction("hsa01100", data, gene, metabolite, 1000)

permutationFunction <- function(pathwayId,data,geneMeasured,metaboliteMeasured,
                     permutation,output=c("medians","pvalue","histogram")){

     # process bar
    # pb <- txtProgressBar(min = 0, max = total, style = 3)
    graphe <-  createGraphFromPathway(pathwayId);
    permutatedMedians<-rep(NA,permutation)



    # change data to eliminate associations not where the gene or the metabolite
    # or both are not in the graph
    geneList<-(data[,1])
    rGeneList<-numberOfReactions(graphe@edgeDF,geneList)
    metaboliteList<-(data[,2])
    rMetaboliteList <- numberOfMetabolites(graphe@nodeDF, metaboliteList)


    tempDf1 <- data.frame(cbind(g1 = as.vector(geneList),
                                g2 = as.vector(as.numeric(rGeneList)),
                                m1 = as.vector(metaboliteList),
                                m2 = as.vector(as.numeric(rMetaboliteList))))

    tempDf1 <- removeNotInGraph(tempDf1)
    rGeneList <- as.numeric(as.vector(tempDf1[,2]))
    rMetaboliteList <- as.numeric(as.vector(tempDf1[,4]))
    print(tempDf1)
    data <- subset(tempDf1, select = c(1,3))

    geneList<-as.vector(unique(data[,1]))
    metaboliteList<-as.vector(unique(data[,2]))


    # making sure there is no doubles in the list of all gene measured
    geneMeasured <- unique(geneMeasured)

    # generate vectors with number of reactions for catalysed by gene
    rPossibleGeneToPermutate<-numberOfReactions(graphe@edgeDF,geneMeasured)

    # making sure there is no doubles in the list of all metabolites measured
    metaboliteMeasured <- as.vector(unique(metaboliteMeasured))

    # generate vectors with number of repetition of each metabolite
    rPossibleMetaboliteToPermutate <- numberOfMetabolites(graphe@nodeDF,
                                                         metaboliteMeasured)




    # remove metabolites not in the graph from list of all metabolites
    tempDf2 <- data.frame(cbind(g1 = as.vector(metaboliteMeasured),
                 g2 = as.vector(as.numeric(rPossibleMetaboliteToPermutate))))

    tempDf2 <- removeNotInGraph(tempDf2)

    metaboliteMeasured <- as.vector(tempDf2[,1])

    rPossibleMetaboliteToPermutate <- as.numeric(as.vector(tempDf2[,2]))



    for(k in 1:permutation)
    {

      #  print(geneList)
      #  print(metaboliteList)

        # remove gene not in the graph from list of all genes
        tempDf <- data.frame(cbind(g1 = as.vector(geneMeasured),
                                   g2 = as.vector(as.numeric(rPossibleGeneToPermutate))))
      #  print(tempDf)
        tempDf <- removeNotInGraph(tempDf)
     #   print(tempDf)

        tempDf <- tempDf[!tempDf$g1 %in% geneList,]
        possibleGeneToPermutate <- as.vector(tempDf[,1])
        rPossibleGeneToPermutate <- as.numeric(as.vector(tempDf[,2]))

        metaboliteShuffled <-sample(metaboliteMeasured,length(metaboliteList))
        print(metaboliteShuffled)
        geneShuffled<-vector()

        # for all associated, replace by permutated genes
        for(i in 1:length(geneList)){
            # pick genes that catalyze (+/-1) the same number of reaction
            if(rGeneList[i]!=1)
            {
                # to get position of possible genes
                possibleGenesReaction<-as.vector(which(abs(rPossibleGeneToPermutate
                                                      - rGeneList[i])<=1))

                # condition if there is no gene that has +/- 1 the number of
                # reaction of the gene to replace then use the same gene..?
                if(length(possibleGenesReaction) == 0){
                   possibleGeneToPermutate <- c(possibleGeneToPermutate,
                                             geneList[i])
                    genePositionShuffled <- length(possibleGeneToPermutate)
                }else{
                    genePositionShuffled <-sample(possibleGenesReaction,1)
                    }
            }
            if(rGeneList[i]==1)
            {
                # to get position of possible genes
                possibleGenesReaction<-which((rPossibleGeneToPermutate
                                              -rGeneList[i])==0)
                # condition if there is no gene that has +/- 1 the number of
                # reaction of the gene to replace then use the same gene..?
                if(length(possibleGenesReaction) == 0){
                    possibleGeneToPermutate <- c(possibleGeneToPermutate,
                                                 geneList[i])
                    genePositionShuffled <- length(possibleGeneToPermutate)
                }else{
                    genePositionShuffled <-sample(possibleGenesReaction,1)
                }
            }
            # add chosen gene to shuffle list
            geneShuffled<-c(geneShuffled,as.character(possibleGeneToPermutate[genePositionShuffled]))

            # take out chosen gene (and his number of reactions) from list of possible genes
            possibleGeneToPermutate<-possibleGeneToPermutate[-genePositionShuffled]
            rPossibleGeneToPermutate<-rPossibleGeneToPermutate[-genePositionShuffled]
        }
      #print(geneShuffled)
      permutatedData <- data.frame(data);

        # create new 'data' where associated genes and metabolites are replace by shuffle genes and metabolites
        # for genes
        for(j in 1:length(geneList)){


              f<-  apply(permutatedData,1, function(x) {

                     permutatedData <- gsub(geneList[j],geneShuffled[j], x)
                     return <- permutatedData;

                                 })
                 permutatedData <- data.frame(t(f));
        }
               for(k in 1:length(metaboliteList)){
                f1<-  apply(permutatedData,1, function(x) {

                     permutatedData <- gsub(metaboliteList[k],metaboliteShuffled[k], x)
                     return <- permutatedData;

                 })
                 permutatedData <- data.frame(t(f1));

         }

      # print(permutatedData)
        # calculate distance with permutated data
        distPermutated<-getDistanceAsso(pathwayId,permutatedData,F)
        print(distPermutated)
        permutatedMedians[k]<-median(distPermutated$distance)
        print(permutatedMedians)
        #process bar
        #setTxtProgressBar(pb, k)
    }
    # get median distance associated
   # distAssociated<-getDistanceAsso(pathwayId,data,F)
   #medianAssociated<-median(distAssociated$distance)

    # output functions
    if(output == "medians"){
        return <- permutatedMedians
    }
    if(output == "pvalue"){
        pvalue<-(sum(permutatedMedians<=medianAssociated)+1)/(permutation+1)
        return <- pvalue
    }
    if(output == "histogram"){
        return <- hist(permutatedMedians, abline(v=medianAssociated,lwd=3, col="red"))
    }
}

numberOfReactions <- function(reactionDF,geneList){

       f <- lapply(geneList, function(x){

           nbreReactions <- length(grep(x, reactionDF$ko))

           return <- nbreReactions;
       })
   f <- do.call(rbind, f)
   f <- as.vector(f)

   return <- f
}

numberOfMetabolites <- function(nodeDF,metaboliteList){

    f <- lapply(metaboliteList, function(x){


        nbreMetabolites <- length(grep(x, nodeDF$keggId))

        return <- nbreMetabolites;
    })
    f <- do.call(rbind, f)
    f <- as.vector(f)

    return <- f
}


removeNotInGraph <- function(df){

    row_sub <- apply(df, 1, function(row) all(row != 0))
    df <- df[row_sub,]

    return <- df;

}






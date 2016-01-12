testIfRowNameExistInDF <- function(dF, rowName){



}

merge2DFWithSmallestValue <- function(df1, df2)

for (row in 1:nrow(df)) {
    r <- mergeVectorsLowerValues(df1[row,], df2[row,]);
    rf <- t(data.frame(r));
    finalDF <- rbind(finalDF,rf);
    return <- finalDF;
}



removeDuplicatedColumnDF <- function(df){

df[is.na(df)] <- Inf;

    df<- t(df)


    df<- aggregate(df,list(row.names(df)),function(x) x[which.min(abs(x))]);

    row.names(df) <- df[,1]
    df <- data.frame(t(df))

    df <- df[-c(1), ]

return <- df;
}


removeRowsDistanceAsso <- function(df){

    # df$geneKEGGId <- factor(df$metabolites, levels=unique(df$metabolites))
    df$geneKEGGId <- I(df$geneKEGGId)
    df$metaboliteKEGGId <- I(df$metaboliteKEGGId)

    aa <- split(df, list(df$metaboliteKEGGId,df$geneKEGGId),drop = TRUE)

    f <- lapply(aa ,function(x){

        tempDF <- data.frame(x)

        tempDF <- tempDF[order(tempDF[,7]),]
        return <-tempDF[1,]
    })

    finalDF <- do.call(rbind.data.frame, f)
    row.names(finalDF) <- row.names(1:length(finalDF[,1]))
    return <- finalDF;

}

removeRowsDistanceAll <- function(df){

    df$geneKEGGId <- I(df$geneKEGGId)
    aa <- split(df, list(df$geneKEGGId),drop = TRUE)
    nbreOfCol <- ncol(df) -1

    f <- lapply(aa ,function(x){

        tempDF <- data.frame(x)
        tempDF <- apply(tempDF,2,sort,decreasing=F)

        if((length(tempDF)/ncol(df)) == 1){
            tempDF <- as.data.frame(tempDF)
            tempDF <- as.data.frame(t(tempDF))
        }else tempDF <- as.data.frame(tempDF)
        return <-tempDF[1,]
    })

    finalDF <- do.call(rbind.data.frame, f)
    return <- finalDF;

}







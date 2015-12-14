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


removeRowsWithConditions <- function(df){

        df %>%
        group_by(df[,1], df[,2]) %>%
        arrange(distance) %>% # in each group, arrange in ascending order by distance
        filter(row_number() == 1)

        return <- df;

}


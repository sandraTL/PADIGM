
#' A KEGG Function
#'
#' This function allows you to get a KGML format of a pathway and save
#' xml format in data directory
#' @param this function parses all data from kegg
#' @keywords  kegg
#' @examples
#' getPathwayKGML("hsa01100")


getPathwayKGML <- function(pathwayId) {

    adressfile <- toStringAdressfile(pathwayId);
    destfile <- toStringDestfile(pathwayId);


    download.file(adressfile, destfile, method = "curl");

}

toStringDestfile <- function(pathwayId){
    #concatenation of pathwayId to set swdir for the xml
    s1 <- "~/";
    setwd(s1);
    s2 <-  toString(pathwayId);
    s3 <- ".xml"
    s4 <- paste(s1,s2, sep= "");
    destfile <- paste(s4, s3, sep="");

    return <- destfile;
}


toStringAdressfile <- function(pathwayId){

    s1 <- "rest.kegg.jp/get/";
    s2 <-  toString(pathwayId);
    s3 <- "/kgml"
    s4 <- paste(s1,s2, sep= "");
    adressfile <- paste(s4, s3, sep="");

    return <- adressfile;
}

isFileInDirectory <- function(pathwayId){
    #concatenation of pathwayId to set swdir for the xml
    bool = FALSE;
    files <- list.files("~/")
    s2 <-  toString(pathwayId);
    s3 <- ".xml"

    file <- paste(s2, s3, sep="");
    m <- match(file, files, nomatch = NA, incomparables = NULL)

    if(is.na(m) == FALSE){
        bool = TRUE;
    }

    return <- bool;
}

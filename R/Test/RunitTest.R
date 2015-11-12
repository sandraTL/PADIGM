#'Fonction qui sépare les listes de réactions : id R1... R2... R3... ko
#'
#' @param data frame des reactions extraites du KGML en listes
#' @keywords  kegg
#' @examples RUnitTest()
RUnitTest <- function(){
    myTestSuite <- RUnit::defineTestSuite("RUnit Example",
                               system.file("R", package = "PADIGM"),
                               testFileRegexp = "runitc2f.R")
    return <- myTestSuite;
}


#'Fonction qui sépare les listes de réactions : id R1... R2... R3... ko
#'
#' @param data frame des reactions extraites du KGML en listes
#' @keywords  kegg
#' @examples c2f(20)
#'
c2f <- function(c) return(9/5 * c + 32)

#'Fonction qui sépare les listes de réactions : id R1... R2... R3... ko
#'
#' @param data frame des reactions extraites du KGML en listes
#' @keywords  kegg
#' @examples test.c2f()
#'
test.c2f <- function() {
    RUnit::checkEquals(c2f(0), 32)
    RUnit::checkEquals(c2f(10), 50)
    RUnit::checkException(c2f("xx"))
}

# test.showNodeDF <- function() {
#     checkEquals(showNodeDF(nodeDF))
#
# }

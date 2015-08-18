#
# setClass("swift",
#             representation = "VIRTUAL"
#             )
# setClass("brob",
#             representation = representation(x="numeric",positive="logical"),
#             prototype = list(x=numeric(),positive=logical()),
#             contains = "swift"
#             )
#
#
# .Brob.valid <- function(object){
#     len <- length(object@positive)
#     if(len != length(object@x)){
#         return("length mismatch")
#     } else {
#         return(TRUE)
#     }
# }
# setValidity("brob", .Brob.valid)




# #
# #
# #
setClass(
    Class= "Grap",
    representation = representation(
         nodeDF = "data.frame"
#         reactionDF = "character",
#         pathwayId = "character"
    )
)
# crGraph <- function(object){
# }

setMethod("length","Grap",function(x){length(x@nodeDF)})

.Grap.crGraph <- function(object) {cat(object@nodeDF)}

#  setMethod("crGraph", "Grap",
#           function(object){  cat(object@nodeDF)  })

setMethod("show", "Grap", function(object){print.b(object)})

#               cat(igraph::graph.data.frame(object@reactionDF,
#                             directed=FALSE, vertices=object@nodeDF));


#  setClass("BMI", representation(weight="numeric", size="numeric"))
# ### a voir quel fonction peuvent être réécrite..... pas n'importe quoi :-/...
#  setMethod("show", "BMI",
#              function(object){cat("BMI=",object@weight/(object@size^2)," \n ")}
#             )


# setClass("Grap", representation(nodeDF="numeric"))
#
# setMethod("slow", "Grap", function(object){cat("Grap=",object@nodeDF," \n ")}
# )

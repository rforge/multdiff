setClass(
         Class = "MultDiffModel",
         representation = representation(
           data = "ContinuousProcess",
           parameters = "list",
           "VIRTUAL"
           )
         )


setMethod(
          "getData",
          "MultDiffModel",
          function(object){
            return(object@data)
          }
          )


setMethod(
          "getParameters",
          "MultDiffModel",
          function(object){
            return(object@parameters)
          }
          )


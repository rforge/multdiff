test.markedPointProcess <- function() {
  thisClass <- "MarkedPointProcess"
  attr(thisClass, "package") <- "processdata"
  checkException(markedPointProcess())
  checkEquals(class(markedPointProcess(data.frame(time = 1:10))), thisClass)
}

test.colNames <- function() {
  checkEquals(colNames(continuousProcess(data.frame())), character())
  CP <- continuousProcess(data.frame(value = 1:10))
  checkEquals(colNames(CP), "value")
  MPP <- markedPointProcess(data.frame(time = c(3,5,8.3)), CP)
  checkEquals(colNames(MPP), c("value", "point"))
  MPP <- markedPointProcess(data.frame(time = c(3,5,8.3),
                                       markType = c("A", "A", "B")),
                            CP)
  checkEquals(colNames(MPP), c("value", "A", "B")) 
}


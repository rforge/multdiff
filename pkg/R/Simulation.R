OUEulerScheme <- function (A, B, D, x0, sqdelta, n) {
   ## TODO: checks!
.Call("cOUEuler", A, B, D, x0, sqdelta, n, PACKAGE = "smd")
}

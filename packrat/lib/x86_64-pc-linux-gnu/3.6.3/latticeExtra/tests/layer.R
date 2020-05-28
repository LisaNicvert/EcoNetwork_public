## Tests for non-standard evaluation in layer()

library(lattice)
library(latticeExtra)
oopt <- lattice.options(panel.error = NULL)

xx <- 1970
xyplot(nhtemp) + layer(panel.abline(v = xx))

layer.test1 <- function() {
    xyplot(nhtemp) + layer(panel.abline(v = xx))
}
layer.test1()

layer.test2 <- function() {
    vline <- 1920
    xyplot(nhtemp) + layer(panel.abline(v = vline))
}
try(plot(layer.test2()), silent = TRUE) # fails

layer.test3 <- function() {
    zz <- 1930
    ls.str(environment())
    xyplot(nhtemp) + layer(panel.abline(v = vline),
                           data = list(vline = zz))
}
layer.test3()

layer.test4 <- function() {
    zz <- 1940
    xyplot(nhtemp) + layer(panel.abline(v = c(zz, xx)),
                           data = environment())
}
layer.test4()

do.call(lattice.options, oopt)




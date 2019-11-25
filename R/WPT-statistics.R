
wptstat <- function(XYP, XU, YU, nperm) {
    N3 <- dim(XYP)[1]
    N1 <- dim(XU)[1]
    N2 <- dim(YU)[1]
    N <- N1 + N2 + N3
    XP <- XYP[, 1]
    YP <- XYP[, 2]
    ## evaluate the weighted test
    variancemlp <- var(XP - YP)
    if (variancemlp == 0) {
        variancemlp <- 0.01
    }
    a <- 2 * N3 / (N + N3)

    Tmlup <- (mean(XU) - mean(YU)) / (sqrt((var(XU) / N1) + (var(YU) / N2)))
    Tmlp <- (mean(XP - YP)) / (sqrt(variancemlp / N3))
    t.ml <- sqrt(a) * Tmlp + (sqrt(1 - a)) * Tmlup


    #-----------------------Permutation----------------------------#

    P <- apply(matrix(rep(1:(2 * N3), nperm), ncol = nperm), 2, sample)
    Px <- matrix(c(XP, YP)[P], ncol = nperm)


    XPstar <- Px[1:N3, ]
    YPstar <- Px[(N3 + 1):(2 * N3), ]



    Pu <- apply(matrix(rep(1:(N1 + N2), nperm), ncol = nperm), 2, sample)
    Pxu <- matrix(c(XU, YU)[Pu], ncol = nperm)
    XUstar <- Pxu[1:N1, ]
    YUstar <- Pxu[(N1 + 1):(N1 + N2), ]



    variancemlpstar <- resample::colVars(XPstar - YPstar)
    variancemlpstar0 <- (variancemlpstar == 0)
    variancemlpstar[variancemlpstar0] <- 0.01

    Tmlpstar <- (colMeans(XPstar - YPstar)) / (sqrt(variancemlpstar / N3))

    Tmlupstar <- (colMeans(XUstar) - colMeans(YUstar)) / (sqrt((resample::colVars(XUstar) / N1) +
                 (resample::colVars(YUstar) / N2)))

    t.mlstar <- sqrt(a) * Tmlpstar + (sqrt(1 - a)) * Tmlupstar


    p1ml <- mean(t.mlstar >= c(t.ml))


    list(t.ml = t.ml, p1ml = p1ml)

}

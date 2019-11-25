MCThmtest <- function(XP, YP, XU, YU, nperm, alpha) {
    nc <- dim(XP)[1]
    nu1 <- dim(XU)[1]
    nu2 <- dim(YU)[1]
    VXYU <- c(XU, YU)
    VXYP <- c(XP, YP)
    #### Paired t-test####
    vart <- var(XP - YP)
    if (vart == 0) {
        vart <- 1 / nc
    }
    Tmlp <- (mean(XP - YP)) / sqrt(vart / nc)

    ## permutation Paired t-test##

    P <- apply(matrix(rep(1:(2 * nc), nperm), ncol = nperm), 2, sample)
    xyps <- VXYP

    nc1 <- nc + 1
    nc2 <- 2 * nc
    xypst <- matrix(xyps[P], nrow = nc2, ncol = nperm)
    xypst1 <- xypst[1:nc, ]
    xypst2 <- xypst[nc1:nc2, ]
    vartst <- resample::colVars(xypst1 - xypst2)
    vartst0 <- (vartst == 0)
    vartst[vartst0] <- 1 / nc
    Tmlpst <- (colMeans(xypst1 - xypst2)) / sqrt(vartst / nc)
    p1pst <- mean(Tmlpst >= c(Tmlp))

    ### Welch Test ###

    itm <- (var(XU) / nu1) + (var(YU) / nu2)
    Tmlup <- (mean(XU) - mean(YU)) / sqrt(itm)

    ## permutation Welch test###

    U <- apply(matrix(rep(1:(nu1 + nu2), nperm), ncol = nperm), 2, sample)

    xyus <- VXYU
    xyust <- matrix(xyus[U], nrow = nu1 + nu2, ncol = nperm)

    xyust1 <- xyust[1:nu1, ]
    xyust2 <- xyust[(nu1 + 1):(nu1 + nu2), ]


    itmst <- (resample::colVars(xyust1) / nu1) +
             (resample::colVars(xyust2) / nu2)
    Tmlupst <- (colMeans(xyust1) - colMeans(xyust2)) / sqrt(itmst)

    p1ust <- mean(Tmlupst >= c(Tmlup))



    list(Tmlp = Tmlp, Tmlup = Tmlup, p1pst = p1pst, p1ust = p1ust)
}

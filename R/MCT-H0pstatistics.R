MCThptest <- function(XP, YP, XU, YU, nperm, alpha) {
    nc <- dim(XP)[1]
    nu1 <- dim(XU)[1]
    nu2 <- dim(YU)[1]
    VXYU <- c(XU, YU)
    VXYP <- c(XP, YP)
    #-------------konitscke pauly test (dependent)-----#

    nc1 <- nc + 1
    nc2 <- 2 * nc

    rvxyp <- rank(VXYP)
    rvxp <- rvxyp[1:nc]
    rvyp <- rvxyp[nc1:nc2]

    rixp <- rank(XP)
    riyp <- rank(YP)
    BM1 <- (1 / nc) * (rvxp - rixp)
    BM2 <- (1 / nc) * (rvyp - riyp)
    BM3 <- BM1 - BM2
    pd <- mean(BM2)
    m <- mean(BM3)

    v <- (sum((BM3 - m)^2)) / (nc - 1)
    v0 <- (v == 0)
    v[v0] <- 1 / nc
    T.kp <- sqrt(nc) * (pd - 1 / 2) / sqrt(v)

    P <- apply(matrix(rep(1:(2 * nc), nperm), ncol = nperm), 2, sample)


    xyps <- VXYP
    rxyps <- rvxyp




    xypstar <- matrix(xyps[P], nrow = nc2, ncol = nperm)
    rxypstar <- matrix(rxyps[P], nrow = nc2, ncol = nperm)
    xypstar1 <- xypstar[1:nc, ]
    xypstar2 <- xypstar[nc1:nc2, ]
    rpstar1 <- rxypstar[1:nc, ]
    rpstar2 <- rxypstar[nc1:nc2, ]
    ripstar1 <- apply(xypstar1, 2, rank)
    ripstar2 <- apply(xypstar2, 2, rank)
    BMstar2 <- (1 / nc) * (rpstar2 - ripstar2)
    BMstar3 <- (1 / nc) * (rpstar1 - ripstar1) - BMstar2
    pdstar <- colMeans(BMstar2)
    mstar3 <- colMeans(BMstar3)
    mstar3 <- matrix(colMeans(BMstar3), ncol = nperm)
    mstar3 <- mstar3[rep(seq_len(nrow(mstar3)), times = nc), ]

    vstar3 <- (colSums((BMstar3 - mstar3)^2)) / (nc - 1)
    vstar30 <- (vstar3 == 0)
    vstar3[vstar30] <- 1 / nc
    T.kpstar <- sqrt(nc) * (pdstar - 0.5) / sqrt(vstar3)
    p1pstar <- mean(T.kpstar >= T.kp)


    #-----estimate  Neubert & Brunner test (Independent samples)----#
    U <- apply(matrix(rep(1:(nu1 + nu2), nperm), ncol = nperm), 2, sample)

    ## the ranks for the whole Independent data together
    rvxyu <- rank(VXYU)

    rvxu <- rvxyu[1:nu1]
    rvyu <- rvxyu[(nu1 + 1):(nu1 + nu2)]

    ## the intern ranks inside each independet sample
    rixu <- rank(XU)
    riyu <- rank(YU)


    Nu <- nu1 + nu2

    mrvxu <- mean(rvxu)
    mrvyu <- mean(rvyu)


    varu1 <- (1 / (nu1 - 1)) * sum((rvxu - rixu - mrvxu + ((nu1 + 1) / 2)) ^ 2)
    varu2 <- (1 / (nu2 - 1)) * sum((rvyu - riyu - mrvyu + ((nu2 + 1) / 2)) ^ 2)

    varu10 <- (varu1 == 0)
    varu20 <- (varu2 == 0)
    varu1[varu10] <- 1 / (4 * nu1)
    varu2[varu20] <- 1 / (4 * nu2)

    Vusq <- Nu * ((1 / nu2) * varu1 + (1 / nu1) * varu2)
    Vu <- sqrt(Vusq)

    T.nb <- ((mean(rvyu) - mean(rvxu)) / Vu) * sqrt(nu1 * nu2 / Nu)


    #---------Permutationstest(Neubert brunner)unpaired----------#

    xyus <- VXYU
    rxyus <- rvxyu


    xyustar <- matrix(xyus[U], nrow = nu1 + nu2, ncol = nperm)
    rxyustar <- matrix(rxyus[U], nrow = nu1 + nu2, ncol = nperm)

    xyustar1 <- xyustar[1:nu1, ]
    xyustar2 <- xyustar[(nu1 + 1):(nu1 + nu2), ]

    ## the ranks for the whole data together
    rustar1 <- rxyustar[1:nu1, ]
    rustar2 <- rxyustar[(nu1 + 1):(nu1 + nu2), ]

    ## the intern ranks inside each sample
    riustar1 <- apply(xyustar1, 2, rank)
    riustar2 <- apply(xyustar2, 2, rank)



    mrustar1 <- matrix(colMeans(rustar1), ncol = nperm)
    mrustar1 <- mrustar1[rep(seq_len(nrow(mrustar1)), times = nu1), ]

    mrustar2 <- matrix(colMeans(rustar2), ncol = nperm)
    mrustar2 <- mrustar2[rep(seq_len(nrow(mrustar2)), times = nu2), ]

    varustar1 <- (1 / (nu1 - 1)) * colSums((rustar1 - riustar1 - mrustar1 + ((nu1 + 1) / 2)) ^ 2)
    varustar2 <- (1 / (nu2 - 1)) * colSums((rustar2 - riustar2 - mrustar2 + ((nu2 + 1) / 2)) ^ 2)
    varustar10 <- (varustar1 == 0)
    varustar20 <- (varustar2 == 0)
    varustar1[varustar10] <- 1 / (4 * nu1)
    varustar2[varustar20] <- 1 / (4 * nu2)

    Vustar <- sqrt(Nu * ((1 / nu2) * varustar1 + (1 / nu1) * varustar2))



    T.nbstar <- ((colMeans(rustar2) - colMeans(rustar1)) / Vustar) * sqrt(nu1 * nu2 / Nu)


    p1ustar <- mean(T.nbstar >= T.nb)


    list(Tmlp = T.kp, Tmlup = T.nb, p1pst = p1pstar, p1ust = p1ustar)
}


PBTstat <- function(XYP, XU, nbsp = 1000) {
    NC <- dim(XYP)[1]
    NU <- dim(XU)[1]
    N <- NC + NU
    XP <- XYP[, 1]
    YP <- XYP[, 2]
    A <- matrix(c(1, -1, 0, 0, -1, 1), ncol = 3, nrow = 2, byrow = TRUE)

    k1 <- NC / N
    k2 <- 1 - k1


    SigmaXSqComp <- var(XP)
    SigmaXSqInComp <- var(XU)
    SigmaYSq <- var(YP)
    CorrXY <- cor(XP, YP)

    if (is.na(var(XU))) {
        SigmaXSqInComp <- 0
        print("varianceXUequalzero")
    }


    CovMat <- matrix(c(SigmaXSqComp * (1 / k1), CorrXY * sqrt(SigmaXSqComp * SigmaYSq) * (1 / k1), 0,
            CorrXY * sqrt(SigmaXSqComp * SigmaYSq) * (1 / k1), SigmaYSq * (1 / k1), 0,
            0, 0, SigmaXSqInComp * (1 / k2)), ncol = 3, nrow = 3, byrow = TRUE)

    CovHatInv <- MASS::ginv(A %*% CovMat %*% t(A))

    fZ <- sqrt(N) * matrix(c(mean(XP) - mean(YP), mean(XU) - mean(YP)))

    t.BL <- t(fZ) %*% CovHatInv %*% fZ

    a.BL <- t(fZ) %*% fZ / (sum(diag(A %*% CovMat %*% t(A))))
    m.BL <- t(fZ) %*% diag(diag(CovHatInv)) %*% fZ




    ############# Parametric Bootstrap ###################

    ##### estimating the variance for Parametric Bootstrap calculations!
    sigmaY <- var(YP)
    sigmaX <- var(c(XP, XU))

    corrVal <- cor(XYP)[1, 2] * sqrt(sigmaX) * sqrt(sigmaY)


    CovData <- matrix(c(sigmaX, corrVal, corrVal, sigmaY), ncol = 2, nrow = 2,
                      byrow = TRUE)


    #### calculating the p values of the parametric bootstrap tests ####

    Par <- sapply(1:nbsp, function(i, ...) {

        XYPboot <- mvtnorm::rmvnorm((NC + NU), mean = c(0, 0), sigma = CovData)
        XPstar <- XYPboot[1:NC, 1]
        YPstar <- XYPboot[1:NC, 2]

        XUstar <- XYPboot[(NC + 1):(NC + NU), 1]

        SigmaXSqCompstar <- var(XPstar)
        SigmaXSqInCompstar <- var(XUstar)
        SigmaYSqstar <- var(YPstar)
        CorrXYstar <- cor(XPstar, YPstar)



        CovMatstar <- matrix(c(SigmaXSqCompstar * (1 / k1), CorrXYstar * sqrt(SigmaXSqCompstar * SigmaYSqstar) * (1 / k1),
            0, CorrXYstar * sqrt(SigmaXSqCompstar * SigmaYSqstar) * (1 / k1), SigmaYSqstar * (1 / k1), 0, 0, 0, SigmaXSqInCompstar *
                (1 / k2)), ncol = 3, nrow = 3, byrow = TRUE)

        CovHatInvstar <- MASS::ginv(A %*% CovMatstar %*% t(A))

        fZstar <- sqrt(N) * matrix(c(mean(XPstar) - mean(YPstar), mean(XUstar) - mean(YPstar)))

        t.BLstar <- t(fZstar) %*% CovHatInvstar %*% fZstar

        a.BLstar <- t(fZstar) %*% fZstar / (sum(diag(A %*% CovMatstar %*% t(A))))
        m.BLstar <- t(fZstar) %*% diag(diag(CovHatInvstar)) %*% fZstar

        return(c(t.BLstar, a.BLstar, m.BLstar))
    })

    ecdf_WTPS <- ecdf(Par[1, ])
    pvalueWTPS <- 1 - ecdf_WTPS(t.BL)

    ecdf_ATPS <- ecdf(Par[2, ])
    pvalueATPS <- 1 - ecdf_ATPS(a.BL)

    ecdf_MTPS <- ecdf(Par[3, ])
    pvalueMTPS <- 1 - ecdf_MTPS(m.BL)

    list(WTS = t.BL, ATS = a.BL, MTS = m.BL, pvWTPS = pvalueWTPS, pvATPS = pvalueATPS, pvMTPS = pvalueMTPS)

}

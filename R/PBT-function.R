#' Parametric bootstrap tests for paired data with missing values
#' in a single arm
#'
#' The PBT function calculates the asymptotic model based bootstrap version
#'  of three different quadratic forms for testing H0m: {mu1 = mu2};
#'  the Wald-type statistic, the ANOVA-type
#'  statistic, and the MATS-type statistic for paired data with
#'  missingness in a single arm.
#'
#' @param x a (non-empty) numeric vector of data values
#' representing the first components of the pairs.
#' @param y a (non-empty) numeric vector of data values
#' representing the second components of the pairs.
#' @param nbsp The number of bootstraps used for calculating the
#' bootstrapped tests. The default option is 1000.
#' @param alpha A number specifying the significance level;
#' the default is 0.05.
#'
#' @return
#' A \code{PBT} object containing the following components:
#' \item{Input}{The data input of the user.}
#' \item{Descriptive}{Some descriptive statistics of the data.
#'   Displayed are the number of individuals per variable,
#'   the mean, and variance.}
#'  \item{PBT}{The value of the Wald-type statistic, the ANOVA-type statistic,
#'   and the MATS-type
#'  staistic and the p-value of the bootstrapped procedure of each test.}
#'
#'
#' @references Amro, L., Pauly, M., and Ramosaj, B. (2019). Asymptotic based
#' bootstrap approach for matched pairs with missingness in a single-arm.
#'
#' @export
PBT <- function(x, y, nbsp = 1000, alpha = 0.05) {

    #------------------------------Basic Checks------------------------------#

    if (alpha >= 1 || alpha <= 0)
        stop("'alpha' must be a single number between 0 and 1!")

    if (is.numeric(x) == FALSE)
        stop("x must be numeric !")

    if (is.vector(x) == FALSE)
        stop("x must be vector !")

    if (is.numeric(y) == FALSE)
        stop("y must be numeric !")

    if (is.vector(y) == FALSE)
        stop("y must be vector !")


    #-----------------------Arrange the data---------------------------------#




    data <- cbind(x, y)
    voll <- (!is.na(data[, 1]) & !is.na(data[, 2]))
    un1 <- (is.na(data[, 1]) & !is.na(data[, 2]))
    un2 <- (!is.na(data[, 1]) & is.na(data[, 2]))
    un3 <- (is.na(data[, 1]) & is.na(data[, 2]))
    nc <- length(voll[voll == 1])
    nu1 <- length(un2[un2 == 1])
    nu2 <- length(un1[un1 == 1])
    nu3 <- length(un3[un3 == 1])
    n <- nc + nu1 + nu2 + nu3
    gruppe <- rep(1, n)
    gruppe[un1] <- 3
    gruppe[un2] <- 2
    gruppe[un3] <- 0



    xp <- data[, 1][gruppe == 1]
    yp <- data[, 2][gruppe == 1]
    xu <- matrix(c(data[, 1][gruppe == 2]))
    yu <- matrix(c(data[, 2][gruppe == 3]))

    if (nrow(xu) != 0) {
        xyp <- matrix(c(xp, yp), ncol = 2)
        xu <- xu
    }
    if (nrow(yu) != 0) {
        xyp <- matrix(c(yp, xp), ncol = 2)
        xu <- yu
    }
    nu <- nrow(xu)


    #---------------Needed Checks for calculating Statistics----------------#
    if (nu1 == 0 & nu2 == 0)
        stop("There are no missing values")
    if (nu1 != 0 & nu2 != 0)
        stop("Missing values exist in both variables")
    if (nu < 2)
        stop("very few missing observations")
    if (nc < 2)
        stop("not enough observations")

    #-------------------Compute the test statistics-------------------#


    PBTresult <- PBTstat(xyp, xu, nbsp)

    PWTPS <- PBTresult$pvWTPS
    PATPS <- PBTresult$pvATPS
    PMTPS <- PBTresult$pvMTPS


    #-----------------------Arrange the output----------------------------#

    cat("Parametric Bootstrap Test\n")


    var_names <- c("x", "y")
    Tp0 <- c(PBTresult$WTS, PBTresult$ATS, PBTresult$MTS)
    PvalTpbt <- c(PWTPS, PATPS, PMTPS)
    PBT_output <- cbind.data.frame(Tp0, PvalTpbt)
    colnames(PBT_output) <- cbind("Statistic", "p.value")
    rownames(PBT_output) <- rbind("Wald", "ANOVA", "MATS")



    var_names <- c("x", "y")
    descriptive <- cbind.data.frame(var_names, c(nc + nu1, nc + nu2),
                   c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE)),
        c(var(x, na.rm = TRUE), var(y, na.rm = TRUE)))
    colnames(descriptive) <- c("Variable", "n", "Means", "Variances")

    output <- list()
    class(output) <- "PBT"
    output$Input <- data
    output$Descriptive <- descriptive
    output$PBT <- PBT_output

    return(output)




}  ##End

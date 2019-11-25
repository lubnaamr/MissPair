#' A weighted permutation test for paired data with missing values in both arms
#'
#' The WPT function calculates the weighted studentized permutation test
#' for testing H0m: {mu1 = mu2} for matched pairs with missingness in both arm.
#' WPT test #' is a combination of the randomization version
#' of the paired t-test (Konietschke and Pauly, 2014)
#' and Janssen's permutation version for the Welch test
#' (Janssen, 1997) for complate and incomplete data respectively.
#'
#' @param x a (non-empty) numeric vector of data values
#' representing the first components of the pairs.
#' @param y a (non-empty) numeric vector of data values
#' representing the second components of the pairs.
#' @param alternative a character string specifying the alternative hypothesis,
#'  must be one of 'two.sided' (default), 'greater' or 'less'.
#'  You can specify just the initial letter.
#' @param nperm The number of permutations used for calculating
#' the weighted permutation test. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#'
#' @return A \code{WPT} object containing the following components:
#' \item{Input}{The data input of the user.}
#' \item{Descriptive}{Some descriptive statistics of the data.
#'   Displayed are the number of individuals per variable,
#'   the mean, and variance.}
#'  \item{WPT}{The value of the weighted test staistic
#'   and the p-value of the permutation procedure, as well as the test
#'   decision of the weighted permutation test.}
#'
#'
#'
#' @references Janssen, A. (1997). Studentized permutation tests
#' for non-iid hypotheses and the generalized Behrens-Fisher problem.
#' Statistics & probability letters, 36(1), 9-21.
#'
#'
#' Konietschke, F., & Pauly, M. (2014). Bootstrapping and permuting paired
#' t-test type statistics. Statistics and Computing, 24(3), 283-296.
#'
#'
#' Amro, L., and Pauly, M. (2017). Permuting incomplete paired
#' data: a novel exact and asymptotic correct randomization test.
#' Journal of Statistical Computation and Simulation, 87(6), 1148-1159.
#'
#' @export
WPT <- function(x,
                y,
                alternative = c("two.sided", "less", "greater"),
                nperm = 10000,
                alpha = 0.05) {



    #-------------------------Basic Checks---------------------------#
    if (alpha >= 1 || alpha <= 0)
        stop("'alpha' must be a single number between 0 and 1!")

    alternative <- match.arg(alternative)

    if (is.numeric(x) == FALSE)
        stop("x must be numeric !")

    if (is.vector(x) == FALSE)
        stop("x must be vector !")

    if (is.numeric(y) == FALSE)
        stop("y must be numeric !")

    if (is.vector(y) == FALSE)
        stop("y must be vector !")


    #-----------------------Arrange the data---------------------------#

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

    xyp <- matrix(c(xp, yp), ncol = 2)


    #-------------Needed Checks for calculating Statistics-----------------#
    if (nu1 == 0 & nu2 == 0)
        stop("There are no missing values")
    if (nu1 == 0 || nu2 == 0)
        stop("Missing values exist only in one variable")
    if (nu1 < 2 || nu2 < 2)
        stop("very few missing observations")
    if (nc < 2)
        stop("not enough observations")

    #-----------------Compute the test statistics------------------#
    Tp0 <- wptstat(xyp, xu, yu, nperm)$t.ml

    P1Tml <- wptstat(xyp, xu, yu, nperm)$p1ml


    switch(alternative, two.sided = {
        AltHypothesis <- paste("alternative hypothesis: true difference in means is not equal to 0")
        PvalTml <- min(2 - 2 * P1Tml, 2 * P1Tml)

    }, less = {
        AltHypothesis <- paste("alternative hypothesis: true difference in means is less than 0")
        PvalTml <- P1Tml

    }, greater = {
        AltHypothesis <- paste("alternative hypothesis: true difference in means is greater than 0")
        PvalTml <- 1 - P1Tml

    })



    if (PvalTml <= alpha) {
        TestDec <- noquote(paste0("WPT rejects the null hypothesis at the ", alpha * 100, "% significance level"))

    } else {
        TestDec <- noquote(paste0("WPT fails to reject the null hypothesis at the ", alpha * 100, "% significance level"))
    }

    #-----------------------Arrange the output----------------------------#

    cat("Weighted Permutation Test\n")
    cat(AltHypothesis, "\n")

    var_names <- c("x", "y")
    WPT_output <- cbind.data.frame(Tp0, PvalTml)
    colnames(WPT_output) <- cbind("Statistic", "p.value")


    descriptive <- cbind.data.frame(var_names, c(nc + nu1, nc + nu2), c(mean(c(xp, xu)), mean(c(yp, yu))), c(var(c(xp,
        xu)), var(c(yp, yu))))
    colnames(descriptive) <- c("Variable", "n", "Means", "Variances")

    output <- list()
    class(output) <- "WPT"
    output$Input <- data
    output$Descriptive <- descriptive
    output$WPT <- WPT_output
    output$TestDecision <- TestDec



    return(output)
}

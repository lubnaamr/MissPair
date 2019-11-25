#' Multiplication-combination tests for paired data with missing values in
#' both arms
#'
#' The MCT function calculates the multiplication combination tests
#' for testig H0m or H0p hypothesis for matched pairs with missingness
#' in both arms. For testing the Behrens-Fisher problem H0p: {p = 1/2},
#' MCT is based upon permutation versions of the Munzel (1999) (paired)
#' and Brunner and Munzel (2002) (unpaired) tests. And, For testing
#' the null hypothesis H0m: {mu1 = mu2}, MCT is based upon the
#' randomization version of the paired t-test (Konietschke
#' and Pauly, 2014) and Janssen's permutation version
#' for the Welch test (Janssen, 1997) for paired and unpaired
#' data respectively.
#'
#' @param x a (non-empty) numeric vector of data values
#'  representing the first components of the pairs.
#' @param y a (non-empty) numeric vector of data values
#' representing the second components of the pairs.
#' @param hypothesis a character string indicating which hypothesis is
#' to be tested, must be either 'h0m' (default) or 'h0p'.
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of 'two.sided' (default), 'greater' or 'less'.
#'   You can specify just the initial letter.
#' @param nperm The number of permutations used for calculating the permutation
#'   tests. The default option is 10,000.
#' @param alpha A number specifying the significance level;
#' the default is 0.05.
#'
#' @return
#' A \code{MCT} object containing the following components:
#' \item{Input}{The data input of the user.}
#' \item{Descriptive}{Some descriptive statistics of the data.
#'   Displayed are the number of individuals per variable,
#'   the mean, and variance.}
#'  \item{MCT}{The value of the test staistics of the completely observed
#'  and the incompletely observed data and the p-value of the permutation
#'  procedure of each test, as well as the overall test decision
#'  of the multiplication combination test.}
#'
#'
#' @references Janssen, A. (1997). Studentized permutation tests
#' for non-iid hypotheses and the generalized Behrens-Fisher problem.
#' Statistics & probability letters, 36(1), 9-21.
#'
#'
#' Munzel, U.(1999). Nonparametric methods for paired samples.
#' Statistica Neerlandica, 53(3), 277-286.
#'
#'
#'
#' Brunner, E., Munzel, U. (2000). The Nonparametric
#' Behrens-Fisher Problem: Asymptotic Theory and a Small Sample Approximation.
#' Biometrical Journal 42, 17 -25.
#'
#'
#' Konietschke, F., & Pauly, M. (2014). Bootstrapping and permuting paired
#' t-test type statistics. Statistics and Computing, 24(3), 283-296.
#'
#'
#' Amro, L., Konietschke, F., and Pauly, M.(2019).
#' Multiplication-combination tests for incomplete paired data.
#' Statistics in Medicince, 38(17), 3243-3255.
#'
#' @export
MCT <- function(x, y, hypothesis = c("h0m", "h0p"),
                alternative = c("two.sided", "less", "greater"),
                nperm = 10000,
                alpha = 0.05) {

    #------------------------------Basic Checks------------------------------#
    if (alpha >= 1 || alpha <= 0)
        stop("'alpha' must be a single number between 0 and 1!")

    alternative <- match.arg(alternative)
    hypothesis <- match.arg(hypothesis)

    if (is.numeric(x) == FALSE)
        stop("x must be numeric !")

    if (is.vector(x) == FALSE)
        stop("x must be vector !")

    if (is.numeric(y) == FALSE)
        stop("y must be numeric !")

    if (is.vector(y) == FALSE)
        stop("y must be vector !")


    #-----------------------Arrange the data--------------------------------#


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


    xp <- matrix(data[, 1][gruppe == 1], ncol = 1)
    yp <- matrix(data[, 2][gruppe == 1], ncol = 1)
    xu <- matrix(c(data[, 1][gruppe == 2]))
    yu <- matrix(c(data[, 2][gruppe == 3]))


    #----------------Needed Checks for calculating Statistics----------------#
    if (nu1 == 0 & nu2 == 0)
        stop("There are no missing values")
    if (nu1 == 0 || nu2 == 0)
        stop("Missing values exist only in one variable")
    if (nu1 < 2 || nu2 < 2)
        stop("very few missing observations")
    if (nc < 2)
        stop("not enough observations")

    #---------------------Compute the test statistics----------------------#


    switch(hypothesis, h0m = {
        TestHypothesis <- paste("Testing h0m hypothesis:")
        MCTresult <- MCThmtest(xp, yp, xu, yu, nperm, alpha)

    }, h0p = {
        TestHypothesis <- paste("Testing h0p hypothesis:")
        MCTresult <- MCThptest(xp, yp, xu, yu, nperm, alpha)

    })

    P1Tpd <- MCTresult$p1pst
    P1Tupd <- MCTresult$p1ust

    switch(alternative, two.sided = {
        AltHypothesis <- paste(if (hypothesis == "h0m") "alternative hypothesis: true difference in means is not equal to 0" else "alternative hypothesis: true relative effect p is not equal to 1/2")
        PvalTpd <- min(2 - 2 * P1Tpd, 2 * P1Tpd)
        PvalTupd <- min(2 - 2 * P1Tupd, 2 * P1Tupd)

    }, less = {
        AltHypothesis <- paste(if (hypothesis == "h0m") "alternative hypothesis: true difference in means is less than 0" else "alternative hypothesis: true relative effect p is less than 1/2")
        PvalTpd <- P1Tpd
        PvalTupd <- P1Tupd

    }, greater = {
        AltHypothesis <- paste(if (hypothesis == "h0m") "alternative hypothesis: true difference in means is greater than 0" else "alternative hypothesis: true relative effect p is greater than 1/2")
        PvalTpd <- 1 - P1Tpd
        PvalTupd <- 1 - P1Tupd

    })



    Decision <- c()
    if (PvalTpd <= sqrt(alpha)) {
        Decision[1] <- noquote("p.value <= sqrt(alpha)")
    } else {
        Decision[1] <- noquote("p.value > sqrt(alpha)")
    }

    if (PvalTupd <= sqrt(alpha)) {
        Decision[2] <- noquote("p.value <= sqrt(alpha)")
    } else {
        Decision[2] <- noquote("p.value > sqrt(alpha)")
    }



    #-----------------------Arrange the output----------------------------#

    cat("Multiplication Combination Test\n")
    cat(TestHypothesis, "\n")
    cat(AltHypothesis, "\n")

    var_names <- c("x", "y")
    Tp0 <- c(MCTresult$Tmlp, MCTresult$Tmlup)
    PvalTmct <- c(PvalTpd, PvalTupd)
    MCT_output <- cbind.data.frame(Tp0, PvalTmct, Decision)
    colnames(MCT_output) <- cbind("Statistic", "p.value", "Result")
    rownames(MCT_output) <- rbind("paired data", "unpaired data")

    if ((PvalTpd <= sqrt(alpha)) * (PvalTupd <= sqrt(alpha))) {
        TestDec <- noquote(paste0("MCT rejects the null hypothesis at the ", alpha * 100, "% significance level"))
    } else {
        TestDec <- noquote(paste0("MCT fails to reject the null hypothesis at the ", alpha * 100, "% significance level"))
    }

    if (hypothesis == "h0p") {

        descriptive <- cbind.data.frame(var_names, c(nc + nu1, nc + nu2))
        colnames(descriptive) <- c("Variable", "n")
    } else {

        descriptive <- cbind.data.frame(var_names, c(nc + nu1, nc + nu2), c(mean(c(xp, xu)), mean(c(yp, yu))), c(var(c(xp,
            xu)), var(c(yp, yu))))
        colnames(descriptive) <- c("Variable", "n", "Means", "Variances")
    }




    output <- list()
    class(output) <- "MCT"
    output$Input <- data
    output$Descriptive <- descriptive
    output$MCT <- MCT_output
    output$TestDecision <- TestDec



    return(output)

}


## G matrices
create_G <- function(p, d) {

    u01 <- with(p, sigma_subject_intercept * sigma_subject_slope * cor_subject)
    v01 <- with(p, sigma_cluster_intercept * sigma_cluster_slope * cor_cluster)

    pp <- prepare_paras(p)

    n1 <- p$n1
    tot_n2 <- get_tot_n(p)
    n2_cc <- unlist(pp$control$n2)
    n2_tx <- unlist(pp$treatment$n2)

    n3 <- get_n3(p)
    n3_cc <- n3$control
    n3_tx <- n3$treatment
    if(p$partially_nested) {
        n2_cc <- tot_n2$control
    } else {
        if(length(n2_cc) == 1) n2_cc <- rep(n2_cc, n3_cc)
    }
    if(length(n2_tx) == 1) {
        n2_tx <- rep(n2_tx, n3_tx)
    }


    time <- get_time_vector(p)

    A <- Matrix::Diagonal(tot_n2$total)
    B <- array(c(rep(1, n1), time), dim=c(n1, 2))
    X <- kronecker(A, B)

    ## Misinng data
    X <- X[!is.na(d$y), ]

    G <- list()
    ## G 1 intercept
    if(p$sigma_subject_intercept != 0) {
        X1 <- X[, c(1, 1 + (1:(tot_n2$total-1)) * 2)]
        G1 <- tcrossprod(X1)

        G <- c(G, G1)
    }

    Z <- matrix(c(0,1,1,0), ncol = 2)
    ## G 2 cov
    if(u01 != 0) {
        G2 <- X %*% kronecker(Diagonal(tot_n2$total), Z) %*% t(X)
        G <- c(G, G2)
    }

    ## G 3 slope
    if(p$sigma_subject_slope != 0) {
        X2 <- X[, c((1:(tot_n2$total))*2)]
        G3 <- tcrossprod(X2)
        G <- c(G, G3)
    }

    Z2 <- bdiag(lapply(c(n2_cc, n2_tx), create_Z_block))
    if(p$partially_nested) {
        n3_cc <- 1
        ind <- c(2, (1:n3_tx)*2 + 2)
        Z2[1:(n2_cc*2), 1:2] <- 0
    } else {
        ind_cc <- (1:n3_cc)*2
        ind_tx <- (1:n3_tx)*2 + max(ind_cc)
        ind <- c(ind_cc, ind_tx)
    }

    ## G4 intercept cluster
    if(p$sigma_cluster_intercept != 0) {
        X3 <- X %*% Z2
        X3i <- X3[, ind - 1]
        G4 <- tcrossprod(X3i)

        G <- c(G, G4)
    }
    X3 <- X %*% Z2

    ## G5 cov cluster
    if(v01 != 0) {
        G5 <- X3 %*% kronecker(Diagonal(n3_cc + n3_tx), Z) %*% t(X3)
        G <- c(G, G5)
    }

    ## G6 slope cluster
    if(p$sigma_cluster_slope != 0) {
        X32 <- X3[, ind]
        G6 <- tcrossprod(X32)

        G <- c(G, G6)
    }

    ## G7 sigma
    tot_n <- nrow(X)
    Ge <- Diagonal(tot_n)
    G <- c(G, Ge)

    G
}

# Approximate asymptotic covariance of random effects
# from package pbkrtest
vcovAdj16_internal <- function (Phi, SigmaG, X)
{
    SigmaInv <- SigmaG$iV
    n.ggamma <- SigmaG$n.ggamma
    TT <- SigmaInv %*% X
    HH <- OO <- vector("list", n.ggamma)
    for (ii in 1:n.ggamma) {
        HH[[ii]] <- SigmaG$G[[ii]] %*% SigmaInv
        OO[[ii]] <- HH[[ii]] %*% X
    }
    PP <- QQ <- NULL
    for (rr in 1:n.ggamma) {
        OrTrans <- t(OO[[rr]])
        PP <- c(PP, list(Matrix::forceSymmetric(-1 * OrTrans %*% TT)))
        for (ss in rr:n.ggamma) {
            QQ <- c(QQ, list(OrTrans %*% SigmaInv %*% OO[[ss]]))
        }
    }
    Ktrace <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)
    for (rr in 1:n.ggamma) {
        HrTrans <- t(HH[[rr]])
        for (ss in rr:n.ggamma) {
            Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(HrTrans *
                                                        HH[[ss]])
        }
    }
    IE2 <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)
    for (ii in 1:n.ggamma) {
        Phi.P.ii <- Phi %*% PP[[ii]]
        for (jj in c(ii:n.ggamma)) {
            www <- .indexSymmat2vec(ii, jj, n.ggamma)
            IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 *
                sum(Phi * QQ[[www]]) + sum(Phi.P.ii * (PP[[jj]] %*%
                                                           Phi))
        }
    }
    eigenIE2 <- eigen(IE2, only.values = TRUE)$values
    condi <- min(abs(eigenIE2))
    WW <- if (condi > 1e-10)
        Matrix::forceSymmetric(2 * solve(IE2))
    else Matrix::forceSymmetric(2 * MASS::ginv(IE2))

    attr(WW, "P") <- PP
    WW
}
# from package pbkrtest
.indexSymmat2vec <- function (i, j, N)
{
    k <- if (i <= j) {
        (i - 1) * (N - i/2) + j
    }
    else {
        (j - 1) * (N - j/2) + i
    }
}

get_balanced_df <- function(object) {
    pp <- prepare_paras(object)

    tot_n2 <- get_tot_n(object)
    n2_cc <- tot_n2$control
    n2_tx <- tot_n2$treatment

    n3 <- get_n3(object)
    n3_cc <- n3$control
    n3_tx <- n3$treatment


    if(object$sigma_cluster_slope == 0) {
        df <-  (n2_tx + n2_cc) - 2
    } else if(object$partially_nested) {
        df <- n3_tx - 1
    } else {
        df <- (n3_cc + n3_tx) - 2

    }
    df
}

get_satterth_df <- function(object, d, pars, Lambdat, X, Zt, L0, Phi, varb) {
    A <- Lambdat %*% Zt
    L <- as(L0, "sparseMatrix")
    pvec <- L0@perm + 1L
    P <- as(pvec, "pMatrix")
    I <- Matrix::Diagonal(ncol(A))
    iL <- solve(L)
    PA <- P %*% A
    iLPA <- iL %*% PA
    sigma2 <- pars["sigma"]
    iV <- (I - crossprod(iLPA))/sigma2
    V <- sigma2 * (crossprod(A) + I)

    SigmaG <- list(G = create_G(object, d = d))
    SigmaG$Sigma <- V
    SigmaG$iV <- iV
    SigmaG$n.ggamma <- length(SigmaG$G)

    ## delta method
    vv <- vcovAdj16_internal(Phi, SigmaG, X)
    Lc <- c(0,0,0,1)
    g <- gradient(function(x)  as.numeric(varb(x = x, Lc)), x = pars[pars != 0], delta = 1e-4)
    df <- 2*(Phi[4,4])^2 / (t(g) %*% vv %*% g)
    df <- as.numeric(df)

    df
}


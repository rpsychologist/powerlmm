
## G matrices
create_G <- function(p) {

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
    d <- simulate_data(p)

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
    if(p$cor_subject != 0) {
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
    if(p$cor_cluster != 0) {
        G5 <- X3 %*% kronecker(Diagonal(n3_cc + n3_tx), Z) %*% t(X3)
        G <- c(G, G5)
    }

    ## G6 slope cluster
    if(p$sigma_cluster_intercept != 0) {
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


#' Median bias reduction in cumulative link models
#'
#' @description  Such a function contains all the needed quantities to fit the model
#' @param x model matrix (without intercept)
#' @param y response matrix (each row with an indicator of observed category)
#' @param start vector of starting values for the algorithm
#' @param type type of adjustment of the score function
#' @param link link function
#' @param maxiter maximum number of iterations for the algorithm
#' @param tolerance value used for the stopping criterion
#'
#' @return Parameter estimates
#' @export
#' @importFrom evd qgev
#'
mbrclm <- function(x, y, start = NULL,
                   type = c("AS_median", "AS_mean", "AS_ml"),
                   link = c("logit", "probit", "cloglog"),
                   maxiter = 100,
                   tolerance = 10 ^ (-10)){


  type <- match.arg(type)
  link <- match.arg(link)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)

  key_quantity <- function(par){
    alpha <- par[1 : (q - 1)]
    alpha <- c(-Inf, alpha, Inf)
    beta <- par[q : (q + p - 1)]

    mat1 <- matrix(rep(alpha[-1], n), ncol = q, byrow = TRUE)
    mat2 <- matrix(rep(alpha[1 : q], n), ncol = q, byrow = TRUE)
    mat3 <- matrix(rep(alpha[3 : (q + 1)], n), ncol = q - 1, byrow = TRUE)

    eta_beta  <- x %*% beta

    eta_ij <-  apply(mat1, 2, function(t) t + eta_beta)
    eta_ijm1 <- apply(mat2, 2, function(t) t + eta_beta)
    eta_ijp1 <- apply(mat3, 2, function(t) t + eta_beta)

    #link <- match.arg(link)
    pfun <- switch(link,
                   logit = make.link(link)$linkinv,
                   probit = make.link(link)$linkinv,
                   cloglog = make.link(link)$linkinv)
    dfun <- switch(link,
                   logit = make.link(link)$mu.eta,
                   probit = make.link(link)$mu.eta,
                   cloglog = make.link(link)$mu.eta)
    d1fun <- switch(link,
                    logit = function(eta) dfun(eta) * (1 - 2 * pfun(eta)),
                    probit = function(eta) -eta * dfun(eta),
                    cloglog = function(eta) dfun(eta) * (1 - exp(eta)))
    d2fun <- switch(link,
                    logit = function(eta) dfun(eta) * (1 - 2 * pfun(eta))^2 - 2 * dfun(eta)^2,
                    probit = function(eta) (eta ^ 2 - 1) * dfun(eta),
                    cloglog = function(eta) d1fun(eta) * (1 + log(1 - pfun(eta))) - dfun(eta)^2/(1 - pfun(eta)))

    rip_ij <- pfun(eta_ij)
    rip_ijm1 <- pfun(eta_ijm1)
    rip_ijp1 <- pfun(eta_ijp1)

    den_ij <- dfun(eta_ij)
    den_ijm1 <- dfun(eta_ijm1)
    den_ijp1 <- dfun(eta_ijp1)

    der_den_ij <- d1fun(eta_ij)
    der_den_ijm1 <- d1fun(eta_ijm1)
    der_den_ijp1 <- d1fun(eta_ijp1)

    der2_den_ij <- d2fun(eta_ij)
    der2_den_ijm1 <- d2fun(eta_ijm1)

    pi_ij <- rip_ij - rip_ijm1
    pi_ijp1 <- rip_ijp1 - rip_ij[, -q]

    #this part try to solve the instability numerical problems
    for(j in 1 : q){
      for(i in 1 : n){
        if(pi_ij[i,j] == 0) pi_ij[i, j] <- 1e-16
        if(j < q){
          if(pi_ijp1[i,j] == 0) pi_ijp1[i, j] <- 1e-16
        }
      }
    }

    if(link == "probit" | link == "cloglog"){
      den_ij[, q] <- rep(1e-12, n)
      den_ijp1[, q - 1] <- rep(1e-12, n)
      der_den_ij[, q] <- rep(1e-12, n)
      der_den_ijp1[, q - 1] <- rep(1e-12, n)
      der_den_ijm1[, 1] <- rep(1e-12, n)
      der2_den_ij[, q] <- rep(1e-12, n)
      der2_den_ijm1[, 1] <- rep(1e-12, n)
    }

    A <- den_ij/pi_ij
    B <- den_ij[, -q]/pi_ijp1
    C <- (den_ij - den_ijm1)/pi_ij
    D <- der_den_ij/pi_ij - (den_ij/pi_ij)^2
    G <- der_den_ij[, -q]/pi_ijp1 + (den_ij[, -q]/pi_ijp1)^2
    H <- den_ij[, -q] * den_ijp1/pi_ijp1^2
    I <- der_den_ij/pi_ij - den_ij * (den_ij - den_ijm1)/pi_ij^2
    J <- der_den_ij[, -q]/pi_ijp1 - den_ij[,-q] * (den_ijp1 - den_ij[, -q])/pi_ijp1^2
    K <- (der_den_ij - der_den_ijm1)/pi_ij - (den_ij - den_ijm1)^2/pi_ij^2
    L <- der2_den_ij/pi_ij - 3 * der_den_ij * den_ij/pi_ij^2 + 2 * den_ij^3/pi_ij^3
    M <- der2_den_ij[, -q]/pi_ijp1 + 3 * der_den_ij[,-q] * den_ij[, -q]/pi_ijp1^2 + 2 * den_ij[, -q]^3/pi_ijp1^3
    N <- der_den_ij[, -q] * den_ijp1/pi_ijp1^2 + 2 * den_ij[, -q]^2 * den_ijp1/pi_ijp1^3
    O <- den_ij[, -q] * der_den_ijp1/pi_ijp1^2 - 2 * den_ij[, -q] * den_ijp1^2/pi_ijp1^3
    P <- der2_den_ij/pi_ij - der_den_ij * (3 * den_ij - den_ijm1)/pi_ij^2 + 2 * den_ij^2*(den_ij - den_ijm1)/pi_ij^3
    Q <- der2_den_ij[, -q]/pi_ijp1 + der_den_ij[, -q] * (3 * den_ij[, -q] - den_ijp1)/pi_ijp1^2 - 2 * den_ij[, -q]^2 * (den_ijp1 - den_ij[, -q])/pi_ijp1^3
    R <- der_den_ij[, -q] * den_ijp1/pi_ijp1^2 + der_den_ijp1 * den_ij[, -q]/pi_ijp1^2 - 2 * den_ij[, -q] * den_ijp1 * (den_ijp1 - den_ij[, -q])/pi_ijp1^3
    S <- der2_den_ij/pi_ij - den_ij * (der_den_ij - der_den_ijm1)/pi_ij^2 - 2 * der_den_ij * (den_ij - den_ijm1)/pi_ij^2 + 2 * den_ij * (den_ij - den_ijm1)^2/pi_ij^3
    T <- der2_den_ij[, -q]/pi_ijp1 - den_ij[, -q] * (der_den_ijp1 - der_den_ij[, -q])/pi_ijp1^2 - 2 * der_den_ij[, -q] * (den_ijp1 - den_ij[, -q])/pi_ijp1^2 + 2 * den_ij[, -q] * (den_ijp1 - den_ij[, -q])^2/pi_ijp1^3
    U <- (der2_den_ij - der2_den_ijm1)/pi_ij - 3 * (der_den_ij - der_den_ijm1) * (den_ij - den_ijm1)/pi_ij^2 + 2 * (den_ij - den_ijm1)^3/pi_ij^3

    return(list(A = A[, -q], B = B, C = C, D = D[, -q], G = G, H = H, I = I[, -q], J = J, K = K,
                L = L[,-q], M = M, N = N, O = O, P = P[, -q], Q = Q, R = R, S = S[, -q], T = T, U = U,
                pi_ij = pi_ij, pi_ijp1 = pi_ijp1))
  }

  ################
  #Score function#
  ################
  score <- function(par, fit = NULL)
  {
    res <- rep(0, length(par))
    if (is.null(fit))
    {
      fit <- key_quantity(par)
    }
    with(fit,{
      res[1 : (q - 1)] <- colSums(y[, -q] * A - y[, -1] * B)
      for(i in 1 : p){
        X <- matrix(rep(x[,i], q), ncol = q, byrow = FALSE)
        res[i + q - 1] <- sum(colSums(X * y * C))
      }
      return(res)
    })
  }


  ######################
  #Fisher's information#
  ######################
  info<- function(par, fit = NULL)
  {
    infoFisher <- matrix(0, nrow = q + p - 1, ncol = q + p - 1)
    if (is.null(fit))
    {
      fit <- key_quantity(par)
    }
    with(fit, {
      #i_aj_aj
      if(q == 2){
        infoFisher[1, 1] <- -colSums(pi_ij[, -q] * D - pi_ijp1 * G)
      } else {
        diag(infoFisher[1 : (q - 1), 1 : (q - 1)]) <- -colSums(pi_ij[, -q] * D - pi_ijp1 * G)
        #i_aj_aj+1
        for(i in 1 : (q - 1)){
          infoFisher[i, i + 1] <-  infoFisher[i + 1, i] <- -colSums(pi_ijp1 * H)[i]
          #infoFisher[i+1,i] <-  infoFisher[i,i+1]
        }
      }

      #i_aj_br
      for(i in 1 : p){
        X <- matrix(rep(x[, i], q - 1), ncol = q - 1, byrow = FALSE)
        if(q == 2){
          infoFisher[q + i - 1, 1 : (q - 1)] <- infoFisher[1 : (q - 1), q + i - 1] <- -sum(X * (pi_ij[,-q] * I - pi_ijp1 * J))
          #infoFisher[1:(q-1),q+i-1] <-  infoFisher[q+i-1,1:(q-1)]
        } else {
          infoFisher[q + i - 1, 1 : (q - 1)] <- infoFisher[1 : (q - 1), q + i - 1] <- -colSums(X * (pi_ij[, -q] * I[,-q] - pi_ijp1 * J))
          #infoFisher[1:(q-1),q+i-1] <-  infoFisher[q+i-1,1:(q-1)]
        }

        #i_br_bs
        for(j in 1 : p){
          if(j <= i){
            X1 <- matrix(rep(x[, i], q), ncol = q, byrow = FALSE)
            X2 <- matrix(rep(x[, j], q), ncol = q, byrow = FALSE)
            infoFisher[i + q - 1, j + q - 1] <- infoFisher[j + q - 1, i + q - 1] <- -sum(colSums(X1 * X2 * pi_ij * K))
          }
          #infoFisher[j+q-1,i+q-1] <- infoFisher[i+q-1,j+q-1]
        }
      }
      return(infoFisher)
    })
  }

  nu_quantities <- function(par, fit = NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantity(par)
    }

    #nu_r,st
    nu_a_aa <- array(0, dim = c(q - 1))
    if(q > 2){
      nu_a_aap1 <- array(0, dim = c(q - 2))
      nu_a_ap1ap1 <- array(0, dim=c(q - 2))
      nu_ap1_aa <- array(0, dim = c(q - 2))
      nu_ap1_aap1 <- array(0, dim = c(q - 2))
      nu_a_ap1b <- array(0, dim = c(p, q - 2))
      nu_ap1_ab <- array(0, dim = c(p, q - 2))
      nu_b_aap1 <- array(0, dim = c(p, q - 2))
      nu_a_a_ap1 <- array(0, dim = c(q - 2))
      nu_a_ap1_ap1 <- array(0, dim = c(q - 2))
      nu_a_ap1_b <- array(0, dim = c(p, q - 2))
    }

    nu_a_ab <- array(0, dim = c(p, q - 1))
    nu_b_aa <- array(0, dim = c(p, q - 1))
    nu_b_ab <- array(0, dim = c(p, q - 1, p))
    nu_a_bb <- array(0, dim = c(p, q - 1, p))
    nu_b_bb <- array(0, dim = c(p, p, p))


    nu_a_a_a <- array(0, dim = c(q - 1))
    nu_a_a_b <- array(0, dim = c(p, q - 1))
    nu_a_b_b <- array(0, dim = c(p, q - 1, p))
    nu_b_b_b <- array(0, dim = c(p, p, p))

    with(fit,{
      nu_a_aa <- colSums(pi_ij[, -q] * A * D + pi_ijp1 * B * G)
      nu_a_a_a <- -colSums(pi_ij[, -q] * L - pi_ijp1 * M) - 3 * nu_a_aa

      if(q>2){
        nu_a_aap1 <- -colSums(as.matrix(pi_ijp1[, -(q - 1)] * B[, -(q - 1)] * H[, -(q - 1)]))
        nu_a_ap1ap1 <- -colSums(as.matrix(pi_ijp1[,-(q-1)] * B[, -(q - 1)] * D[, -1]))
        nu_ap1_aa <- -colSums(as.matrix(pi_ijp1[, -(q - 1)] * A[,-1] * G[, -(q - 1)]))
        nu_ap1_aap1 <- colSums(as.matrix(pi_ijp1[, -(q - 1)] * A[, -1] * H[, -(q - 1)]))
        nu_a_a_ap1 <- -colSums(as.matrix(pi_ijp1[, -(q - 1)] * N[, -(q - 1)])) - 2 * nu_a_aap1 - nu_ap1_aa
        nu_a_ap1_ap1 <- -colSums(as.matrix(pi_ijp1[, -(q - 1)] * O[, -(q - 1)])) - 2 * nu_ap1_aap1 - nu_a_ap1ap1
      }


      for(i in 1 : p){
        X1 <- matrix(rep(x[,i], q - 1), ncol = q - 1, byrow = FALSE)

        nu_a_ab[i,] <- colSums(X1 * (pi_ij[, -q] * A * I + pi_ijp1 * B * J))
        nu_b_aa[i,] <- colSums(X1 * (pi_ij[, -q] * C[, -q] * D - pi_ijp1 * C[, -1] * G))
        nu_a_a_b[i,] <- -colSums(X1 * (pi_ij[, -q] * P - pi_ijp1 * Q[, -q])) - 2 * nu_a_ab[i,] - nu_b_aa[i,]

        if(q > 2){
          X2 <- matrix(rep(x[,i], q - 2), ncol = q - 2, byrow = FALSE)
          nu_a_ap1b[i,] <-  -colSums(X2 * pi_ijp1[, -(q - 1)] * B[, -(q - 1)] * I[, -1])
          nu_ap1_ab[i,] <-  -colSums(X2 * pi_ijp1[, -(q - 1)] * A[, -1] * J[, -(q - 1)])
          nu_b_aap1[i,] <- colSums(X2 * pi_ijp1[, -(q - 1)] * C[, -c(1, q)] * H[, -(q - 1)])
          nu_a_ap1_b[i,] <- -colSums(X2 * pi_ijp1[, -(q - 1)] * R[, -(q - 1)]) - nu_b_aap1[i,] - nu_a_ap1b[i,] - nu_ap1_ab[i,]
        }

        for(j in 1 : p){
          X1 <- matrix(rep(x[,i], q - 1), ncol = q - 1, byrow = FALSE)
          X2 <- matrix(rep(x[,j], q - 1), ncol = q - 1, byrow = FALSE)
          nu_b_ab[j, , i] <-  colSums(X1 * X2 * (pi_ij[, -q] * C[, -q] * I - pi_ijp1 * C[, -1] * J))
          nu_a_bb[j, , i] <- colSums(X1 * X2 * (pi_ij[, -q] * A * K[, -q] - pi_ijp1 * B * K[, -1]))
          nu_a_b_b[j, , i] <- -colSums(X1 * X2 * (pi_ij[, -q] * S - pi_ijp1 * T)) - 2 * nu_b_ab[j,,i] - nu_a_bb[j,,i]

          for(k in 1:p){
            X1 <- matrix(rep(x[, i], q), ncol = q, byrow = FALSE)
            X2 <- matrix(rep(x[, j], q), ncol = q, byrow = FALSE)
            X3 <- matrix(rep(x[, k], q), ncol = q, byrow = FALSE)

            nu_b_bb[i, j, k] <-  sum(colSums(X1 * X2 * X3 * pi_ij * C * K))
            nu_b_b_b[i, j, k] <- -sum(colSums(X1 * X2 * X3 * pi_ij * U)) - 3 * nu_b_bb[i, j, k]
          }
        }
      }

      if(q > 2){
        return(list(nu_a_aa = nu_a_aa,
                    nu_a_aap1 = nu_a_aap1,
                    nu_a_ap1ap1 = nu_a_ap1ap1,
                    nu_ap1_aa = nu_ap1_aa,
                    nu_ap1_aap1 = nu_ap1_aap1,
                    nu_a_ab = nu_a_ab,
                    nu_a_ap1b = nu_a_ap1b,
                    nu_ap1_ab = nu_ap1_ab,
                    nu_b_aa = nu_b_aa,
                    nu_b_aap1 = nu_b_aap1,
                    nu_b_ab = nu_b_ab,
                    nu_a_bb = nu_a_bb,
                    nu_b_bb = nu_b_bb,
                    nu_a_a_a = nu_a_a_a,
                    nu_a_a_ap1 = nu_a_a_ap1,
                    nu_a_ap1_ap1 = nu_a_ap1_ap1,
                    nu_a_a_b = nu_a_a_b,
                    nu_a_ap1_b = nu_a_ap1_b,
                    nu_a_b_b = nu_a_b_b,
                    nu_b_b_b = nu_b_b_b))}
      if(q == 2){
        return(list(nu_a_aa = nu_a_aa,
                    nu_a_ab = nu_a_ab,
                    nu_b_aa = nu_b_aa,
                    nu_b_ab = nu_b_ab,
                    nu_a_bb = nu_a_bb,
                    nu_b_bb = nu_b_bb,
                    nu_a_a_a = nu_a_a_a,
                    nu_a_a_b = nu_a_a_b,
                    nu_a_b_b = nu_a_b_b,
                    nu_b_b_b = nu_b_b_b))
      }
    })
  }

  mean_adjustment <- function(par, fit = NULL, information, invInformation){
    nu <- nu_quantities(par, fit)
    Q <- P <- array(0, dim = c((p + q - 1), (p + q - 1), (p + q - 1)))

    for(i in 1 : (q - 1)){
      Q[i, i, i] <- nu$nu_a_aa[i]
      P[i, i, i] <- nu$nu_a_a_a[i]
      Q[i, q : (q + p - 1), i] <-  Q[q : (q + p - 1), i, i] <- nu$nu_a_ab[, i]
      P[i, q : (q + p - 1), i] <-  P[q : (q + p - 1), i, i] <- nu$nu_a_a_b[, i]

      for(j in 1 : p){
        Q[q + j - 1, q : (q + p - 1), i] <- nu$nu_a_bb[, i, j]
        P[q + j - 1, q : (q + p - 1), i] <- nu$nu_a_b_b[, i, j]
      }
    }

    if(q > 2){
      for(i in 1 : (q - 2)){
        Q[i, i + 1, i] <- Q[i + 1, i, i] <- nu$nu_a_aap1[i]
        P[i, i + 1, i] <- P[i + 1, i, i] <- P[i, i, i + 1] <- nu$nu_a_a_ap1[i]
        Q[i + 1, i + 1, i] <- nu$nu_a_ap1ap1[i]
        P[i + 1, i + 1, i] <- nu$nu_a_ap1_ap1[i]
        Q[i, i, i + 1] <- nu$nu_ap1_aa[i]
        Q[i, i + 1, i + 1] <- Q[i + 1, i, i + 1] <- nu$nu_ap1_aap1[i]
        P[i, i + 1, i + 1] <- P[i + 1, i, i + 1] <- nu$nu_a_ap1_ap1[i]
        Q[i + 1, q : (q + p - 1), i] <- Q[q : (q + p - 1), i + 1, i] <- nu$nu_a_ap1b[, i]
        Q[i, q : (q + p - 1), i + 1] <- Q[q : (q + p - 1), i, i + 1] <- nu$nu_ap1_ab[, i]
        P[i, q : (q + p - 1), i + 1] <- P[q : (q + p - 1), i, i + 1] <- nu$nu_a_ap1_b[, i]
        P[i + 1, q : (q + p - 1), i] <- P[q : (q + p - 1), i + 1, i] <- nu$nu_a_ap1_b[, i]
      }
    }

    for(i in 1 : p){
      diag(Q[, , i + q - 1])[1 : q - 1] <- nu$nu_b_aa[i, ]
      diag(P[, , i + q - 1])[1 : q - 1] <- nu$nu_a_a_b[i, ]
      Q[1 : (q - 1), q : (p + q - 1), i + q - 1] <- t(nu$nu_b_ab[, , i])
      Q[q : (p + q - 1), 1 : (q - 1), i + q - 1] <- nu$nu_b_ab[, , i]
      P[1 : (q  - 1), q : (p + q - 1), i + q - 1] <- t(nu$nu_a_b_b[, , i])
      P[q : (p + q - 1), 1 : (q - 1), i + q - 1] <- nu$nu_a_b_b[, , i]
      Q[q : (q + p - 1), q : (q + p - 1), i + q - 1] <- nu$nu_b_bb[, , i]
      P[q : (q + p - 1), q : (q + p - 1), i + q - 1] <- nu$nu_b_b_b[, , i]
      if(q > 2){
        for(j in 1 : (q - 2)){
          Q[j, j + 1, i + q - 1] <- Q[j + 1, j, i + q - 1] <- nu$nu_b_aap1[i, j]
          P[j, j + 1, i + q - 1] <- P[j + 1, j, i + q - 1] <- nu$nu_a_ap1_b[i, j]
        }
      }
    }

    Astar1 <- rep(0, q + p - 1)
    for(i in 1 : (q + p - 1)){
      Astar1[i] <- 0.5 * sum(diag(invInformation %*% (Q[, , i] + P[, , i])))
    }
    Astar1
  }

  median_adjustment <- function(par, fit = NULL, information, invInformation){
    nu <- nu_quantities(par, fit)
    Q <- P <- array(0, dim = c((p + q - 1), (p + q - 1), (p + q - 1)))

    for(i in 1 : (q - 1)){
      Q[i, i, i] <- nu$nu_a_aa[i]
      P[i, i, i] <- nu$nu_a_a_a[i]
      Q[i, q : (q + p - 1), i] <-  Q[q : (q + p - 1), i, i] <- nu$nu_a_ab[, i]
      P[i, q : (q + p - 1), i] <-  P[q : (q + p - 1), i, i] <- nu$nu_a_a_b[, i]

      for(j in 1 : p){
        Q[q + j - 1, q : (q + p - 1), i] <- nu$nu_a_bb[, i, j]
        P[q + j - 1, q : (q + p - 1), i] <- nu$nu_a_b_b[, i, j]
      }
    }

    if(q > 2){
      for(i in 1 : (q - 2)){
        Q[i, i + 1, i] <- Q[i + 1, i, i] <- nu$nu_a_aap1[i]
        P[i, i + 1, i] <- P[i + 1, i, i] <- P[i, i, i + 1] <- nu$nu_a_a_ap1[i]
        Q[i + 1, i + 1, i] <- nu$nu_a_ap1ap1[i]
        P[i + 1, i + 1, i] <- nu$nu_a_ap1_ap1[i]
        Q[i, i, i + 1] <- nu$nu_ap1_aa[i]
        Q[i, i + 1, i + 1] <- Q[i + 1, i, i + 1] <- nu$nu_ap1_aap1[i]
        P[i, i + 1, i + 1] <- P[i + 1, i, i + 1] <- nu$nu_a_ap1_ap1[i]
        Q[i + 1, q : (q + p - 1), i] <- Q[q : (q + p - 1), i + 1, i] <- nu$nu_a_ap1b[, i]
        Q[i, q : (q+p-1),i+1] <- Q[q : (q + p - 1), i, i + 1] <- nu$nu_ap1_ab[, i]
        P[i, q : (q+p-1),i+1] <- P[q : (q + p - 1), i, i + 1] <- nu$nu_a_ap1_b[, i]
        P[i + 1, q:(q+p-1),i] <- P[q : (q + p - 1), i + 1, i] <- nu$nu_a_ap1_b[, i]
      }
    }

    for(i in 1 : p){
      diag(Q[, , i + q - 1])[1 : q - 1] <- nu$nu_b_aa[i, ]
      diag(P[, , i + q - 1])[1 : q - 1] <- nu$nu_a_a_b[i, ]
      Q[1 : (q - 1), q : (p + q - 1), i + q - 1] <- t(nu$nu_b_ab[, , i])
      Q[q : (p + q - 1), 1 : (q - 1), i + q - 1] <- nu$nu_b_ab[, , i]
      P[1 : (q - 1), q : (p + q - 1), i + q - 1] <- t(nu$nu_a_b_b[, , i])
      P[q : (p + q - 1), 1 : (q - 1), i + q - 1] <- nu$nu_a_b_b[, , i]
      Q[q : (q + p - 1), q : (q + p - 1), i + q - 1] <- nu$nu_b_bb[, , i]
      P[q : (q + p - 1), q : (q + p - 1), i + q - 1] <- nu$nu_b_b_b[, , i]
      if(q > 2){
        for(j in 1 : (q - 2)){
          Q[j, j + 1, i + q - 1] <- Q[j + 1, j, i + q - 1] <- nu$nu_b_aap1[i, j]
          P[j, j + 1, i + q - 1] <- P[j + 1, j, i + q - 1] <- nu$nu_a_ap1_b[i, j]
        }
      }
    }

    Astar1 <- rep(0, q + p - 1)
    for(i in 1 : (q + p - 1)){
      Astar1[i] <- 0.5 * sum(diag(invInformation %*% (Q[, , i] + P[, , i])))
    }

    h <- array(0, dim = c((p + q - 1), (p + q - 1), (p + q - 1)))
    for(i in 1 : (q + p - 1)){
      h[, , i] <- invInformation[, i] %*% t(invInformation[, i])/diag(invInformation)[i]
    }
    F_theta <- rep(0, q + p - 1)
    for(i in 1 : (q + p - 1)){
      F_tilde <- rep(0, q + p - 1)
      for(j in 1 : (q + p - 1)){
        F_tilde[j] <-  sum(diag(h[, , i] %*% (P[, , j]/3 + Q[, , j]/2)))
      }
      F_theta[i] <- t(invInformation[,i]) %*% F_tilde
    }

    Astar1 - information %*% F_theta
  }


  adjustment_function <- switch(type,
                                AS_mean = mean_adjustment,
                                AS_median = median_adjustment,
                                AS_ml = function(par, ...) 0)


  if(is.null(start)){
    start <- switch(link,
                    logit = c(qlogis((1 : (q - 1)) / q), rep(0, p)) ,
                    probit = c(qnorm((1 : (q - 1)) / q), rep(0, p)) ,
                    cloglog = c(-qgev((1 : (q - 1)) / q, lower.tail = FALSE), rep(0, p)))
  }

  par <- start
  for (k in 1 : maxiter){
    fit <- key_quantity(par)
    information <- info(par, fit)
    inv <- try(solve(information), TRUE)
    if(!is.character(inv)) info_inv <- inv
    adjscore <- score(par,fit) +  adjustment_function(par, fit, information, info_inv)
    if(all(abs(adjscore) <= tolerance)) break
    par <- par + info_inv %*% adjscore
  }
  if (k < maxiter) converged <- TRUE
  else converged <- FALSE

  res <- list(alphas = par[1 : (q - 1)],
              betas = par[q : (q + p - 1)],
              coefficients = par,
              std.err = sqrt(diag((info_inv))),
              iter = k,
              start = start,
              convergence = converged,
              vcov = info_inv,
              link = link,
              type = type)
  class(res) <- c("mbrclm")
  res
}





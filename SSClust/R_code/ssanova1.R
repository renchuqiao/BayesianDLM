##############################################################
# Smoothing spline method used in                            #
# Ma, Castillo-Davis, Zhong, and Liu (2006)                  #
# "A data-driven clustering method for Time Course Gene      #
# Expression Data"  Nucleic Acids Research, 34(4), 1261-1269.# 
#                                                            #
# Author: Chong Gu  adapted by Ping Ma                       #
# Email:  chong@stat.purdue.edu, pingma@uiuc.edu             #
##############################################################



mkrandom <- function (formula, randata) {
    attach(randata)
    form.wk <- terms.formula(formula)[[2]]
    if (!("|" %in% strsplit(deparse(form.wk), "")[[1]])) 
        stop("gss error in mkran: missing | in grouping formula")
    term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
    z2.wk <- eval(parse(text = term.wk[2]))  
    z <- NULL
    lvl.z2 <- levels(z2.wk)
    for (i in lvl.z2) z <- cbind(z, as.numeric(z2.wk == i))
    init <- -1
    env <- length(levels(z2.wk))
    fun <- function(zeta, env) diag(10^(-zeta), env)
    sigma <- list(fun = fun, env = env)    
#    detach(randata)	
    list(z = z, sigma = sigma, init = init)
}

chol<- function (x, pivot = FALSE, LINPACK = pivot) 
{
    if (is.complex(x)) 
        stop("complex matrices not permitted at present")
    else if (!is.numeric(x)) 
        stop("non-numeric argument to 'chol'")
    if (is.matrix(x)) {
        if (nrow(x) != ncol(x)) 
            stop("non-square matrix in 'chol'")
        n <- nrow(x)
    }
    else {
        if (length(x) != 1) 
            stop("non-matrix argument to 'chol'")
        n <- as.integer(1)
    }
    if (!pivot && !LINPACK) 
        return(.Call("La_chol", as.matrix(x), PACKAGE = "base"))
    if (!is.double(x)) 
        storage.mode(x) <- "double"
    if (pivot) {
        xx <- x
        xx[lower.tri(xx)] <- 0
        z <- .Fortran("dchdc",
                      x = xx,
                      n,
                      n,
                      double(n),
                      piv = integer(n),
                      as.integer(pivot),
                      rank = integer(1L),
                      DUP = FALSE, PACKAGE = "lapack")
        if (z$rank < n) 
            if (!pivot) 
                stop("matrix not positive definite")
        robj <- z$x
        if (pivot) {
            attr(robj, "pivot") <- z$piv
            attr(robj, "rank") <- z$rank
            if (!is.null(cn <- colnames(x))) 
                colnames(robj) <- cn[z$piv]
        }
        robj
    }
    else {
        z <- .Fortran("chol", x = x, n, n, v = matrix(0, nr = n, 
            nc = n), info = integer(1), DUP = FALSE, PACKAGE = "base")
        if (z$info) 
            stop("non-positive definite matrix in 'chol'")
        z$v
    }
}



sspreg1 <- function (s, r, q, y, method, alpha, varht, random, smpar=NULL) 
{
    qr.trace <- FALSE
    if ((alpha < 0) & (method %in% c("u", "v", "n"))) 
        qr.trace <- TRUE
    alpha <- abs(alpha)
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) 
            nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) 
        nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    cv <- function(lambda) {
        if (is.null(random)) {
            q.wk <- 10^(lambda + theta) * q
            q.wk=(q.wk+t(q.wk))/2
          }
        else {
            q.wk <- matrix(0, nxiz, nxiz)
            q.wk[1:nxi, 1:nxi] <- 10^(lambda[1] + theta) * q
            q.wk[(nxi + 1):nxiz, (nxi + 1):nxiz] <- 10^(2 * ran.scal) * 
                random$sigma$fun(lambda[-1], random$sigma$env)
            q.wk=(q.wk+t(q.wk))/2
        }
        if (qr.trace) {
            qq.wk <- chol(q.wk+diag(.Machine$double.eps, dim(q.wk)[1]), pivot = TRUE)
            sr <- cbind(s, 10^theta * r[, attr(qq.wk, "pivot")])
            sr <- rbind(sr, cbind(matrix(0, nxiz, nnull), qq.wk))
            sr <- qr(sr, tol = 0)
            rss <- mean(qr.resid(sr, c(y, rep(0, nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs, ]^2)/nobs
            if (method == "u") 
                score <- rss + alpha * 2 * varht * trc
            if (method %in% c("v", "n")) 
                score <- rss/(1 - alpha * trc)^2
            alpha.wk <- max(0, log.la0 - lambda[1] - 5) * (3 - 
                alpha) + alpha
            alpha.wk <- min(alpha.wk, 3)
            if (alpha.wk > alpha) {
                if (method == "u") 
                  score <- score + (alpha.wk - alpha) * 2 * varht * 
                    trc
                if (method %in% c("v","n")) 
                  score <- rss/(1 - alpha.wk * trc)^2
            }
            if (return.fit) {
                z <- .Fortran("reg", as.double(cbind(s, 10^theta * 
                  r)), as.integer(nobs), as.integer(nnull), as.double(q.wk), 
                  as.integer(nxiz), as.double(y), as.integer(switch(method, 
                    u = 1, v = 2, m = 3, n=2)), as.double(alpha), 
                  varht = as.double(varht), score = double(1), 
                  dc = double(nn), as.double(.Machine$double.eps), 
                  chol = double(nn * nn), double(nn), jpvt = as.integer(c(rep(1, 
                    nnull), rep(0, nxiz))), wk = double(nobs + 
                    nnull + nz), rkv = integer(1), info = integer(1), 
                  PACKAGE = "gss")[c("score", "varht", "dc", 
                  "chol", "jpvt", "wk", "rkv", "info")]
                z$score <- score
                assign("fit", z[c(1:5, 7)], inherit = TRUE)
            }
        }
        else {
            z <- .Fortran("reg", as.double(cbind(s, 10^theta * 
                r)), as.integer(nobs), as.integer(nnull), as.double(q.wk), 
                as.integer(nxiz), as.double(y), as.integer(switch(method, 
                  u = 1, v = 2, m = 3, n=2)), as.double(alpha), varht = as.double(varht), 
                score = double(1), dc = double(nn), as.double(.Machine$double.eps), 
                chol = double(nn * nn), double(nn), jpvt = as.integer(c(rep(1, 
                  nnull), rep(0, nxiz))), wk = double(nobs + 
                  nnull + nz), rkv = integer(1), info = integer(1), 
                PACKAGE = "gss")[c("score", "varht", "dc", "chol", 
                "jpvt", "wk", "rkv", "info")]
            if (z$info) 
                stop("gss error in ssanova: evaluation of GML score fails")
            assign("fit", z[c(1:5, 7)], inherit = TRUE)
            score <- z$score
            alpha.wk <- max(0, log.la0 - lambda[1] - 5) * (3 - 
                alpha) + alpha
            alpha.wk <- min(alpha.wk, 3)
            if (alpha.wk > alpha) {
                if (method == "u") 
                  score <- score + (alpha.wk - alpha) * 2 * varht * 
                    z$wk[2]
                if (method %in% c("v", "n")) 
                  score <- z$wk[1]/(1 - alpha.wk * z$wk[2])^2
            }
        }
        score
    }

    cv.wk <- function(lambda) cv.scale * cv(lambda) + cv.shift
    tmp <- sum(r^2)
    if (is.null(s)) 
        theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp * nxi)/2
    log.la0 <- log10(tmp/sum(diag(q))) + theta

    if (!is.null(random)) {
        ran.scal <- theta - log10(sum(random$z^2)/nz/tmp * nxi)/2
        r <- cbind(r, 10^(ran.scal - theta) * random$z)
    }
    else ran.scal <- NULL
    return.fit <- FALSE
    fit <- NULL
    if (is.null(random)) 
        la <- log.la0
    else la <- c(log.la0, random$init)

    if (method %in% c("u", "v", "m")){
        
    if (length(la) - 1) {
        counter <- 0
        tmp <- abs(cv(la))
        cv.scale <- 1
        cv.shift <- 0
        if (tmp < 1 & tmp > 10^(-4)) {
            cv.scale <- 10/tmp
            cv.shift <- 0
        }
        if (tmp < 10^(-4)) {
            cv.scale <- 10^2
            cv.shift <- 10
        }
        repeat {
            zz <- nlm(cv.wk, la, stepmax = 1, ndigit = 7)
            if (zz$code <= 3) 
                break
            la <- zz$est
            counter <- counter + 1
            if (counter >= 5) {
                break
            }
        }
    }
    else {
        repeat {
            mn <- la - 1
            mx <- la + 1
            zz <- nlm0(cv, c(mn, mx))
            if (min(zz$est - mn, mx - zz$est) >= 0.001) 
                break
            else la <- zz$est
        }
    }
    return.fit <- TRUE
    jk1 <- cv(zz$est)
  }
 else{
     return.fit <- TRUE
     zz=NULL
     if (is.null(smpar)) {
       zz$est=la
       jk1 <- cv(la)
     }
     else {
       zz$est=smpar
       jk1 <- cv(smpar)
     }  
  }   
    if (is.null(random)) {
        q.wk <- 10^theta * q
        q.wk=(q.wk+t(q.wk))/2
      }
    else {
        q.wk <- matrix(0, nxiz, nxiz)
        q.wk[1:nxi, 1:nxi] <- 10^theta * q
        q.wk[(nxi + 1):nxiz, (nxi + 1):nxiz] <- 10^(2 * ran.scal - 
            zz$est[1]) * random$sigma$fun(zz$est[-1], random$sigma$env)
        q.wk=(q.wk+t(q.wk))/2
    }
            qq.wk <- chol(q.wk+diag(.Machine$double.eps, dim(q.wk)[1]), pivot = TRUE)
            sr <- cbind(s, 10^theta * r[, attr(qq.wk, "pivot")])
            sr <- rbind(sr, cbind(matrix(0, nxiz, nnull), qq.wk))
            sr <- qr(sr, tol = 0)
            rss <- mean(qr.resid(sr, c(y, rep(0, nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs, ]^2)/nobs

 
    qinv <- eigen(q.wk, TRUE)
    se.aux <- t(cbind(s, 10^theta * r)) %*% (10^theta * r) %*%
        qinv$vec
    c <- fit$dc[nnull + (1:nxi)]
    if (nnull)
        d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz)
        b <- 10^(ran.scal) * fit$dc[nnull + nxi + (1:nz)]
    else b <- NULL    
    c(list(method = method, theta = theta, trc=trc*nobs, mgcv=jk1, ran.scal = ran.scal, 
        c = c, d = d, b = b, nlambda = zz$est[1], zeta = zz$est[-1]), 
        fit[-3], list(qinv = qinv, se.aux = se.aux))
}


ssanova1 <- function (formula, type = NULL, data = list(), weights, subset, 
    offset, na.action = na.omit, partial = NULL, method = "v", smpar=NULL,
    alpha = 1.4, varht = 1, id.basis = NULL, nbasis = NULL, seed = NULL, 
    random = NULL) 
{
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- mf$partial <- mf$smpar <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$random <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    wt <- model.weights(mf)
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis)) 
            nbasis <- max(30, ceiling(10 * nobs^(2/9)))
        if (nbasis >= nobs) 
            nbasis <- nobs
        if (!is.null(seed)) 
            set.seed(seed)
        id.basis <- sample(nobs, nbasis, prob = wt)
    }
    else {
        if (max(id.basis) > nobs | min(id.basis) < 1) 
            stop("gss error in ssanova: id.basis out of range")
        nbasis <- length(id.basis)
    }
    term <- mkterm(mf, type)
    if (!is.null(random)) {
        if (class(random) == "formula")
         random <- mkrandom(random, data)
#            random <- mkran(random, data)
    }
    s <- r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label == "1") {
            s <- cbind(s, rep(1, len = nobs))
            next
        }
        x <- mf[, term[[label]]$vlist]
        x.basis <- mf[id.basis, term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) s <- cbind(s, phi$fun(x, nu = i, 
                env = phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq + 1
                r <- array(c(r, rk$fun(x, x.basis, nu = i, env = rk$env, 
                  out = TRUE)), c(nobs, nbasis, nq))
            }
        }
    }
    if (is.null(r)) 
        stop("gss error in ssanova: use lm for models with only fixed effects")
    else q <- r[id.basis, , , drop = FALSE]
    if (!is.null(partial)) {
        if (is.vector(partial)) 
            partial <- as.matrix(partial)
        if (dim(partial)[1] != dim(mf)[1]) 
            stop("gss error in ssanova: partial data are of wrong size")
        term$labels <- c(term$labels, "partial")
        term$partial <- list(nphi = dim(partial)[2], nrk = 0, 
            iphi = ifelse(is.null(s), 0, dim(s)[2]) + 1)
        s <- cbind(s, partial)
        mf$partial <- partial
    }
#    if (qr(s)$rank < dim(s)[2]) 
#        stop("gss error in ssanova: fixed effects are linearly dependent")
    y <- model.response(mf, "numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels, "offset")
        term$offset <- list(nphi = 0, nrk = 0)
        y <- y - offset
    }
    if (!is.null(wt)) {
        wt <- sqrt(wt)
        y <- wt * y
        s <- wt * s
        r <- wt * r
        if (!is.null(random)) 
            random$z <- wt * random$z
    }
#    if (qr(s)$rank < dim(s)[2]) 
#        stop("gss error in ssanova: fixed effects are linearly dependent")
    if (nq == 1) {
        r <- r[, , 1]
        q <- q[, , 1]
        z <- sspreg1(s, r, q, y, method, alpha, varht, random, smpar)
    }
    else z <- mspreg1(s, r, q, y, method, alpha, varht, random)
    desc <- NULL
    for (label in term$labels) desc <- rbind(desc, as.numeric(c(term[[label]][c("nphi", 
        "nrk")])))
    desc <- rbind(desc, apply(desc, 2, sum))
    rownames(desc) <- c(term$labels, "total")
    colnames(desc) <- c("Unpenalized", "Penalized")
    obj <- c(list(call = match.call(), mf = mf, terms = term, 
        desc = desc, alpha = alpha, id.basis = id.basis, random = random), 
        z)
    class(obj) <- c("ssanova")
    obj
}








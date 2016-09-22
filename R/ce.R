#fms2-----------------------------------
##' Covariance Estimation by modified PCA
##' 
##' @param x matrix or dataframe of timeseries returns
##' @param weight weights in estimation
##' @param center flag to center 
##' @param frac.var controls auto-selection of number of factord
##' @param iter.max maximum number of iterations
##' @param nfac.miss number of factors to estimate if data is missing
##' @param full.min minimum acceptable number of NA-free columns
##' @param reg.min minimum dates to do regression
##' @param sd.min minimum dates to estimate vol
##' @param quan.sd missing vol assigned this quantile
##' @param tol estimation tolerance
##' @param zero.load flag to use zero loadings for columns with missing
##' @param range.factors range of factors to estimate, as a function of valid data length
##' @param lambda exponent on eigenvalue for shrinkage
##' @param minunique minimum uniqueness
##' @param shrinkb shrinkage for factor 1
##' @param shrinkv shrinkage for vol
##' @param shrinkr shrinkage for regressed loadings
##' @return list(loadings fmp hpl method full uniqueness sdev qua weight call)
##' @section Details: more detail on the underlying algorithm may be found in documentation for BurStFin
##' @author Giles Heywood from Pat Burns original
##' @export fms2
`fms2` <-function (x, weight=seq(1, 3, length=nobs), center=TRUE, 
    frac.var=.5, iter.max=1, nfac.miss=1, full.min=20, reg.min=40, 
    sd.min=20, quan.sd=.90, tol=1e-3, zero.load=FALSE, 
    range.factors=c(20,20),
    lambda=0,               #lambda=0 for PB original; lambda=0.15 OK
    minunique=0.02,         #minunique=0 for PB original
    shrinkb=0.3,            #factor 1 shrinkage to mean, shrink=1 for equal factor 1 loadings
    shrinkv=shrinkb,        #vol shrinkage to mean
    shrinkr=0.9,             #regressed loadings shrinkage to mean
    ...
    ) {     
    fun.copyright <- "Placed in the public domain in 2006 by Burns Statistics"
    fun.version <- "factor.model.stat 006"
    stopifnot(lambda>=0 & lambda<=1)
    stopifnot(shrinkb>=0 && shrinkv>=0 && shrinkr>=0 && shrinkb<=1 && shrinkv<=1 && shrinkr<=1)
    `subfun.ssd` <- function(z, weight, sd.min) {
        stopifnot(all.equal(sum(weight),1))
        nas <- is.na(z)
        if(any(nas)) {
            if(sum(!nas) < sd.min) return(NA)
            sum(weight[!nas] * z[!nas]^2) / sum(weight[!nas])
        } else {
            sum(weight * z^2)
        }
    }
    if(is.data.frame(x)) {
        x <- as.matrix(x)
    } else {
        if(!is.matrix(x)) stop("x needs to be a matrix")
    }
    if(!is.numeric(x)) stop("x needs to be numeric")
    x[!is.finite(x)] <- NA

    xna <- is.na(x)
    allmis <- rowSums(xna) == ncol(x)
    if(any(allmis)) {
        x <- x[!allmis, , drop=FALSE]
        xna <- is.na(x)
    }
    num.mis <- colSums(xna)

    if(any(num.mis > 0)) {
        if(sum(num.mis == 0) < full.min) 
            stop("not enough columns without missing values")
        if(!length(dimnames(x)[[2]]) | length(dimnames(x)[[2]])!=unique(length(dimnames(x)[[2]]))) 
            stop("x needs unique column names when missing values exist")

        max.miss <- max(num.mis)
        lnfm <- length(nfac.miss)
        if(lnfm == 0) stop("nfac.miss must have positive length")
        nfac.miss <- round(nfac.miss)
        if(any(nfac.miss < 0)) 
            stop("negative values in nfac.miss")
        if(lnfm < max.miss) {
            nfac.miss <- c(nfac.miss, rep(nfac.miss[lnfm],
                max.miss - lnfm))
        }
    }

    nassets <- ncol(x)
    nobs <- nrow(x)
    if(length(weight) != nobs) {
        if(length(weight) == nobs + sum(allmis)) {
            weight <- weight[!allmis]
        } else {
            stop("weight vector is the wrong length")
        }
    }
    if(any(weight < 0)) stop("negative weights not allowed")
    weight <- weight / sum(weight)

    if(is.logical(center)) {
        if(center) {
            center <- colMeans(x, na.rm=TRUE)
        } else {
            center <- rep(0, nassets)
        }
    } else if(length(center) != nassets) stop("center is the wrong length")
    
    x <- sweep(x=x, MARGIN=2, STATS=center, FUN="-")

    sdev <- sqrt(apply(X=x, MARGIN=2, FUN=subfun.ssd, weight=weight, sd.min=sd.min))

    if(any(sdev <= 0, na.rm=TRUE)) {
        stop(paste(sum(sdev <= 0, na.rm=TRUE), 
            "asset(s) with constant returns"))
    }
    if(any(is.na(sdev))) {
        sdev[is.na(sdev)] <- quantile(x=sdev, probs=quan.sd, na.rm=TRUE)
    }
    x <- scale(x=x, scale=sdev, center=FALSE)
    xw <- sqrt(weight) * x
    fullxw <- xw[, num.mis == 0, drop=FALSE]
    nrank <- min(nrow(fullxw),ncol(fullxw))
    stopifnot(all.equal(as.numeric(diag(cov.wt(x[, num.mis == 0, drop=FALSE],wt=weight,method="ML",center=FALSE)$cov)),rep(1,ncol(fullxw)))) #x is unit weighted variance
    fx.svd <- svd(x=fullxw, nu=0)
    cumvar <- cumsum(fx.svd$d^2) / sum(fx.svd$d^2)
    nfac <- sum(cumvar < frac.var) + 1
    if(nfac > max(range.factors)) nfac <- max(range.factors)
    if(nfac < min(range.factors)) nfac <- min(range.factors)
    if(nfac > length(cumvar)) nfac <- length(cumvar)
    fseq <- 1:nfac
    loadings <- scale(x=fx.svd$v, scale=1./fx.svd$d, center=FALSE)[, fseq, drop=FALSE]
    method <- matrix("D",nassets,nfac,dimnames=list(dimnames(x)[[2]], NULL))
    method[num.mis == 0, ] <- "S"
    if(iter.max > 0) {
        method[num.mis == 0,] <- "P"
        cormat <- cov.wt(
                        x=x[, num.mis == 0, drop=FALSE], 
                        wt = weight,                        #weighted estimation
                        cor = TRUE, 
                        center = FALSE,
                        method = "unbiased"
                        )$cor
        uniqueness <- 1 - rowSums(loadings^2)
        uniqueness[uniqueness < 0] <- 0
        uniqueness[uniqueness > 1] <- 1
        start <- uniqueness
        converged <- FALSE
            for(i in 1:iter.max) {
                cor.red <- cormat
                diag(cor.red) <- diag(cor.red) - uniqueness
                t.eig <- eigen(x=cor.red)
                t.val <- t.eig$value[fseq]
                t.val[t.val < 0] <- 0                       #could want to know how many of these there are... 
                loadings1 <- scale(x=t.eig$vectors[,fseq],scale=1./sqrt(t.val),center=FALSE)
                loadings <- scale(x=loadings1,scale=(t.val[1]/t.val)**lambda,center=FALSE) #shrink factor 2-k loadings
                loadings[,1] <- loadings[,1]*(1-shrinkb)+mean(loadings[,1]*shrinkb,na.rm=TRUE) #shrink factor 1 loadings
                comm <- rowSums(loadings^2)
                if(any(comm > 1-minunique)) {   # adjust loadings where communalities exceed threshold
                    toobig <- comm > 1-minunique
                    loadings[toobig,] <- sweep(x=loadings[toobig,,drop=FALSE],STATS=sqrt(comm[toobig]/(1-minunique)),FUN="/",MARGIN=1)
                    comm[toobig] <- rowSums(loadings[toobig,,drop=FALSE]^2)
                }
                uniqueness <- 1 - rowSums(loadings^2)

                if(all(abs(uniqueness - start) < tol)) {
                        converged <- TRUE
                        break
                }
                start <- uniqueness
            }
        }
    dimnames(loadings) <- list(dimnames(fullxw)[[2]], NULL)
    if(sum(uniqueness==0)>0) {stop("Heywood cases")}

    zeros <- rep(0,ncol(loadings))
    Dmat <- diag(uniqueness*sdev[num.mis == 0]**2)
    dvec <- rep(0,nrow(loadings))
    Amat <- loadings
    meq <- ncol(Amat)
    gma <- delta <- loadings*NA
    for(j in 1:ncol(loadings)) {
        bvec <- zeros
        bvec[j] <- 1
        gma[,j] <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq, factorized=FALSE)$solution
        if(all(c(1,-1)%in%sign(loadings[,j]))) {
            delta[,j] <-  solve.QP(Dmat=Dmat, dvec=dvec, Amat=cbind(loadings,abs(loadings)[,j]), bvec=c(zeros,1), meq=meq+1, factorized=FALSE)$solution 
        } else {
            delta[,j] <-  0
        }
    }
    fmp <- hpl <- matrix(0,nrow=nassets,ncol=nfac,dimnames=list(dimnames(x)[[2]], NULL))
    fmp[rownames(gma),] <- gma     #this relies on x having colnames
    hpl[rownames(delta),] <- delta
    qua <- array(0, c(nassets, 1))       #loading on factor1^2
    dimnames(qua) <- list(dimnames(x)[[2]], NULL)
    unwscores <- x[, num.mis == 0, drop=FALSE]%*%gma   
    imissing <- (1:nassets)[num.mis > 0 & nobs - num.mis > reg.min]
    icomplete <- (1:nassets)[num.mis == 0]
    if(any(num.mis>0)) {
        # calculate loadings for columns with NAs, by regression
        floadings <- loadings
        if(zero.load) {
            loadings <- array(0, c(nassets, nfac))
        } else {
            meanload <- colMeans(floadings)
            loadings <- t(array(meanload, c(nfac, nassets)))
        }
        dimnames(loadings) <- list(dimnames(x)[[2]], NULL)
        loadings[dimnames(floadings)[[1]], ] <- floadings
        #scores <- fullxw %*% floadings #this commented out because we use weighted estimation later
        #dsquare <- fx.svd$d[1:nfac]^2
        nfac.miss[nfac.miss > nfac] <- nfac
        priorloadings <- apply(loadings,2,mean,na.rm=TRUE)
        for(i in imissing) {
            t.nfac <- nfac.miss[ num.mis[i] ]
            if(t.nfac == 0) next
            t.okay <- !is.na(xw[, i])
            t.seq <- 1:t.nfac
            #t.load <- lsfit(xw[t.okay, i], scores[t.okay, t.seq],intercept=FALSE)$coef / dsquare[t.seq]
            regs <- lsfit(                   #amended version, identical for complete data, if scores calculated from loadings (not gma)
                        y=x[t.okay, i, drop=FALSE], 
                        x=addq(unwscores[t.okay, t.seq, drop=FALSE]), #add regressor for qua
                        wt = weight[t.okay], 
                        intercept = FALSE, 
                        tolerance = 1e-07,
                        yname = NULL
                        )
            loadings[i, t.seq] <- regs$coefficients[t.seq]*(1-shrinkr)+priorloadings[t.seq]*shrinkr
            qua[i,1] <- regs$coefficients[t.nfac+1]
            method[i, t.seq] <- "R"
            NULL
        }
    }
    scores1 <- addq(unwscores)
    scores1w <- 
        sweep(x=scores1, 
        MARGIN=1, 
        STATS=weight, 
        FUN="*"
        )
    xtxixt <- solve(t(scores1w)%*%scores1)[nfac+1,,drop=FALSE]%*%t(scores1w) #WLS
    qua[icomplete,1] <- xtxixt%*%x[,icomplete,drop=FALSE]   #the following is equivalent
    comm <- rowSums(loadings^2)
    if(any(comm > 1-minunique)) {   # adjust loadings where communalities exceed threshold
        toobig <- comm > 1-minunique
        loadings[toobig,] <- sweep(x=loadings[toobig,,drop=FALSE],STATS=sqrt(comm[toobig]/(1-minunique)),FUN="/",MARGIN=1)
        comm[toobig] <- rowSums(loadings[toobig,,drop=FALSE]^2)
    }
    pol <- sign(apply(loadings,2,sum))
    ans <- list(
                loadings=sweep(loadings,2,pol,"*"),     #loadings for correlation
                fmp=sweep(fmp,2,pol,"*"),               #fmp in units of standardised data, full=TRUE columns only
                hpl=hpl,
                method=method,                          #method for estimation of each loading
                full=as.matrix(num.mis==0),             #'complete data' (after holiday elimination) flag
                uniqueness=as.matrix(1 - comm),
                sdev=as.matrix(sdev)*(1-shrinkv)+mean(sdev)*shrinkv,
                qua=qua,                                #the regressed-in coefficent on centered factor1**2
                weight=weight,
                call=match.call()
                )
    colnames(ans$loadings) <- psz("loadings",1:ncol(ans$loadings))
    colnames(ans$fmp) <- psz("fmp",1:ncol(ans$fmp))
    colnames(ans$hpl) <- psz("hpl",1:ncol(ans$hpl))
    colnames(ans$method) <- psz("method",1:ncol(ans$method))
    colnames(ans$full) <- "full"
    colnames(ans$uniqueness) <- "uniqueness"
    colnames(ans$sdev) <- "sdev"
    colnames(ans$qua) <- "qua"
    class(ans) <- "ce"
    ans
}
`addq` <- function(x)
{
    stopifnot(is.matrix(x))
    cbind(x,(x[,1,drop=FALSE]-mean(x[,1,drop=FALSE],na.rm=TRUE))**2)
}
#vcvce-----------------------------------
##' Extractor for covariance matrix
##' 
##' @param x object of class ce
##' @param out component to return: M=factor1, S=factors 2:k, R=residual, T=M+S+R
##' @param units variance or correlation
##' @param po row/column identifiers
##' @return named list of the components specified in argument: out
##' @section Details: po is a vector of identifiers, a subset of the original data columnames
##' @author Giles Heywood 
##' @family extractors
##' @export vcvce
`vcvce` <- function(x,out=c("M","S","R","T"),units=c("variance","correlation"),po=buice(x)) {
    stopifnot(is(x,"ce"))
    stopifnot(all(out %in% c("M","S","R","T")))
    stopifnot(all(po %in% buice(x)))
    stopifnot(!any(is.na(po)))
    units <- match.arg(units)
    r <- list(M=NULL,S=NULL,R=NULL,T=NULL)
    if("M" %in% out) r$M <- tcrossprod(x$loadings[po,1,drop=FALSE])  
    if("S" %in% out) r$S <- tcrossprod(x$loadings[po,-1,drop=FALSE])
    if("R" %in% out) r$R <- diag(as.numeric(x$uniqueness[po,]))
    if(all(c("M","S","R","T") %in% out)) r$T <- r$M+r$S+r$R
    if(units=="variance") {
        v <- crossprod(t(x$sdev[po,]))
        for(i in out) r[[i]] <- r[[i]]*v
    }
    r
}
#ldggmace-----------------------------------
##' Extractor for loadings * fmp, (colnames=T, rownames=T)
##' 
##' @param x object of class ce
##' @return named list(M,S,R)
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export ldggmace
`ldggmace` <- function(x) {
    stopifnot(is(x,"ce"))
    r <- list(M=NULL,S=NULL,R=NULL)
    r$M <- ldgce(x)[,1,drop=FALSE] %*% t(fmpce(x)[,1,drop=FALSE])
    r$S <- ldgce(x)[,-1,drop=FALSE] %*% t(fmpce(x)[,-1,drop=FALSE])
    r$R <- diag(length(x$uniqueness))-r$M-r$S
    r
}
#vfuce-----------------------------------
##' Convenience wrapper on vcvce(x)$T -the total covariance estimate for columns with full history
##' 
##' @param x object of class ce
##' @return VCV matrix
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export vfuce
`vfuce` <- function(x){
    stopifnot(is(x,"ce"))
    vcvce(x)$T[x$full,x$full]
}
#ldgce-----------------------------------
##' Extractor for loadings
##' 
##' @param x object of class ce
##' @return loadings matrix
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export ldgce
`ldgce` <- function(x,type=c("fmp","hpl")){
    type <- match.arg(type)
    stopifnot(is(x,"ce"))
    if(type=="fmp") {ldg <- x$loadings} else {ldg <- abs(x$loadings)}    
    ldg*as.numeric(x$sdev)
}
`pococe` <- function(x){
    stopifnot(is(x,"ce"))
    x$poco
}
#quace-----------------------------------
##' Extractor for score on quadratic component
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export quace
`quace` <- function(x,ret){
    stopifnot(is(x,"ce"))
    mz((scoce(x,ret)[,1,drop=FALSE]**2) %*% t(x$qua)*as.numeric(x$sdev))
}
#devce-----------------------------------
##' Extractor for score on: deviation of rem from x-section mean
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export devce
`devce` <- function(x,ret){
    stopifnot(is(x,"ce"))
    xx <- mktce(x,ret)
    xx-apply(xx,1,mean,na.rm=TRUE)
}
#mktce-----------------------------------
##' Extractor for market component of return (factor 1)
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix of attributed return
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export mktce
`mktce` <- function(x,ret){
    stopifnot(is(x,"ce") && buialigned(x,ret))
    mz(scoce(x,ret)[,1,drop=FALSE]%*%t(ldgce(x)[,1,drop=FALSE]))
}
#sysce-----------------------------------
##' Extractor for systematic component of return (factors 2:k)
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix of attributed return
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export sysce
`sysce` <- function(x,ret){
    stopifnot(is(x,"ce") && buialigned(x,ret))
    mz(scoce(x,ret)[,-1,drop=FALSE]%*%t(ldgce(x)[,-1,drop=FALSE]))
}
#msrtce-----------------------------------
##' Extractor for components Market, Systematic, Residual, Total
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix of attributed return
##' @section Details: the output has 4x columns of the input, grouped by component
##' @author Giles Heywood 
##' @family extractors
##' @export msrtce
`msrtce` <- function(x,ret) {
    stopifnot(is(x,"ce"))
    if(any(is.na(ret))) coredata(ret)[which(is.na(ret))] <- coredata(msce(ce,natox(ret,0)))[which(is.na(ret))]
    sco <- scoce(x,ret)
    cbind(
        mz(sco[,1,drop=FALSE]%*%t(ldgce(x)[,1,drop=FALSE])),
        mz(sco[,-1,drop=FALSE]%*%t(ldgce(x)[,-1,drop=FALSE])),
        ret-mz(sco%*%t(ldgce(x))),
        ret
        )
}
#msce-----------------------------------
##' Extractor for components Market+Systematic
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return matrix of attributed return
##' @section Details:
##' @author Giles Heywood 
##' @family extractors
##' @export msce
`msce` <- function(x,ret){
    stopifnot(is(x,"ce") && buialigned(x,ret))
    mz(scoce(x,ret)%*%t(ldgce(x)))
}
#stce-----------------------------------
##' Market,Systematic,Residual,Total returns, divided by estimated vol
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return named list of attributed and normalised return
##' @section Details: components / vol estimate, in-sample
##' @author Giles Heywood 
##' @family extractors
##' @export stce
`stce` <- function(x,ret){
    stopifnot(is(x,"ce"))
    if(!vz(ret)) ret <- mz(ret)
    varmsrt <- prdce(x,po=poce(x),security=TRUE)[,c("vam","vas","var","vat")]
    varfix <- sweep(x=varmsrt,MARGIN=2,STATS=apply(varmsrt,2,median)*1.e-3+1.e-7,FUN="+") #stabilise denominator
    M <- sweep(x=mktce(x,ret),MARGIN=2,STATS=sqrt(varfix[,"vam"]),FUN="/")
    S <- sweep(x=sysce(x,ret),MARGIN=2,STATS=sqrt(varfix[,"vas"]),FUN="/")
    R <- sweep(x=ret-msce(x,ret),MARGIN=2,STATS=sqrt(varfix[,"var"]),FUN="/")
    T <- sweep(x=ret,MARGIN=2,STATS=sqrt(varfix[,"vat"]),FUN="/")
    list(remnor=mz(M),resnor=mz(S),rernor=mz(R),retnor=mz(T))
}
#face-----------------------------------
##' Returns decomposed by all factors 1:k, + residual, + quadratic
##' 
##' @param x object of class ce
##' @param ret returns matrix
##' @return named list of attributed return: 1:k, ret, rer, qua
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export face
`face` <- function(x,ret){
    stopifnot(is(x,"ce"))
    if(!vz(ret)) ret <- mz(ret)
    sco <- scoce(x,ret)
    ldg <- t(ldgce(x))
    qua <- quace(x,ret) 
    n <- nrow(ldg)
    res <- vector("list",n+3)
    res[[n+1]] <- ret
    res[[n+2]] <- ret
    res[[n+3]] <- qua*NA    
    for(i in 1:n) {
        res[[i]] <- mz(sco[,i,drop=FALSE]%*%ldg[i,,drop=FALSE])
        res[[n+2]] <- res[[n+2]] - res[[i]]
    }
    names(res) <- c(psz("re",1:n),"ret","rer","qua")
    res
}
#fmpce-----------------------------------
##' Factor mimicking portfolio weights (but still in reduced dimension ie non-NA columns only)
##' 
##' @param x object of class ce
##' @param type fmp for factor portfolio weights
##' @return matrix
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export fmpce
`fmpce` <- function(x,type=c("fmp","hpl")){
    type <- match.arg(type)
    stopifnot(is(x,"ce"))
    if(type=="fmp") {
        x$fmp/as.numeric(x$sdev)
    } else {
        x$hpl/as.numeric(x$sdev)
    }
}
#metce-----------------------------------
##' Extractor for method
##' 
##' @param x object of class ce
##' @return method
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export metce
#metce - method, (colnames=F, rownames=T)
`metce` <- function(x){
    stopifnot(is(x,"ce"))
    x$method
}
#scoce-----------------------------------
##' Scores
##' 
##' @param x object of class ce
##' @param ret returns
##' @param type fmp for scores
##' @return scores
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export scoce
`scoce` <- function(x,ret,type=c("fmp","hpl")){
    type <- match.arg(type)
    stopifnot(is(x,"ce") && buialigned(x,ret))
    stopifnot(all(fulce(x)%in%colnames(ret)))
    fmp <- fmpce(x,type=type)
    mz(coredata(ret)[,fulce(x),drop=FALSE]%*%coredata(fmp)[fulce(x),,drop=FALSE])
}
#resce-----------------------------------
##' Residual return
##' 
##' @param x object of class ce
##' @param ret returns
##' @return residuals
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export resce
`resce` <- function(x,ret){
    stopifnot(is(x,"ce"))
    mz(ret-scoce(x,ret)%*%t(ldgce(x)))
}
#spvce-----------------------------------
##' Extractor for pecific vol
##' 
##' @param x object of class ce
##' @return specific vol
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export spvce
`spvce` <- function(x){
    stopifnot(is(x,"ce"))
    x$sdev*sqrt(x$uniqueness)
}
#sdvce-----------------------------------
##' Extractor for standard deviation
##' 
##' @param x object of class ce
##' @return standard deviation (vol)
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export sdvce
`sdvce` <- function(x){
    stopifnot(is(x,"ce"))
    x$sdev
}
#unqce-----------------------------------
##' Extractor for uniqueness
##' 
##' @param x object of class ce
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export unqce
`unqce` <- function(x){
    stopifnot(is(x,"ce"))
    x$uniqueness
}
#fulce-----------------------------------
##' Extractor for identifiers of columns with no NA (full data)
##' 
##' @param x object of class ce
##' @return identifiers
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export fulce
`fulce` <- function(x){
    stopifnot(is(x,"ce"))
    rownames(x$full)[x$full]
}
#buice-----------------------------------
##' Extractor for identifiers (all)
##' 
##' @param x object of class ce
##' @return identifiers : colnames from original data
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export buice
`buice` <- function(x){
    stopifnot(is(x,"ce"))
    rownames(x$full) #sort(rownames(x$full))
}
#prdce-----------------------------------
##' Portfolio risk decomposition
##' 
##' @param x object of class ce
##' @param po column matrix of portfolio weights with rownames=identifiers
##' @param scaletovol flag : scale variance contributions by 1/sqrt(total variance), so 'units of variance are rescaled to add up to portfolio vol'
##' @param security flags use of variance only
##' @return identifiers : matrix, rows are securities, columns are M=Market, S=Systematic, R=Residual, TOT=total, F1=factor1 vol
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export prdce
`prdce` <- function(
                x=getce(),
                po=poce(x),
                scaletovol=FALSE,
                security=FALSE) {#if security=TRUE, only diagonal elements are used, but weights po are still used
    stopifnot(is(x,"ce"))
    aggregvar <- if(security) diag else rowSums
    n <- length(po)
    stopifnot(all(rownames(po) %in% buice(x)))
    bui <- rownames(po)
    wldgs <- sweep(x=x$loadings[bui,,drop=FALSE], MARGIN=1, STATS=x$sdev[bui,,drop=FALSE]*as.numeric(po), FUN="*")
    M <- aggregvar(tcrossprod(wldgs[,1,drop=FALSE]))
    S <- aggregvar(tcrossprod(wldgs[,-1,drop=FALSE]))
    R <- (x$sdev[bui,,drop=FALSE]**2)*x$uniqueness[bui,,drop=FALSE]*po**2
    TOT <- M+S+R
    if(scaletovol & !security) {
        pvol <- sqrt(sum(TOT))
        M <- M/pvol
        S <- S/pvol
        R <- R/pvol
        TOT <- TOT/pvol
    }
    F1 <- wldgs[,1,drop=FALSE]
    matrix(c(M,S,R,TOT,F1),n,5,dimnames=list(bui,c("vam","vas","var","vat","f1")))
}
#genfrdce-----------------------------------
##' Factor risk decomposition by security
##' 
##' @param x object of class ce
##' @return matrix, rows=security, columns = (total, residual, factors 1:k)
##' @section Details: 
##' @author Giles Heywood 
##' @family extractors
##' @export genfrdce
`genfrdce` <- function(x=getce()) {
    stopifnot(is(x,"ce"))
    va <- sweep(x=cbind(1,x$uniqueness,x$loadings**2), MARGIN=1, STATS=x$sdev**2, FUN="*")
    colnames(va) <- c("vat","var",psz("va",1:ncol(x$loadings)))
    rownames(va) <- rownames(x$loadings)
    va
}
`buialigned` <- function(x,ret) { all(colnames(ret)==rownames(x$loadings)) }

library(nlme)
library(mvtnorm)
require(MASS)

tr = function(M)  sum(diag(M))

vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

cor.fn = function(phi, ga=NULL, dim, type = c("UNC", "DEC", "CAR1", "ARp"), Ti=NULL)
{
  if(type[1] == "UNC") return (diag(dim))
  if(type[1] == "DEC") return(phi^((abs(outer(Ti, Ti, '-')))^ga))
  if(type[1] == "CAR1") return(phi^(abs(outer(Ti, Ti, '-'))))
  if(type[1] == "ARp") return(arp.Ci(phi, P = length(phi), dim = dim))
  Cor = NULL
  for(i in 1:dim) Cor = rbind(Cor, tmp[1:dim + (dim - i)])
  return(Cor)
}

mse = function(v1, v2) mean((v1 - v2)^2)

# dot DEC:
DEC.dot.phi = function(phi, ga, Ti)
{
  tem = abs(outer(Ti, Ti,'-'))^ga
  dot.CAR = tem*(phi^(tem-1))
  return(dot.CAR)
}
DEC.dot.ga = function(phi, ga, Ti)
{
  tem = abs(outer(Ti, Ti,'-'))
  dot.ga.CAR = log(phi^tem)*cor.fn(phi=phi, ga=ga, type="DEC", Ti=Ti)
  return(dot.ga.CAR)
}

# Autocorrelation matrix:
arp.Ci= function(phi, P = length(phi), dim = P + 1, sigma.wn = 1)
{
  if(dim <= 0){return(cat('Waring: dim must be positive integer. \n'))}
  n = dim
  if(n <= P){n = P+1}
  posi.mt = abs(1:P-outer(rep(1,P), 0:P))
  coef.mt = outer(rep(1,P), c(-1,phi))
  coef.rho = matrix(0, P, P)
  for(i in 1:P){
    for(j in 1:P){
      coef.rho[i,j] = sum(coef.mt[i,][posi.mt[i,]==j])
    }
  }
  rho =  - solve(coef.rho) %*% phi
  gamma0 = sigma.wn^2/(1-sum(rho*phi))
  rho = c(rho, rep(NA, (n-P)))
  for(i in (P+1):n){
    rho[i] = sum(rev(phi)*rho[(i-P):(i-1)])
  }
  rho = c(rev(rho), 1, rho)
  corr.mt = matrix(NA, n, n)
  for(i in 1:n){
    corr.mt[i,] = rho[(n+2-i):(2*n-i+1)]
  }
  return(gamma0*corr.mt[1:dim,1:dim])
}

HL = function(phi, dim=ni)
{
  P = length(phi)
  tmp.para = c(1, - phi)
  HL.para = matrix(0, dim, P + dim)
  for(i in 1: dim)
  {
    HL.para[i,(i + P):i] = tmp.para
  }
  return(HL.para)
}

d.inv.Ci.fn = function(phi, dim=ni)
{
  P = length(phi)
  HL.phi = HL(phi, dim)
  Hp = HL.phi[,1:P]
  Lp = HL.phi[,-(1:P)]

  dS = as.list(numeric(P))
  posi = outer(1:dim, 1:(P+dim), function(i,j,P) i + P - j, P = P)
  for(k in 1:P)
  {
    d.HL.k = matrix(0, dim, P + dim)
    d.HL.k[posi == k] = - 1
    dH = d.HL.k[,1:P]
    dL = d.HL.k[,-(1:P)]
    dS[[k]] = t(dL)%*%Lp+t(Lp)%*%dL-dH%*%t(Hp)-Hp%*%t(dH)
  }
  return(dS)
}

arp.Ci.dot=function(phi, dim, k) -arp.Ci(phi, dim=dim)%*%d.inv.Ci.fn(phi, dim=dim)[[k]]%*%arp.Ci(phi, dim=dim)

pacf.to.phi=function(pacf)
{
  P=length(pacf)
    if(P==1) phi=pacf
      if(P>1){
        Phi=matrix(diag(pacf), P, P)
        for (i in 2: P)
           for (j in 1:(i-1))
             Phi[i, j]=Phi[i-1, j]-Phi[i,i]*Phi[i-1, i-j]
        phi=Phi[P,]
      }
  return(phi)
}

# classError function in package: mclust
classError = function(classification, class)
{
    q <- function(map, len, x) {
        x <- as.character(x)
        map <- lapply(map, as.character)
        y <- sapply(map, function(x) x[1])
        best <- y != x
        if (all(len) == 1)
            return(best)
        errmin <- sum(as.numeric(best))
        z <- sapply(map, function(x) x[length(x)])
        mask <- len != 1
        counter <- rep(0, length(len))
        k <- sum(as.numeric(mask))
        j <- 0
        while (y != z) {
            i <- k - j
            m <- mask[i]
            counter[m] <- (counter[m]%%len[m]) + 1
            y[x == names(map)[m]] <- map[[m]][counter[m]]
            temp <- y != x
            err <- sum(as.numeric(temp))
            if (err < errmin) {
                errmin <- err
                best <- temp
            }
            j <- (j + 1)%%k
        }
        best
    }
    if (any(isNA <- is.na(classification))) {
        classification <- as.character(classification)
        nachar <- paste(unique(classification[!isNA]), collapse = "")
        classification[isNA] <- nachar
    }
    MAP <- mapClass(classification, class)
    len <- sapply(MAP[[1]], length)
    if (all(len) == 1) {
        CtoT <- unlist(MAP[[1]])
        I <- match(as.character(classification), names(CtoT),
            nomatch = 0)
        one <- CtoT[I] != class
    }
    else {
        one <- q(MAP[[1]], len, class)
    }
    len <- sapply(MAP[[2]], length)
    if (all(len) == 1) {
        TtoC <- unlist(MAP[[2]])
        I <- match(as.character(class), names(TtoC), nomatch = 0)
        two <- TtoC[I] != classification
    }
    else {
        two <- q(MAP[[2]], len, classification)
    }
    err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
        as.vector(one)
    else as.vector(two)
    bad <- seq(along = classification)[err]
    list(misclassified = bad, errorRate = length(bad)/length(class))
}

mapClass = function(a, b)
{
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
        warning("unequal lengths")
        return(x)
    }
    if (is.factor(a) & is.factor(b) & nlevels(a) == nlevels(b)) {
        aTOb <- as.list(levels(b))
        names(aTOb) <- levels(a)
        bTOa <- as.list(levels(a))
        names(bTOa) <- levels(b)
        out <- list(aTOb = aTOb, bTOa = bTOa)
        return(out)
    }
    if (is.character(a) & is.character(b) & length(unique(a)) ==
        length(unique(b))) {
        aTOb <- as.list(unique(b))
        names(aTOb) <- unique(a)
        bTOa <- as.list(unique(a))
        names(bTOa) <- unique(b)
        out <- list(aTOb = aTOb, bTOa = bTOa)
        return(out)
    }
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
        I <- match(Max[i], Tab[i, ], nomatch = 0)
        aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b))
        aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
        J <- match(Max[j], Tab[, j])
        bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a))
        bTOa <- lapply(bTOa, as.numeric)
    out <- list(aTOb = aTOb, bTOa = bTOa)
    return(out)
}

# adjustedRandIndex function in package: mclust
adjustedRandIndex = function(x, y)
{
    x <- as.vector(x)
    y <- as.vector(y)
    if (length(x) != length(y))
        stop("arguments must be vectors of the same length")
    tab <- table(x, y)
    if (all(dim(tab) == c(1, 1)))
        return(1)
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c
    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    return(ARI)
}

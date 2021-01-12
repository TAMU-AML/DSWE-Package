# MIT License
# 
# Copyright (c) 2020 Nitesh Kumar, Abhinav Prakash, and Yu Ding
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

get.dpill = function (cov, y)
{
  bw <- KernSmooth::dpill(cov, y)
  if (is.nan(bw)) {
    par <- 0.06
    while (is.nan(bw)) {
      bw <- KernSmooth::dpill(cov, y, proptrun = par)
      par <- par + 0.01
    }
  }
  return(bw)
}

cut.pts = function (cov, circ = TRUE, .n.reg)
{
  hobj <- graphics::hist(cov, breaks = "FD", plot = FALSE)
  brks <- hobj$breaks
  cnts <- hobj$counts
  id.cand <- which(cnts == 0)
  if (length(id.cand) > 0) {
    loc.cut <- c()
    if (circ & brks[1] <= 0 & brks[length(brks)] >= 360) {
      if (length(id.cand) == 1)
        loc.cut <- (brks[id.cand] + brks[id.cand + 1])/2
      else {
        id.diff <- diff(id.cand)
        n.cons <- rle(id.diff)
        vals <- n.cons$values
        lens <- n.cons$lengths
        if (any(vals == 1)) {
          id.pos <- which(lens == max(lens[which(vals ==
                                                   1)]))
          if (id.pos == 1)
            id.rng <- c(id.cand[1], id.cand[sum(lens[1:(id.pos)]) +
                                              1])
          else id.rng <- id.cand[c(sum(lens[1:(id.pos -
                                                 1)]), sum(lens[1:(id.pos)])) + 1]
        }
        else id.rng <- rep(id.cand[1], 2)
        loc.cut <- (brks[id.rng[1]] + brks[(id.rng[2] +
                                              1)])/2
      }
      cov[which(cov < loc.cut)] <- cov[which(cov < loc.cut)] +
        360
    }
    gmm <- mixtools::normalmixEM(cov, k = .n.reg, verb = FALSE)
    if (length(which(gmm$posterior[, 1] > gmm$posterior[,
                                                        2])) > 0 & length(which(gmm$posterior[, 1] < gmm$posterior[,
                                                                                                                   2])) > 0) {
      if (gmm$mu[1] < gmm$mu[2])
        id.clust <- which(gmm$posterior[, 1] > gmm$posterior[,
                                                             2])
      else id.clust <- which(gmm$posterior[, 1] < gmm$posterior[,
                                                                2])
      bnd1 <- max(cov[id.clust])
      bnd2 <- min(cov[-id.clust])
      loc.cut <- c(loc.cut, ((bnd1 + bnd2)/2))
    }
    if (circ & any(loc.cut > 360))
      loc.cut[which(loc.cut > 360)] <- loc.cut[which(loc.cut >
                                                       360)] - 360
    loc.cut <- sort(loc.cut)
    loc.cut <- c(brks[1], loc.cut, brks[length(brks)])
    return(loc.cut)
  }
  else {
    return(NA)
  }
}

find.bw = function (y.tr, X.tr, X.ts, bw)
{
  X.tr <- as.matrix(X.tr)
  X.ts <- as.matrix(X.ts)
  y.tr <- as.matrix(y.tr)
  p <- ncol(X.tr)
  n.tr <- nrow(X.tr)
  n.ts <- nrow(X.ts)
  id.adp <- bw$id.adp
  id.fix <- 1:p
  if (!is.na(id.adp)){
    id.fix <- id.fix[-which(id.fix == id.adp)]
  }
  h.adp <- bw$bw.adp
  h.fix <- bw$bw.fix
  cutpt <- bw$cutpt
  bins <- c()
  if (!is.na(id.adp)) {
    bins <- as.data.frame(matrix(NA, nrow = n.ts, ncol = length(id.adp)))
    for (k in 1:length(id.adp)) {
      bins[, k] <- .bincode(X.ts[, id.adp[k]], cutpt[[k]],
                            include.lowest = TRUE)

    }
  }

  diff <- get.diff(X.tr, X.ts)
  h <- rep(NA, p)
  h[id.fix] <- h.fix
  for (q in 1:length(id.adp)) h[id.adp[q]] <- h.adp[[q]][bins[, q]]

  return(h)

}

get.diff = function (X.tr, x.TS)
{
  X.tr <- as.matrix(X.tr)
  n.TR <- dim(X.tr)[1]
  q.TR <- dim(X.tr)[2]
  x.TS <- matrix(x.TS, 1, q.TR)
  oneV <- matrix(1, n.TR, 1)
  diff <- X.tr - (oneV %*% x.TS)
  return(diff)
}

bw.gap = function (y, x, id.dir = NA, id.adp = id.dir) {

  if (is.na(id.adp)) {
    bw.fix <- sapply(1:ncol(x), function(p) get.dpill(x[,
                                                        p], y))
    list(bw.fix = bw.fix, bw.adp = NA, id.adp = NA, cutpt = NA)
  }
  else {
    cutpt <- lapply(id.adp, function(p) cut.pts(cov = x[,
                                                        p], circ = (p %in% id.dir), .n.reg = 2))
    if (all(sapply(cutpt, function(pts) all(is.na(pts))))) {
      bw.fix <- sapply(c(1:ncol(x)), function(p) get.dpill(x[,
                                                             p], y))
      list(bw.fix = bw.fix, bw.adp = NA, id.adp = NA,
           cutpt = NA)
    }
    else {
      bins <- sapply(1:length(id.adp), function(p) .bincode(x[,
                                                              id.adp[p]], breaks = cutpt[[p]], include.lowest = TRUE))
      bw.fix <- sapply(c(1:ncol(x))[-id.adp], function(p) get.dpill(x[,
                                                                      p], y))
      bw.adp <- lapply(1:length(id.adp), function(p) {
        n.bin <- length(unique(bins[, p]))
        if (id.adp[p] %in% id.dir & cutpt[[p]][1] <=
            0 & cutpt[[p]][length(cutpt[[p]])] >= 360) {
          id.adj <- which(bins[, p] == n.bin)
          bins[id.adj, p] <- 1
          n.bin <- n.bin - 1
        }
        bw <- sapply(1:n.bin, function(b) {
          id.bin <- which(bins[, p] == b)
          cov.sub <- x[id.bin, id.adp[p]]
          if (id.adp[p] %in% id.dir & cutpt[[p]][1] <=
              0 & cutpt[[p]][length(cutpt[[p]])] >= 360)
            cov.sub[which(cov.sub < cutpt[[p]][2])] <- cov.sub[which(cov.sub <
                                                                       cutpt[[p]][2])] + 360
          KernSmooth::dpill(cov.sub, y[id.bin])
        })
      })
      list(bw.fix = bw.fix, bw.adp = bw.adp, id.adp = id.adp,
           cutpt = cutpt)
    }
  }
}


computePredGap = function(trainX, trainY, testX, bandwidth, nMultiCov, fixedCov, cirCov ){
  if(!is.na(cirCov)){
    for (i in cirCov) {
      trainX[,i] = trainX[,i]*pi/180
      testX[,i] = testX[,i]*pi/180
    }
  }

  weights = calculateWeights(trainX,testX,bandwidth,nMultiCov,fixedCov,cirCov)
  pred = (rowSums(weights)%*%trainY)/ncol(weights)

  return(pred)
}

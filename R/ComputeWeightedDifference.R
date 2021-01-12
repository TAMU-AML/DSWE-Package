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

#' @title Percentage weighted difference between power curves
#' @description Computes percentage weighted difference between power curves based on user provided weights instead of the weights computed from the data. Please see \code{details} for more information.  
#'
#' @param muDiff a vector of pointwise difference between two power curves on a testset as obtained from \code{ComparePCurve()} or \code{funGP()} function.
#' @param weights a vector of user specified weights for each element of \code{muDiff}. It can be based on any probability distribution of user's choice. The weights must sum to 1.
#' @param base a vector of predictions from a power curve; to be used as the denominator in computing the percentage difference. It can be either \code{mu1} or \code{mu2} as obtained from \code{ComparePCurve()} or \code{funGP()} function. 
#' @param statDiff a boolean specifying whether to compute the statistical significant difference or not. Default is set to \code{FALSE}, i.e. statistical significant difference is not computed. If set to \code{TRUE}, \code{confBand} must be provided.
#' @param confBand a vector of pointwise confidence band for all the points in the testset as obtained from \code{ComparePCurve()} or \code{funGP()} function, named as \code{band}. Should only be provided when \code{statDiff} is set to \code{TRUE}. Default value is \code{NULL}.
#'
#' @return a numeric percentage weighted difference or statistical significant percetage weighted difference based on whether statDiff is set to \code{FALSE} or \code{TRUE}. 
#' @examples 
#' 
#' ws_test = as.matrix(seq(4.5,8.5,length.out = 10))
#' 

#' userweights = dweibull(ws_test, shape = 2.25, scale = 6.5) 
#' userweights = userweights/sum(userweights) 


#' data1 = data1[1:100, ]
#' data2 = data2[1:100, ]
#' datalist = list(data1, data2)
#' xCol = 2
#' xCol.circ = NULL
#' yCol = 7
#' testCol = 2

#' output = ComparePCurve(data = datalist, xCol = xCol, yCol = yCol, 
#' testCol = testCol, testSet = ws_test) 

#' weightedDiff = ComputeWeightedDifference(output$muDiff, userweights, output$mu1)

#' weightedStatDiff = ComputeWeightedDifference(output$muDiff, userweights, output$mu1, 
#' statDiff = TRUE, confBand = output$band)
#' 
#' @details The function is a modification to the percentage weighted difference defined in Ding et. al. (2020). It computes a weighted difference between power curves on a testset, where the weights have to be provided by the user based on any probability distribution of their choice rather than the weights being computed from the data. The weights must sum to 1 to be valid.
#' @references For details, see Ding et. al. (2020) available on \code{arxiv} at this \href{https://arxiv.org/abs/2005.08652}{link}.
#' @export

ComputeWeightedDifference = function(muDiff, weights, base, statDiff = FALSE, confBand = NULL){
  if (!is.numeric(muDiff)){
    stop('muDiff must be a numeric vector' )
  }
  if (!is.numeric(weights)){
    stop('weights must be a numeric vector')
  }
  if (!is.numeric(base)){
    stop('base must be a numeric vector')
  }
  if (length(muDiff) != length(weights) || length(muDiff) != length(base)){
    stop('length of muDiff, weights and base must be the same')
  }
  if(abs(sum(weights) - 1) > .Machine$double.eps){
    stop('the weights must sum to 1 for being a valid weights')
  }
  if (!is.logical(statDiff)){
    stop('statDiff must be either TRUE or FALSE')
  } else if (statDiff){
    if (is.null(confBand)){
      stop('confBand must be provided when statDiff is set to TRUE')
    }
    if (!is.numeric(confBand)){
      stop('confBand must be a numeric vector when statDiff is set to TRUE')
    }
    if (length(muDiff) != length(confBand)){
      stop('length of confBand must be the same as muDiff')
    }
  }
  weightedDiff = computeWeightedDiffExtern(muDiff, weights, base, statDiff, confBand)
  return(as.numeric(weightedDiff))
}

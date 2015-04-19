rm(list=ls())
setInternet2(TRUE)
con = gzcon(url('http://www.systematicportfolio.com/sit.gz', 'rb'))
source(con)
close(con)
load.packages("TTR,PerformanceAnalytics,quantmod,lattice")
source("C:/Users/rf6994/Documents/R/FAA03/FAAreturns.R")

require(quantmod)
require(PerformanceAnalytics)
  
  mutualFunds <- c("CSD", "EDV", "VNQ", "BND", "MUB", "IGOV", "VWO", "ILF", "EWJ","VFISX") 
  #mutualFunds <- c("VTSMX", #Vanguard Total Stock Market Index
  #                 "FDIVX", #Fidelity Diversified International Fund
  #                 "VEIEX", #Vanguard Emerging Markets Stock Index Fund
  #                 "VFISX", #Vanguard Short-Term Treasury Fund
  #                 "VBMFX", #Vanguard Total Bond Market Index Fund
  #                 "QRAAX", #Oppenheimer Commodity Strategy Total Return
  #                 "VGSIX" #Vanguard REIT Index Fund
  #)
  
  #mid 1997 to end of 2012
  #getSymbols(mutualFunds, from="1997-06-30", to="2015-03-20")
  getSymbols(mutualFunds, from="2009-06-30", to="2015-04-17")
  tmp <- list()
  for(fund in mutualFunds) {
    tmp[[fund]] <- Ad(get(fund))
  }
  
  #always use a list when intending to cbind/rbind large quantities of objects
  adPrices <- do.call(cbind, args = tmp)
  colnames(adPrices) <- gsub(".Adjusted", "", colnames(adPrices))

':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}

  FAA <- function(prices, monthLookback = 4,
                  weightMom = 1, weightVol = .5, weightCor = .5,
                  riskFreeName = NULL, bestN = 3,
                  stepCorRank = FALSE, stepStartMethod = c("best", "default"),
                  geometric = TRUE, ...) {
    stepStartMethod <- stepStartMethod[1]
    if(is.null(riskFreeName)) {
      prices$zeroes <- 0
      riskFreeName <- "zeroes"
      warning("No risk-free security specified. Recommended to use one of: quandClean('CHRIS/CME_US'), SHY, or VFISX.
              Using vector of zeroes instead.")
    }
    returns <- Return.calculate(prices)
    monthlyEps <- endpoints(prices, on = "months")
    riskFreeCol <- grep(riskFreeName, colnames(prices))
    tmp <- list()
    dates <- list()
    
    for(i in 2:(length(monthlyEps) - monthLookback)) {
      #subset data
      priceData <- prices[monthlyEps[i]:monthlyEps[i+monthLookback],]
      returnsData <- returns[monthlyEps[i]:monthlyEps[i+monthLookback],]
      
      #perform computations
      momentum <- data.frame(t(t(priceData[nrow(priceData),])/t(priceData[1,]) - 1))
      momentum <- momentum[,!is.na(momentum)]
      #momentum[is.na(momentum)] <- -1 #set any NA momentum to negative 1 to keep R from crashing
      priceData <- priceData[,names(momentum)]
      returnsData <- returnsData[,names(momentum)]
      
      momRank <- rank(momentum)
      vols <- data.frame(StdDev(returnsData))
      volRank <- rank(-vols)
      cors <- cor(returnsData, use = "complete.obs")
      if (stepCorRank) {
        if(stepStartMethod=="best") {
          compositeMomVolRanks <- weightMom*momRank + weightVol*volRank
          maxRank <- compositeMomVolRanks[compositeMomVolRanks==max(compositeMomVolRanks)]
          corRank <- stepwiseCorRank(corMatrix=cors, startNames = names(maxRank),
                                     bestHighestRank = TRUE, ...)
          
        } else {
          corRank <- stepwiseCorRank(corMatrix=cors, bestHighestRank=TRUE, ...)
        }
      } else {
        corRank <- rank(-rowSums(cors))
      }
      
      totalRank <- rank(weightMom*momRank + weightVol*volRank + weightCor*corRank)
      
      upper <- length(names(returnsData))
      lower <- max(upper-bestN+1, 1)
      topNvals <- sort(totalRank, partial=seq(from=upper, to=lower))[c(upper:lower)]
      
      #compute weights
      longs <- totalRank %in% topNvals #invest in ranks length - bestN or higher (in R, rank 1 is lowest)
      longs[momentum < 0] <- 0 #in previous algorithm, removed momentums < 0, this time, we zero them out at the end.
      longs <- longs/sum(longs) #equal weight all candidates
      longs[longs > 1/bestN] <- 1/bestN #in the event that we have fewer than top N invested into, lower weights to 1/top N
      names(longs) <- names(totalRank)
      
      
      #append removed names (those with momentum < 0)
      removedZeroes <- rep(0, ncol(returns)-length(longs))
      names(removedZeroes) <- names(returns)[!names(returns) %in% names(longs)]
      longs <- c(longs, removedZeroes)
      
      #reorder to be in the same column order as original returns/prices
      longs <- data.frame(t(longs))
      longs <- longs[, names(returns)]
      
      #append lists
      tmp[[i]] <- longs
      dates[[i]] <- index(returnsData)[nrow(returnsData)]
    }
    
    weights <- do.call(rbind, tmp)
    dates <- do.call(c, dates)
    weights <- xts(weights, order.by=as.Date(dates))
    weights[, riskFreeCol] <- weights[, riskFreeCol] + 1-rowSums(weights)
    strategyReturns <- Return.rebalancing(R = returns, weights = weights, geometric = geometric)
    colnames(strategyReturns) <- paste(monthLookback, weightMom, weightVol, weightCor, sep="_")
    return(strategyReturns)
    }
  
  c(ra_wts, ra_ret) := FAAreturns(adPrices)
  c(n4_wts, n4_ret) := FAAreturns(adPrices, bestN=4)
  c(n3_wts, n3_ret) := FAAreturns(adPrices, bestN=3)
  c(volcor_wts, volcor_ret) := FAAreturns(adPrices, weightVol = 1, weightCor = 1)
  c(minrisk_wts, minrisk_ret) := FAAreturns(adPrices, weightMom = 0, weightVol=1, weightCor=1)
  c(puremom_wts, puremom_ret) := FAAreturns(adPrices, weightMom=1, weightVol=0, weightCor=0)
  c(maxdecor_wts, maxdecor_ret) := FAAreturns(adPrices, weightMom=0, weightVol=0, weightCor=1)
  c(momdecor_wts, momdecor_ret) := FAAreturns(adPrices, weightMom=1, weightVol=0, weightCor=1)
  
  all <- cbind(ra_ret, n4_ret, n3_ret, volcor_ret, minrisk_ret, puremom_ret, maxdecor_ret, momdecor_ret)
  colnames(all) <- c("Replica Attempt", "N4", "N3", "N3vol1cor1", "minRisk", "pureMomentum", "maxDecor", "momDecor")
  charts.PerformanceSummary(all, colorset=c("black", "red", "blue", "green", "brown", "darkgrey", "purple", "orange"))
  
  stats <- data.frame(t(rbind(Return.annualized(all)*100,
                              maxDrawdown(all)*100,
                              SharpeRatio.annualized(all))))
  stats$Return_To_Drawdown <- stats[,1]/stats[,2]
  
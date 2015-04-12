require(quantmod)
require(PerformanceAnalytics)

stepwiseCorRank <- function(corMatrix, startNames=NULL, stepSize=1, bestHighestRank=FALSE) {
  #edge cases
  if(dim(corMatrix)[1] == 1) {
    return(corMatrix)
  } else if (dim(corMatrix)[1] == 2) {
    ranks <- c(1.5, 1.5)
    names(ranks) <- colnames(corMatrix)
    return(ranks)
  }
  
  if(is.null(startNames)) {
    corSums <- rowSums(corMatrix)
    corRanks <- rank(corSums)
    startNames <- names(corRanks)[corRanks <= stepSize]
  }
  nameList <- list()
  nameList[[1]] <- startNames
  rankList <- list()
  rankCount <- 1
  rankList[[1]] <- rep(rankCount, length(startNames))
  rankedNames <- do.call(c, nameList)
  
  while(length(rankedNames) < nrow(corMatrix)) {
    rankCount <- rankCount+1
    subsetCor <- corMatrix[, rankedNames]
    if(class(subsetCor) != "numeric") {
      subsetCor <- subsetCor[!rownames(corMatrix) %in% rankedNames,]
      if(class(subsetCor) != "numeric") {
        corSums <- rowSums(subsetCor)
        corSumRank <- rank(corSums)
        lowestCorNames <- names(corSumRank)[corSumRank <= stepSize]
        nameList[[rankCount]] <- lowestCorNames
        rankList[[rankCount]] <- rep(rankCount, min(stepSize, length(lowestCorNames)))
      } else { #1 name remaining
        nameList[[rankCount]] <- rownames(corMatrix)[!rownames(corMatrix) %in% names(subsetCor)]
        rankList[[rankCount]] <- rankCount
      }
    } else {  #first iteration, subset on first name
      subsetCorRank <- rank(subsetCor)
      lowestCorNames <- names(subsetCorRank)[subsetCorRank <= stepSize]
      nameList[[rankCount]] <- lowestCorNames
      rankList[[rankCount]] <- rep(rankCount, min(stepSize, length(lowestCorNames)))
    }   
    rankedNames <- do.call(c, nameList)
  }
  
  ranks <- do.call(c, rankList)
  names(ranks) <- rankedNames
  if(bestHighestRank) {
    ranks <- 1+length(ranks)-ranks
  }
  ranks <- ranks[colnames(corMatrix)] #return to original order
  return(ranks)
}


FAAreturns <- function(prices, monthLookback = 4,
                       weightMom=1, weightVol=.5, weightCor=.5,
                       riskFreeName="VFISX", bestN=3,
                       stepCorRank = FALSE, stepStartMethod=c("best", "default")) {
  stepStartMethod <- stepStartMethod[1]
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
    priceData <- priceData[, momentum > 0] #remove securities with momentum < 0
    returnsData <- returnsData[, momentum > 0]
    momentum <- momentum[momentum > 0]
    names(momentum) <- colnames(returnsData)
    vol <- as.numeric(-sd.annualized(returnsData))
    
    if(length(momentum) > 1) {
      
      #perform ranking
      if(!stepCorRank) {
        sumCors <- -colSums(cor(returnsData, use="complete.obs"))
        stats <- data.frame(cbind(momentum, vol, sumCors))
        ranks <- data.frame(apply(stats, 2, rank))
        weightRankSum <- weightMom*ranks$momentum + weightVol*ranks$vol + weightCor*ranks$sumCors
        names(weightRankSum) <- rownames(ranks)
      } else {
        corMatrix <- cor(returnsData, use="complete.obs")
        momRank <- rank(momentum)
        volRank <- rank(vol)
        compositeMomVolRanks <- weightMom*momRank + weightVol*volRank
        maxRank <- compositeMomVolRanks[compositeMomVolRanks==max(compositeMomVolRanks)]
        if(stepStartMethod=="default") {
          stepCorRanks <- stepwiseCorRank(corMatrix=corMatrix, startNames = NULL,
                                          stepSize = 1, bestHighestRank = TRUE)
        } else {
          stepCorRanks <- stepwiseCorRank(corMatrix=corMatrix, startNames = names(maxRank),
                                          stepSize = 1, bestHighestRank = TRUE)
        }
        weightRankSum <- weightMom*momRank + weightVol*volRank + weightCor*stepCorRanks
      }
      
      totalRank <- rank(weightRankSum)
      
      #find top N values, from http://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
      #thanks to Dr. Rob J. Hyndman
      upper <- length(names(returnsData))
      lower <- max(upper-bestN+1, 1)
      topNvals <- sort(totalRank, partial=seq(from=upper, to=lower))[c(upper:lower)]
      
      #compute weights
      longs <- totalRank %in% topNvals #invest in ranks length - bestN or higher (in R, rank 1 is lowest)
      longs <- longs/sum(longs) #equal weight all candidates
      longs[longs > 1/bestN] <- 1/bestN #in the event that we have fewer than top N invested into, lower weights to 1/top N
      names(longs) <- names(totalRank)
      
    } else if(length(momentum) == 1) { #only one security had positive momentum
      longs <- 1/bestN
      names(longs) <- names(momentum)
    } else { #no securities had positive momentum
      longs <- 1
      names(longs) <- riskFreeName
    }
    
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
  strategyReturns <- Return.rebalancing(R = returns, weights = weights, geometric = FALSE)
  colnames(strategyReturns) <- paste(monthLookback, weightMom, weightVol, weightCor, sep="_")
  #  return(strategyReturns)
  out <- list(weights, strategyReturns)
  return(out)
}

[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  
  
  require(PerformanceAnalytics)
  require(quantmod)
  
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
  getSymbols(mutualFunds, from="2009-06-30", to="2015-03-20")
  tmp <- list()
  for(fund in mutualFunds) {
    tmp[[fund]] <- Ad(get(fund))
  }
  
  #always use a list hwne intending to cbind/rbind large quantities of objects
  adPrices <- do.call(cbind, args = tmp)
  colnames(adPrices) <- gsub(".Adjusted", "", colnames(adPrices))
  
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
  
  replicaAttempt <- FAAreturns(adPrices)
  bestN3 <- FAAreturns(adPrices, bestN=3)
  bestN2 <- FAAreturns(adPrices, bestN=2)
  N3vol1cor1 <- FAAreturns(adPrices, weightVol = 1, weightCor = 1)
  minRisk <- FAAreturns(adPrices, weightMom = 0, weightVol=1, weightCor=1)
  pureMomentum <- FAAreturns(adPrices, weightMom=1, weightVol=0, weightCor=0)
  maxDecor <- FAAreturns(adPrices, weightMom=0, weightVol=0, weightCor=1)
  momDecor <- FAAreturns(adPrices, weightMom=1, weightVol=0, weightCor=1)
  
  all <- cbind(replicaAttempt, bestN4, bestN3, minRisk, pureMomentum, maxDecor, momDecor)
  colnames(all) <- c("Replica Attempt", "N4", "N3", "minRisk", "pureMomentum", "maxDecor", "momDecor")
  charts.PerformanceSummary(all, colorset=c("black", "red", "blue", "green", "darkgrey", "purple", "orange"))
  
  stats <- data.frame(t(rbind(Return.annualized(all)*100,
                              maxDrawdown(all)*100,
                              SharpeRatio.annualized(all))))
  stats$Return_To_Drawdown <- stats[,1]/stats[,2]
  
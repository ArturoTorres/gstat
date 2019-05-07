## turning bands:
# for cov-functions valid in 3-dim 

# adjust the covariance matrix (Schlahter, Chap. 2, eq 2.25)
# simulate several layers (up-right hyper-planes in the 3D+T-cube) with random directions on a unit hemi-sphere
# orthogonally project target points onto rndm lyrs and average the values for all rndm lyrs

# -> one simulation of a 3d+T Gaussian random field, for 2D use a hyperplane in 3D

# the tunring bands operator (numerical, 3-dimensional case)
# @input:
#  model target covariance function 
#  dist_grid: data.frame wtih columns spacelag and timelag for which the covariances are calculated
#
# @value: data.frame with columns spacelag, timelag and gamma where the latter contains the covariances for the 1-dim spatial and temporal layers

# model <- separableModel
# dist_grid <- as.data.frame(cbind("spacelag" = rep(1:150*1., each=4),
#                                  "timelag" = rep(1:4, 150)))

# CAVE: distgrid must have the correct spatial and temporal metrics

tbOperator <- function(model, dist_grid) {
  r <- dist_grid$spacelag
  
  derFun <- function(r) r * variogramSurface(model, dist_grid = data.frame(spacelag=r, timelag=dist_grid$timelag), 
                                             covariance = TRUE)$gamma
  
  cbind(dist_grid, "gamma" = diag(attr(numericDeriv(quote(derFun(r)), "r"), "gradient")))
}

randomDirections <- function(n) {
  u <- runif(n, 0, 2*pi)
  v <- runif(n, 0, 1)
  s <- sqrt(1-v^2)
  
  cbind(s*cos(u), s*sin(u), v)
}

# library(rgl)
# plot3d(randomDirections(100), aspect = c(1,1,0.5))

# computes the covariance matrixes and weights once, applied to series of
# variables/simulations where each variable/simulation is stored in one column of
# the multiVarMatrix copied from krigeST to avoid repeated calls to krige with
# multiple, identical inversions of the weights matrix

# TODO: add functionality for temporal sandwich-wise processing: i.e. use the 
# +/- nmaxTime time slices to predict one time slice

krigeSTMultiple <- function(formula, from, to, modelList, multiVarMatrix, nmax=Inf, nmaxTime=Inf, tUnit=NULL, progress=FALSE) {
  stopifnot(all(sapply(list(from, to, multiVarMatrix), function(x) inherits(x, "ST"))))
  
  if (any(c(nmax, nmaxTime) < Inf)) {
    res <- krigeSTMultiple.local(formula = formula, from = from, to = to, 
                                 modelList = modelList, 
                                 multiVarMatrix = multiVarMatrix, 
                                 nmax = nmax, nmaxTime = nmaxTime,
                                 tUnit = tUnit, progress = progress)
    
    addAttrToGeom(to, as.data.frame(res))
  } else {
    res <- krigeSTMultiple.global(formula = formula, from = from, to = to, 
                                  modelList = modelList, 
                                  multiVarMatrix = multiVarMatrix)
    addAttrToGeom(to, as.data.frame(res))
  }
}

krigeSTMultiple.global <- function(formula, from, to, modelList, multiVarMatrix) {
  lst = extractFormula(formula, from, to)
  
  separate <- length(from) > 1 && length(to) > 1 &&
    inherits(from, "STF") && inherits(to, "STF")
  
  X = lst$X
  x0 = lst$x0
  
  V = covfn.ST(from, model = modelList, separate=separate)
  v0 = covfn.ST(from, to, modelList)
  
  if (modelList$stModel == "separable" & separate)
    skwts <- STsolve(V, v0, X) # use Kronecker trick
  else 
    skwts <- CHsolve(V, cbind(v0, X))
  
  npts = length(to)
  ViX = skwts[,-(1:npts)]
  skwts = skwts[,1:npts]
  
  idPredFun <- function(sim) {
    sim <- matrix(sim, ncol = 1)
    beta = solve(t(X) %*% ViX, t(ViX) %*% sim)
    x0 %*% beta + t(skwts) %*% (sim - X %*% beta)
  }
  
  apply(multiVarMatrix@data, 2, idPredFun)
}

krigeSTMultiple.local <- function(formula, from, to, modelList, multiVarMatrix, nmax, nmaxTime, tUnit, progress) {
  
  timeInd <- index(to@time)
  nmaxTime <- switchTimeUnitToSecs(nmaxTime, tUnit)
  
  classNewdata <- class(to)
  coerceNewdataToIrregular <- !(classNewdata %in% c("STI", "STIDF"))
  newDataHasData <- "data" %in% slotNames(to)
  
  coerceDataToIrregular <- !(class(from) %in% c("STI", "STIDF"))
  
  res <- NULL
  
  # loop over all times in to
  if (progress)
    pb <- txtProgressBar(min = 0, max = length(timeInd), initial = 0, style = 3)
  
  for (i in 1:length(timeInd)) {
    if (progress)
      setTxtProgressBar(pb, i)
    
    time <- timeInd[i]
    
    timeBlock <- paste0(time + nmaxTime, collapse = "/")
    
    dataBlock <- from[ , timeBlock, drop=FALSE]
    multiVarBlock <- multiVarMatrix[ , timeBlock, drop=FALSE]
    newdataSlice <- to[ , time, drop=FALSE]
    
    # check whether timeblock is empty (sparse data, diverging alignment of to@time and from@time)
    if (nrow(dataBlock@data) == 0 || nrow(multiVarBlock@data) == 0) {
      warning(paste("Empty time block around:", time, 
                    "Predictions are set to 'NA'. Increase the tmeporal window."))
      
      emptVal <- matrix(NA, 
                        nrow = length(newdataSlice@sp),
                        ncol = ncol(multiVarMatrix@data))
      eemptVal <- as.data.frame(emptVal)
      colnames(emptVal) <- colnames(multiVarMatrix@data)
      
      res <- rbind(res, emptVal)
      
      next;
    }
    
    if (nmax < Inf) { # local ST-neighbourhood kriging
      
      # krigeST.local expect STI*s
      if(coerceDataToIrregular)
        dataBlock <- as(dataBlock, "STIDF")
      
      if(coerceNewdataToIrregular) {
        if (newDataHasData) {
          newdataSlice <- as(newdataSlice, "STIDF")
        } else {
          newdataSlice <- as(newdataSlice, "STI")
        }
      }
      
      res <- rbind(res, 
                   krigeST.local(formula = formula, 
                                 data = dataBlock, newdata = newdataSlice, 
                                 modelList = modelList,
                                 nmax = nmax, stAni = stAni, bufferNmax = bufferNmax))
      
    } else { # full ST kriging
      res <- rbind(res, 
                   krigeSTMultiple.global(formula = formula, 
                                          from = dataBlock, to = newdataSlice, 
                                          multiVarMatrix = multiVarBlock, 
                                          modelList = modelList))
    }  
  }
  
  return(res)
}

# krigeST <- function(formula, data, newdata, modelList, y, beta, nmax=Inf, stAni=NULL,
#                     computeVar = FALSE, fullCovariance = FALSE,
#                     bufferNmax=2, progress=TRUE)

## 
krigeSTSimTB <- function(formula, data, newdata, modelList, nsim, 
                         progress=TRUE, nLyrs=500, tGrid=NULL, sGrid=NULL, ceExt=2,
                         nmax=Inf, nmaxTime=Inf) {
  stopifnot(zoo::is.regular(newdata@time))
  
  tUnitModel <- attr(modelList, "temporal unit")
  tUnitData <- units(diff(c(index(newdata@time[1]), newdata@endTime[1])))
  if (is.null(tUnitModel)) {
    warning("The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.\n The unit '", tUnitData,
            "' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.")
    tUnit <- tUnitData
    attr(modelList, "temporal unit") <- tUnit
  } else {
    tUnit <- tUnitModel
    message("[Using the following time unit: ", tUnit, "]")
  }
  
  condSim <- TRUE
  
  if (missing(data)) {
    condSim <- FALSE
    message("[No data provided: performing unconditional simulation.]")
  } else {
    message("[Performing conditional simulation.]")
  }
  
  pb <- txtProgressBar(0,nsim,style=3)
  
  # ST-simulation grid
  if (is.null(tGrid)) {
    tDis <- diff(c(index(newdata@time[1]), newdata@endTime[1]))
    units(tDis) <- tUnit
    
    tGrid <- c(as.numeric(tDis), length(newdata@time))
    attr(tGrid, "units") <- c("", tUnit)
  }
  
  # ST-simulation grid
  if (is.null(sGrid)) {
    if (gridded(newdata@sp)) {
      # SpatialPixels/SpatialGrid:
      # based on GridTopology: use minimal cellsize of both directions; take enough to cover the diagonal
      sDis <- min(newdata@sp@grid@cellsize)
      sDim <- ceiling(sqrt(sum((newdata@sp@grid@cellsize * newdata@sp@grid@cells.dim)^2))/sDis)
      sGrid <-  c(sDis, sDim)
    } else {
      # treat (as) SpatialPoints:
      # Average area per location --assuming-a-regular-squared-outline-taking-the-sqrt--> length/location --take-later-twice-as-many-->
      bboxExt <- apply(newdata@sp@bbox, 1, diff)
      sDis <- sqrt(prod(bboxExt)/length(newdata@sp))
      sDim <- ceiling(sqrt(sum((bboxExt/sDis)^2)))
      sGrid <- c(sDis, sDim)
    }
  }
  
  # random directions in 3D
  rndDir <- randomDirections(nLyrs)
  
  # coordinates 
  coordMat <- coordinates(newdata@sp)
  
  # coordinates (embedded in 3D) shifted + scaled to the grid index of the spatial gridding sGrid
  coordMat[,1] <- (coordMat[,1] - newdata@sp@bbox[1,1])/sGrid[1]
  coordMat[,2] <- (coordMat[,2] - newdata@sp@bbox[2,1])/sGrid[1]
  
  if (ncol(coordMat) == 2) {
    coordMat <- cbind(coordMat,1)
  } else {
    coordMat[,3] <- (coordMat[,3] - newdata@sp@bbox[3,1])/sGrid[1]
  }
  
  # how much does each direction contribute to the point
  cntrbtn <- rndDir %*% t(coordMat)
  # -> each column corresponds to one location, each row i to the "contribution" of the i-th rndDir
  
  # ceiling and floor + shifted to avoid negative indices
  cl_cntrbtn <- ceiling(cntrbtn) + sGrid[2]
  fl_cntrbtn <- floor(cntrbtn) + sGrid[2]
  
  covRow1 <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(sGrid[1], 2*sGrid[2]), tGrid, modelList, turningLayers = TRUE, ext=ceExt)
  
  origDim <- c(2*sGrid[2], tGrid[2])
  sims <- list()
  
  ##
  for (i in 1:nsim) {
    setTxtProgressBar(pb, i)
    simLyrs <- ceSim(covRow1, nLyrs, origDim)
    simLyrs <- lapply(1:nLyrs, function(col) matrix(simLyrs[,col], 
                                                    nrow=origDim[1], 
                                                    ncol=origDim[2]))
    lambda <- cl_cntrbtn - cntrbtn - sGrid[2]
    
    cntrbSngSimLyr <- function(lyrId) {
      simLyrs[[lyrId]][cl_cntrbtn[lyrId,],] * (1-lambda)[lyrId,] + simLyrs[[lyrId]][fl_cntrbtn[lyrId,],] * lambda[lyrId,]
    }
    
    # reduce to the one realisation based on the combination of nLyrs turning bands
    sims[[paste0("sim",i)]] <- Reduce('+', lapply(1:nLyrs, cntrbSngSimLyr))/sqrt(nLyrs)
  }
  close(pb)
  
  sims <- do.call(cbind, lapply(sims, as.numeric))
  
  newSimData <- newdata
  
  # bind simulations to newdata geometry
  if ("data" %in% slotNames(newSimData))
    newSimData@data <- cbind(newSimData@data, sims)
  else
    newSimData <- addAttrToGeom(newSimData, as.data.frame(sims))
  
  # if unconditional, stop function by returning newSimData STxDF 
  if (!condSim) return(newSimData)
  
  ######################
  ## conditional case ##
  ######################
  
  varName <- all.vars(formula[[2]])
  
  ## conditioning
  # interpolate the observations to the simulation grid
  obsMeanField <- krigeST(formula = formula, 
                          data = data, newdata = newdata,
                          modelList = modelList,
                          nmax = nmax, nmaxTime = nmaxTime)
  
  # interpolate to observation locations from the simulated grids for each simulation
  simMeanObsLoc <- krigeSTMultiple(as.formula(paste0("sim1 ~", formula[[3]])),
                                   newSimData, data, 
                                   modelList = modelList, 
                                   multiVarMatrix = newSimData,
                                   nmax = nmax, nmaxTime = nmaxTime, tUnit = tUnit)
  
  # interpolate from kriged mean sim at observed locations back to the grid for mean surface of the simulations
  simMeanFields <- krigeSTMultiple(as.formula(paste0("sim1 ~", formula[[3]])), 
                                   simMeanObsLoc, newdata, 
                                   modelList = modelList,
                                   multiVarMatrix =  simMeanObsLoc,
                                   nmax = nmax, nmaxTime = nmaxTime, tUnit = tUnit)
  
  # add up the mean field and the corrected data
  sims <- obsMeanField@data$var1.pred + sims - simMeanFields@data
  
  # bind simulations to newdata geometry
  if ("data" %in% slotNames(newdata)) {
    newdata@data <- cbind(newdata@data, sims)
    return(newdata)
  }
  
  addAttrToGeom(newdata, as.data.frame(sims))
}

# 
# sTime <- Sys.time()
# krigedSim <- krigeSTUncSimTB(stf, metricModel, 100)
# Sys.time() - sTime
# 
# # 27 secs for 100 simulated ST fields of 155 locations and 21 time steps: 325500 values
# 
# # plot one simulation along time
# stplot(krigedSim[,1:12])
# 
# # plot one simulation along time
# stplot(krigedSim[1:12,,"sim1"], mode="ts")
# 
# # plot the ten simulations of the first day
# spplot(krigedSim[,1], paste0("sim",1:10), as.table=TRUE)
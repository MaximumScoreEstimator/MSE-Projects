# Companion code for the paper "Stars In Their Constellations: Great person or great team?" 
# by Mindruta, Bercovitz, Mares, and Feldman.
# 
# This code generates the contribution intervals (max and min) of an agent in a match
# The code reads two pieces of information
# First, the code requies a data input file containing 
# a) the index of matching markets (1 to 33)
# b) the market-specific indices of upstream (PIs) and downstream (constellations) agents in each market (2 columns)
# c) the matching covariates 
# d) a variable taking the value of 1 for matched agents in a market (and 0 otherwise)
# e) two additional variables taking value of 1 (and 0 otherwise) for the upstream and downstream agents that need to be removed
# Second, the code requires the list of point estimates of the matching function. 
# Here, we use the estimates produced in Mathematica: see Table 2a in the paper.
# Users can apply the code to their input data, but need to follow the same data file format.
# The input file should always be sorted by marketid, upstreamid, downstreamid. The code does not support missing data.

library(maxscoreest)

shiftPayoffMatrices <- function(payoffMatrices) {
    results <- lapply(
        payoffMatrices,
        function(payoffMatrix) {
            payoffMin <- min(payoffMatrix)
            return(list(payoffMatrix = payoffMatrix - payoffMin, offset = payoffMin))
        })
    newPayoffMatrices <- lapply(results, "[[", "payoffMatrix")
    offsets <- lapply(results, "[[", "offset")
    return(list(payoffMatrices = newPayoffMatrices, offsets = unlist(offsets)))
}

modifyR <- function(uIdxs, dIdxs, payoffMatrix, quotaU, quotaD) {
    quotaU[uIdxs] <- 0
    quotaD[dIdxs] <- 0
    matchMatrix <- generateAssignmentMatrix(payoffMatrix, quotaU, quotaD)
    return(list(
        matchMatrix = matchMatrix, quotaU = quotaU, quotaD = quotaD))
}

# TODO: Let importMatched take an already imported data.table as input.
# Maybe use the unexported developer functions.
# TODO: Use read.table?
# TODO: Allow user to pass extra args to reader functions in importMatched, etc.

ifname <- "~/Documents/Work/MSE-R/modifydata/stars_replication_remove.dat"

# @P: add a short descriptino of the precomputed file. Warn the user the file will be generated and saved each time the code runs

ofname <- "~/Documents/Work/MSE-R/modifydata/precomputed.dat"

# @P Warn the user that reading input data is file-specific. Numbers below will change depending on the number of columns

inputData <- read.csv(ifname, header = TRUE, sep = "\t")
outputData <- inputData[1:11]
removeData <- inputData[c(1:3, 11:13)]
write.table(outputData, file = ofname, sep = "\t", quote = FALSE, row.names = FALSE)
matchedData <- importMatched(ofname)
#@P cmarket is also file specific. The current name is "marketid". The agent indices are called "upstreamid" "downstreamid"
removeu <- unique(removeData[removeData$removeup == 1, c("cmarket", "cid_u")])
removed <- unique(removeData[removeData$removedown == 1, c("cmarket", "cid_d")])
removes <- removeData[, c("removeup", "removedown")]
ineqmembers <- Cineqmembers(matchedData$mate)
dataArray <- CdataArray(matchedData$distanceMatrices, ineqmembers)
# Input point estimates of the matching function below
pointEstimate <- c(
    4.981905468469215,
    -2.4016513033944378,
    0.06642024971697014,
    1.2889501155821126,
    1.3983116175519597,
    2.76088191399337)
obj <- makeScoreObjFun(dataArray, 1, 1)
print(obj(pointEstimate) * dim(dataArray)[2])
payoffMatrices <- evaluatePayoffMatrices(matchedData$distanceMatrices, pointEstimate)
shiftResults <- shiftPayoffMatrices(payoffMatrices)
payoffMatricesOld <- payoffMatrices
payoffMatrices <- shiftResults$payoffMatrices
quotasU <- lapply(matchedData$matchMatrices, colSums)
quotasD <- lapply(matchedData$matchMatrices, rowSums)
totalPayoffs <- mapply(function(m, p) { sum(m*p) }, matchedData$matchMatrices, payoffMatrices)
umatchesHeader <- list()
dmatchesHeader <- list()
umatches <- list()
dmatches <- list()
ures1 <- list()
dres1 <- list()
for (i in seq_len(dim(removeu)[1])) {
    mIdx <- removeu[i, "cmarket"]
    uIdx <- removeu[i, "cid_u"]
    modifyResult <- modifyR(c(uIdx), c(), payoffMatrices[[mIdx]], quotasU[[mIdx]], quotasD[[mIdx]])
    col <- list(stream = "U", mIdx = mIdx, uIdx = uIdx)
    newMatchMatrices <- matchedData$matchMatrices
    newMatchMatrices[[mIdx]] <- modifyResult$matchMatrix
    newTotalPayoff <- sum(modifyResult$matchMatrix * payoffMatrices[[mIdx]])
    totalPayoffDiff <- totalPayoffs[mIdx] - newTotalPayoff
    umatchesHeader <- append(umatchesHeader, list(col))
    # umatches <- append(umatches, newMatchMatrices)
    umatches <- append(umatches, list(modifyResult$matchMatrix))
    out <- list(mIdx = mIdx, uIdx = uIdx, totalPayoffDiff = totalPayoffDiff)
    ures1 <- append(ures1, list(out))
}
for (i in seq_len(dim(removed)[1])) {
    mIdx <- removed[i, "cmarket"]
    dIdx <- removed[i, "cid_d"]
    modifyResult <- modifyR(c(), c(dIdx), payoffMatrices[[mIdx]], quotasU[[mIdx]], quotasD[[mIdx]])
    col <- list(stream = "D", mIdx = mIdx, dIdx = dIdx)
    newMatchMatrices <- matchedData$matchMatrices
    newMatchMatrices[[mIdx]] <- modifyResult$matchMatrix
    newTotalPayoff <- sum(modifyResult$matchMatrix * payoffMatrices[[mIdx]])
    totalPayoffDiff <- totalPayoffs[mIdx] - newTotalPayoff
    dmatchesHeader <- append(dmatchesHeader, list(col))
    # dmatches <- append(dmatches, newMatchMatrices)
    dmatches <- append(dmatches, list(modifyResult$matchMatrix))
    out <- list(mIdx = mIdx, dIdx = dIdx, totalPayoffDiff = totalPayoffDiff)
    dres1 <- append(dres1, list(out))
}
k <- 1
v1 <- list() # maybe use a data.frame
for (mIdx in seq_len(matchedData$noM)) {
    for (uIdx in seq_len(matchedData$noU[mIdx])) {
        for (dIdx in seq_len(matchedData$noD[mIdx])) {
            payoff <- payoffMatrices[[mIdx]][dIdx, uIdx]
            originalMatch <- matchedData$matchMatrices[[mIdx]][dIdx, uIdx]
            match <- matchedData$matchMatrices[[mIdx]][dIdx, uIdx]
            idxs <- which(sapply(ures1, function(v) { v$mIdx == mIdx && v$uIdx == uIdx }))
            urmax <- if (length(idxs) == 0) NA else ures1[[idxs[1]]]$totalPayoffDiff
            # urmin <- if (length(idxs) == 0) NA else payoff - urmax
            idxs <- which(sapply(dres1, function(v) { v$mIdx == mIdx && v$dIdx == dIdx }))
            drmax <- if (length(idxs) == 0) NA else dres1[[idxs[1]]]$totalPayoffDiff
            # drmin <- if (length(idxs) == 0) NA else payoff - drmax
            urmin <- max(0, payoff - drmax)
            drmin <- max(0, payoff - urmax)
            row <- list(
                mIdx = mIdx, uIdx = uIdx, dIdx = dIdx,
                payoff = payoff, originalMatch = originalMatch, match = match,
                urmax = urmax, drmax = drmax, urmin = urmin, drmin = drmin,
                removeU = removes[k, 1], removeD = removes[k, 2])
            k <- k + 1
            v1 <- append(v1, list(row))
        }
    }
}
v1 <- data.frame(do.call(rbind, v1))
colnames(v1)[1:3] <- matchedData$header[1:3]
colnames(v1)[5:6] <- c("originalmatch", "storedmatch")
v1u <- v1[v1$originalmatch == 1 & v1$removeU == 1, ]
v1d <- v1[v1$originalmatch == 1 & v1$removeD == 1, ]
# FIXME: Columns of v1 and v2 are lists, why?
# TODO: Surround field names in header with quotes?
# Create columns for v2
exportMIdx <- 33
v2 <- v1
idxs <- which(sapply(umatchesHeader, function(x) { x$mIdx == exportMIdx }))
for (i in seq_along(idxs)) {
    idx <- idxs[i]
    colDescription <- umatchesHeader[[idx]]
    col <- integer(dim(v1)[1])
    k <- 1
    for (mIdx in seq_len(matchedData$noM)) {
        matchMatrix <- if (mIdx == colDescription$mIdx) umatches[[idx]] else matchedData$matchMatrices[[mIdx]]
        for (uIdx in seq_len(matchedData$noU[mIdx])) {
            for (dIdx in seq_len(matchedData$noD[mIdx])) {
                col[k] <- matchMatrix[dIdx, uIdx]
                k <- k + 1
            }
        }
    }
    colName <- sprintf("{%s, %d, %d}", colDescription$stream, colDescription$mIdx, colDescription$uIdx)
    v2[colName] <- col
}
idxs <- which(sapply(dmatchesHeader, function(x) { x$mIdx == exportMIdx }))
for (i in seq_along(idxs)) {
    idx <- idxs[i]
    colDescription <- dmatchesHeader[[idx]]
    col <- integer(dim(v1)[1])
    k <- 1
    for (mIdx in seq_len(matchedData$noM)) {
        matchMatrix <- if (mIdx == colDescription$mIdx) dmatches[[idx]] else matchedData$matchMatrices[[mIdx]]
        for (uIdx in seq_len(matchedData$noU[mIdx])) {
            for (dIdx in seq_len(matchedData$noD[mIdx])) {
                col[k] <- matchMatrix[dIdx, uIdx]
                k <- k + 1
            }
        }
    }
    colName <- sprintf("{%s, %d, %d}", colDescription$stream, colDescription$mIdx, colDescription$dIdx)
    v2[colName] <- col
}

# The code exports 3 files: ".removeu.originalmatch.dat" ".removed.originalmatch.dat" ".all.dat"
# The first two files keep only the agents in the market for which we calculated the contribution intervals
# These files have been used in the paper to generate the summary statistics reported in Table 3a
# The code also allows exporting a file (".all.dat") that contains all agents in a market (see section 5.3 of the paper)
# This file contains additional columns that mark down "who matches with whom" after an agent is removed from the market. 
# For example, the column called "{U, 33, 6}" takes the value of 1 
# for matches that are formed in market 33 after upstream 6 was "removed" from the market (i.e. its quota was set to zero) 
# The column called "{D, 33, 2}" takes the value of 1 for all matches formed in market 33 after downstream 2 was "removed" 


ofPrefix <- file.path(dirname(ifname), paste("modifyR-output-", basename(ifname), sep = ""))
ofnamev1u <- paste(ofPrefix, ".removeu.originalmatch.dat", sep = "")
ofnamev1d <- paste(ofPrefix, ".removed.originalmatch.dat", sep = "")
ofnamev2 <- paste(ofPrefix, ".all.dat", sep = "")
# write.table(v1u, file = ofnamev1u, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(v1u, file = ofnamev1u, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(v1d, file = ofnamev1d, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(v2, file = ofnamev2, sep = "\t", quote = FALSE, row.names = FALSE)

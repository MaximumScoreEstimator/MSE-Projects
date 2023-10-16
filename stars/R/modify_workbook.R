# Companion code for the paper "Stars In Their Constellations: Great person or great team?"
# by Mindruta, Bercovitz, Mares, and Feldman.
#
# This code generates the contribution intervals (max and min) of an agent in a match
# The code reads two pieces of information
# First, a data input file containing
# a) the index of the matching markets (1 to 33)
# b) the market-specific indices of upstream (PIs) and downstream (constellations) agents in each market (2 columns)
# c) the matching covariates
# d) a variable taking the value of 1 for matched agents in a market (and 0 otherwise)
# e) two additional variables taking value of 1 (and 0 otherwise) for the upstream and downstream agents that need to be "removed"
# Second, the code requires the list of point estimates of the matching function.
# Here, we use the estimates produced in Mathematica: see Table 2a in the paper.
# Users can apply the code to their input data, but need to follow the same data file format.
# The code does not support missing data.

library(maxscoreest)

################################################################################
# Inputs for the workbook.
# Can be replaced by the user.

# The input file must be a delimiter-separated file with a header, and must
# have the same structure as described in the `maxscoreest::importMatched`
# function, but with two additional columns at the end. Fields on these columns
# take values `0` or `1`, depending on whether the corresponding row has to be
# unmatched.
# More specifically:
# * The first three columns specify the market index, the upstream index, and
#   the downstream index, in this order.
# * Each market/upstream/downstream triple should appear only once. These
#   indices should have sorted and consecutive values, starting from `1`.
# * The last three columns specify whether this triple is a match, and whether
#   this upstream or downstream should be unmatched.
#   These fields take values `0` or `1`.
# * Each of the remaining columns in the middle corresponds to a matching covariate.
#   Fields on these columns should be numerical.
# * The input file should always be sorted by marketid, upstreamid, downstreamid.
#   Note that the actual header names are not taken into account.

# Path to the input file.
ifname <- "../data/stars_replication_step2.dat"

# * The optimal parameters *
# Attention, for reproducibility reasons, we provide the estimates obtained in Mathematica
pointEstimate <- c(
  0.8647185894805215,
  5.779401410940838,
  -2.6456032704736994,
  3.3901711331785256,
  1.6199736480992488,
  2.0775031924093352)

# *The index of the market to be exported after removing some agents and rematching the existing ones*
# This allows the user to inspect how the remaining agents rematch after some disappear from the market
exportMIdx <- 33
# Whether to add the `storedmatch` column to the ".all.dat" output file.
useStoredMatch <- FALSE

################################################################################

################################################################################
# Function definitions.

#' Shift payoff matrices to zero
#'
#' Transform each payoff matrix by subtracting from each element the minimal
#' entry of that matrix (called "offset"). The new payoff matrices have minimal
#' value equal to zero Return the shifted payoff matrices and the calculated
#' offsets.
#'
#' This part of the code executes an affine transformation on the payoff matrix
#' (the matrix that contains the payoffs of all real and counterfactual
#' matches). This procedure does not interfere with the calculation of max and
#' min, but speeds up the execution.
#'
#' @param payoffMatrices A list of payoff matrices, one for each market.
#' @return A list with members:
#' \tabular{ll}{
#'   `$payoffMatrices` \tab The shifted payoff matrices, as a list. \cr
#'   `$offsets` \tab The offsets, as a list.
#' }
shiftPayoffMatrices <- function(payoffMatrices) {
    results <- lapply(
        payoffMatrices,
        function(payoffMatrix) {
            payoffMin <- min(payoffMatrix)
            return(list(
                payoffMatrix = payoffMatrix - payoffMin, offset = payoffMin))
        })
    newPayoffMatrices <- lapply(results, "[[", "payoffMatrix")
    offsets <- lapply(results, "[[", "offset")
    return(list(payoffMatrices = newPayoffMatrices, offsets = unlist(offsets)))
}

#' Modify market data by unmatching and rematching
#'
#' Unmatch the given upstreams and downstreams by setting their respective
#' quotas to zero. These agents cannot participate in subsequent matches. After unmatching,
#' calculate the new optimal matches in each market, given the original payoff matrix.
#'
#' @param uIdxs A vector of upstream indices to be unmatched.
#' @param dIdxs A vector of downstream indices to be unmatched.
#' @param payoffMatrix The payoff matrix of the market.
#' @param quotaU The vector of quotas for the upstreams.
#' @param quotaD The vector of quotas for the downstreams.
#' @return A list with members:
#' \tabular{ll}{
#'   `$matchMatrix` \tab The match matrix after rematching. \cr
#'   `$quotaU` \tab The vector of upstream quotas after unmatching. \cr
#'   `$quotaD` \tab The vector of downstream quotas after unmatching.
#' }
modifyR <- function(uIdxs, dIdxs, payoffMatrix, quotaU, quotaD) {
    quotaU[uIdxs] <- 0
    quotaD[dIdxs] <- 0
    matchMatrix <- generateAssignmentMatrix(payoffMatrix, quotaU, quotaD)
    return(list(
        matchMatrix = matchMatrix, quotaU = quotaU, quotaD = quotaD))
}

################################################################################

################################################################################
# Create output paths.
# Output files are created in the same directory as the input file.

# `precomputed.dat` is created by removing the last two columns from the input
# file. This file is generated locally each time this code runs. Note that this is not
# an output file, but a helper file for the rest of the procedure, and it is safe
# to ignore. Do not delete this file until all results are obtained.

ofname <- file.path(dirname(ifname), "precomputed.dat")
ofPrefix <- file.path(
    dirname(ifname), paste("modifyR-output-", basename(ifname), sep = ""))
ofnamev1u <- paste(ofPrefix, ".removeu.originalmatch.dat", sep = "")
ofnamev1d <- paste(ofPrefix, ".removed.originalmatch.dat", sep = "")
ofnamev2 <- paste(ofPrefix, ".all.dat", sep = "")

################################################################################

################################################################################
# Extract precomputed and removal data.
# See section "Inputs for the workbook" for information on the structure of the
# input file.

inputData <- read.csv(ifname, header = TRUE, sep = "\t")
numInputCols <- length(colnames(inputData))
outputData <- inputData[1:(numInputCols - 2)]
removeData <- inputData[c(1:3, (numInputCols - 2):numInputCols)]
write.table(
    outputData, file = ofname, sep = "\t", quote = FALSE, row.names = FALSE)

marketColName <- colnames(removeData)[1]
upstreamColName <- colnames(removeData)[2]
downstreamColName <- colnames(removeData)[3]
removeupColName <- colnames(removeData)[5]
removedownColName <- colnames(removeData)[6]
# `removeu` is a table of pairs of market and upstream indices to be unmatched.
# `removed` is similar, but for downstreams.
removeu <- unique(removeData[
    removeData[[removeupColName]] == 1,
    c(marketColName, upstreamColName)])
removed <- unique(removeData[
    removeData[[removedownColName]] == 1,
    c(marketColName, downstreamColName)])
removes <- removeData[, c(removeupColName, removedownColName)]

################################################################################

################################################################################
# Calculate associated data (match and payoff matrices, quotas, and total
# payoffs). Shift payoff matrices to zero minimum.
# See the `maxscoreest` documentation for more information on their format.
matchedData <- importMatched(ofname)
ineqmembers <- Cineqmembers(matchedData$mate)
dataArray <- CdataArray(matchedData$distanceMatrices, ineqmembers)
payoffMatrices <- evaluatePayoffMatrices(
    matchedData$distanceMatrices, pointEstimate)
shiftResults <- shiftPayoffMatrices(payoffMatrices)
payoffMatrices <- shiftResults$payoffMatrices
quotasU <- lapply(matchedData$matchMatrices, colSums)
quotasD <- lapply(matchedData$matchMatrices, rowSums)
# The total payoff of a market is calculated by summing the payoffs of those
# upstream-downstream pairs that are matched.
totalPayoffs <- mapply(
    function(m, p) { sum(m*p) }, matchedData$matchMatrices, payoffMatrices)
# Match matrices obtained after breaking all existing matches and calculating
# the matchings that maximize the total payoffs given the payoff matrices
# obtained from the optimal parameters.
storedMatchMatrices <- CmatchMatrices(payoffMatrices, quotasU, quotasD)

################################################################################

################################################################################
# Perform unmatchings and accumulate results.

#' Unmatch and rematch all pairs
#'
#' Perform a modification for each row in the `removeu` or `removed` tables.
#' When the upstream functionality is chosen, each modification unmatches only
#' the single upstream described in the `removeu` row, and unmatches no
#' downstreams. When the downstream functionality is chosen, each modification
#' unmatches the single downstream described in the `removed` row, and unmatches
#' no upstreams. In both cases, the new match matrix and the difference in total
#' payoffs is recorded.
#'
#' @param stream `"U"` for upstreams, `"D"` for downstreams.
#' @return A list with one element for each row of `removeu` or `removed`,
#'   respectively. Each element is a list with members:
#' \tabular{ll}{
#'   `$stream` \tab The value of the `stream` parameter. \cr
#'   `$mIdx` \tab The index of the market to be modified. \cr
#'   `$uIdx` or `$dIdx` \tab The index of the upstream or downstream,
#'     respectively, that is unmatched. \cr
#'   `$matches` \tab The new match matrix, after rematching. \cr
#'   `$totalPayoffDiff` \tab The difference in total payoffs before and after
#'     rematching. This is referred to in the paper as the maximum of the contribution interval.
#' }
calcRemoveResults <- function(stream) {
    stopifnot(stream == "U" || stream == "D")
    removeTable <- if (stream == "U") removeu else removed
    streamColName <- if (stream == "U") upstreamColName else downstreamColName
    sIdxName <- if (stream == "U") "uIdx" else "dIdx"
    lapply(seq_len(dim(removeTable)[1]), function(removeTableRowIdx) {
        mIdx <- removeTable[removeTableRowIdx, marketColName]
        sIdx <- removeTable[removeTableRowIdx, streamColName]
        uIdxs <- if (stream == "U") c(sIdx) else c()
        dIdxs <- if (stream == "U") c() else c(sIdx)
        modifyResult <- modifyR(
            uIdxs, dIdxs,
            payoffMatrices[[mIdx]], quotasU[[mIdx]], quotasD[[mIdx]])
        newTotalPayoff <- sum(modifyResult$matchMatrix * payoffMatrices[[mIdx]])
        totalPayoffDiff <- totalPayoffs[mIdx] - newTotalPayoff
        out <- list(
            stream = stream, mIdx = mIdx, matches = modifyResult$matchMatrix,
            totalPayoffDiff = totalPayoffDiff)
        out[[sIdxName]] <- sIdx
        return(out)
    })
}
removeResultsU <- calcRemoveResults("U")
removeResultsD <- calcRemoveResults("D")

################################################################################

################################################################################
# Create output tables.
# `v1` consists of one row for each input row, with entries for market,
# upstream, and downstream indices, the triple's payoff value, and whether it
# was originally matched. If that upstream or downstream was removed above,
# there is extra information regarding the differences in total payoffs.
# The last two columns are the same as in the input data.
# From `v1` we select those rows that were originally matched and were then
# removed, and we construct the tables `v1u` and `v1d`.
# We finally construct `v2` by appending columns to `v1` corresponding to the
# removed upstreams and downstreams from the chosen market (see `exportMIdx`).

rowIdx <- 1
v1 <- list()
for (mIdx in seq_len(matchedData$noM)) {
    for (uIdx in seq_len(matchedData$noU[mIdx])) {
        for (dIdx in seq_len(matchedData$noD[mIdx])) {
            payoff <- payoffMatrices[[mIdx]][dIdx, uIdx]
            originalMatch <- matchedData$matchMatrices[[mIdx]][dIdx, uIdx]
            storedMatch <- storedMatchMatrices[[mIdx]][dIdx, uIdx]
            # Check if that upstream was previously removed.
            idxs <- which(sapply(
                removeResultsU,
                function(v) { v$mIdx == mIdx && v$uIdx == uIdx }))
            urmax <- if (length(idxs) == 0) {
                NA
            } else {
                removeResultsU[[idxs[1]]]$totalPayoffDiff
            }
            # Check if that downstream was previously removed.
            idxs <- which(sapply(
                removeResultsD,
                function(v) { v$mIdx == mIdx && v$dIdx == dIdx }))
            drmax <- if (length(idxs) == 0) {
                NA
            } else {
                removeResultsD[[idxs[1]]]$totalPayoffDiff
            }
            urmax <- max(0, urmax)
            drmax <- max(0, drmax)
            urmin <- max(0, payoff - drmax)
            drmin <- max(0, payoff - urmax)
            row <- list(
                mIdx = mIdx, uIdx = uIdx, dIdx = dIdx,
                payoff = payoff, originalMatch = originalMatch,
                storedMatch = storedMatch,
                urmax = urmax, drmax = drmax, urmin = urmin, drmin = drmin,
                removeU = removes[rowIdx, 1], removeD = removes[rowIdx, 2])
            if (!useStoredMatch) {
                # This actually drops the list element.
                row$storedMatch <- NULL
            }
            rowIdx <- rowIdx + 1
            v1 <- append(v1, list(row))
        }
    }
}
v1 <- data.frame(rbindlist(v1))
colnames(v1)[1:3] <- matchedData$header[1:3]
colnames(v1)[5] <- c("originalmatch")
if (useStoredMatch) {
    colnames(v1)[6] <- c("storedmatch")
}
v2 <- v1
if (useStoredMatch) {
    # Only add the `storedmatch` column to ".all.dat".
    v1 <- v1[, -6]
}
v1u <- v1[v1$originalmatch == 1 & v1$removeU == 1, ]
v1d <- v1[v1$originalmatch == 1 & v1$removeD == 1, ]

#' Append columns to output table
#'
#' When the upstream functionality is chosen, the rows of `removeu` with
#' market index equal to the chosen market index `exportMIdx` are selected. For
#' each upstream index in those rows, a new column for `v2` is created. Each
#' entry of that column is the value of the new match matrix, after rematching,
#' for the corresponding market-upstream-downstream index. The column header is
#' of the form `"{U, mIdx, uIdx}"`, where `mIdx` and `uIdx` are replaced by the
#' market and upstream indices of the removed upstream. The procedure for the
#' downstream functionality is similar.
#'
#' @param stream `"U"` for upstreams, `"D"` for downstreams.
#' @return `NULL`
appendColumns <- function(stream) {
    stopifnot(stream == "U" || stream == "D")
    removeResults <- if (stream == "U") removeResultsU else removeResultsD
    sIdxName <- if (stream == "U") "uIdx" else "dIdx"
    idxs <- which(sapply(removeResults, function(x) { x$mIdx == exportMIdx }))
    lapply(idxs, function(idx) {
        colDescription <- removeResults[[idx]][c("stream", "mIdx", sIdxName)]
        col <- integer(dim(v1)[1])
        rowIdx <- 1
        for (mIdx in seq_len(matchedData$noM)) {
            # The match matrix is only modified for the chosen market. For other
            # markets we can retrieve the original match matrices.
            matchMatrix <- if (mIdx == colDescription$mIdx) {
                removeResults[[idx]]$matches
            } else {
                matchedData$matchMatrices[[mIdx]]
            }
            for (uIdx in seq_len(matchedData$noU[mIdx])) {
                for (dIdx in seq_len(matchedData$noD[mIdx])) {
                    col[rowIdx] <- matchMatrix[dIdx, uIdx]
                    rowIdx <- rowIdx + 1
                }
            }
        }
        colName <- sprintf(
            "{%s, %d, %d}", colDescription$stream,
            colDescription$mIdx, colDescription[[sIdxName]])
        v2[colName] <<- col
        return(NULL)
    })
    return(NULL)
}
appendColumns("U")
appendColumns("D")

################################################################################

################################################################################
# Report number of instances where min exceeds max.

print(sprintf(
    "Number of instances where urmin > drmax: %d",
    sum(v1[, "originalmatch"] == 1 & v1[, "removeU"] == 1 & v1[, "urmin"] > v1[, "drmax"], na.rm = TRUE)))
print(sprintf(
    "Number of instances where drmin > urmax: %d",
    sum(v1[, "originalmatch"] == 1 & v1[, "removeD"] == 1 & v1[, "drmin"] > v1[, "urmax"], na.rm = TRUE)))

################################################################################

################################################################################
# Write output tables.

# The code exports 3 files: ".removeu.originalmatch.dat" ".removed.originalmatch.dat" ".all.dat"
# The first two files keep only the agents in the market for which we calculated the contribution intervals
# These files have been used in the paper to generate the summary statistics reported in Table 3a
# The code also allows exporting a file (".all.dat") that contains all agents in a market (see section 5.3 of the paper)
# This file contains additional columns that mark down "who matches with whom" after an agent is removed from the market.
# For example, the column called "{U, 33, 6}" takes the value of 1
# for matches that are formed in market 33 after upstream 6 was "removed" from the market (i.e. its quota was set to zero)
# The column called "{D, 33, 2}" takes the value of 1 for all matches formed in market 33 after downstream 2 was "removed"

fwrite(v1u, file = ofnamev1u, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(v1d, file = ofnamev1d, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(v2, file = ofnamev2, sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################

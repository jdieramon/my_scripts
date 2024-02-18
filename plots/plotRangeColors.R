# plotRangesColors.R -- a function for plotting IRanges data with different colors
# Copyright (C) 2017 Jose V. Die  <jodiera@upv.es>
# Distributed under terms of the MIT license.


## You may want to read a detailed tutorial on this function on Rpubs: 
## http://rpubs.com/JoseVDie/322171

# Date = date()

## --- Dependencies ------------------------------------------------------------
library(IRanges)

## --- Function to plot IRanges-------------------------------------------------
# Source: https://github.com/genomicsclass/ph525x/blob/master/R/plotRanges.R

## =============================================================================
## Example 1: A single IRange object containing exons and primer coordinates
## =============================================================================

##--- Gen1: PPIB -----------------------------------------------------------------
# PPIB <- IRanges(start=c(15,808,1571,4643,5574, 4806, 5617),
#                end=c(125,921,1664,4827,5696, 4825, 5636), 
#                names=c("exon1", "exon2", "exon3", "exon4", "exon5",
#                        "Forward", "Reverse"))
# PPIB


## --- Create a vector that indexes the IRanges object containing the primer ---
# target <- c("Forward", "Reverse")

color_ind <- function(ir, target) {
  val <- which(names(ir) %in% target)
  val
}

# c_ind <- color_ind(PPIB, target)

## --- Define a new function that reads the IRanges object and color_ind output
getColor <- function(ir, ind, color1, color2) {
  
  color_output <- vector("double", length(ir))
  
  for(i in seq(length(ir))) {
    if(i %in% ind) {
      color_output[i] = color2}
    else {color_output[i] = color1}
  }
  color_output
  
}


# Call that function. 
# color = getColor(PPIB, c_ind, 1, 2)

## --- Modify ´plotRanges´ to color the selected ranges ---
plotRangesColor <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = color, ## modified
                       sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, 
       border = col, ...)   ## modified
  title(main)
  axis(1)
  par(las = 1)
}

## Plot Iranges object
# plotRangesColor(PPIB, col = color, xlim=c(0,5906), main="PPIB")


## =============================================================================
## Example 2: Two IRange objects containing exons and primer coordinates
## =============================================================================

##  ---Gen2: DBNDD2 --------------------------------------------------------------
# DBNDD2 <- IRanges(start=c(1889,2217,3680), 
#                  end=c(2026,2354,3888), 
#                  names = c(paste(rep("exon"),1:3, sep="")))

# primer.DBNDD2 = IRanges(start=c(2270,2324), end=c(2289,2343), 
#                        names = c("Forward", "Reverse"))

## Create a sorted IRanges combaining both objects. Then, the rest of the procedure is similar as Example1.   
# sort_ir <- sort(c(DBNDD2, primer.DBNDD2))

## Indexes for the ranges to color
# c_ind <- color_ind(sort_ir, names(primer.DBNDD2))

## Assign color2 to those indexes
# color = getColor(sort_ir, c_ind, 1, 2)

## Plot the IRAnges
# plotRangesColor(sort_ir, col = color, xlim=c(1500,4100), main="DBNDD2 (AC_000170)")

## SessionInfo
# sessionInfo()

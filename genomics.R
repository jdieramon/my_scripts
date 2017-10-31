# genomics.R -- a set of the most common functions that I use during my genomics analyses
# Copyright (C) 2017 Jose V. Die  <jodiera@upv.es>
# Distributed under terms of the MIT license.

## Define function to get the two variants from a heterozigous SNP

snp_variants <- function(dna) {
    # remove special symbol []
    dna <- gsub("\\[|\\]", "", dna)
    
    # get the position of symbol /
    pos <-  gregexpr2("/", gsub("\\[|\\]", "", dna))[[1]][1]
    
    #var1: base1 + exclude base2
    var1 <-  paste(substr(dna, 1, pos-1), substr(dna, pos+2, nchar(dna)))
    
    # remove empty spaces from the string
    var1 <- gsub(" ", "", var1)
    
    #var2: exclude base1 + base2 
    var2 = paste(substr(dna, 1, pos-2), substr(dna, pos+1, nchar(dna)))
    
    # remove empty spaces from the string
    var2 <- gsub(" ", "", var2)
    
    # report variants
    variants  <-  list(var1 = var1, var2 = var2)
    return(variants)
}

## Usage:
## sequence = "TACCGTCC[G/A]GCCTTC"
## snp_variants(sequence)

# Define a set of functions to convert RPM --> RFC and vice versa

getRFC <- function(rpm, r) {
  # ''' Return the RCF (g-force) '''
  ## rpm, RPM to convert into RFC
  ## r, radius (mm) of centrifuge rotor
  
  g = 1.12 * r * (rpm/1000)**2
  g
}


getRPM <- function(rfc, r) {
  # ''' Return the RPM '''
  ## rfc, RFC (g-force) to convert into RPM
  ## r, radius (mm) of centrifuge rotor
  
    rpm = sqrt(rfc/(1.12 * r))*1000
  rpm
}

## Usage:
## getRFC(1000, 60)
## getRPM(7500, 60)

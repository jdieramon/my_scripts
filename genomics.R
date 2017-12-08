# genomics.R -- a set of the most common functions that I use in my genomics analyses
# Copyright (C) 2017 Jose V. Die  <jodiera@upv.es>
# Distributed under terms of the MIT license.


##------------------------------------------------------------------------------
## Define function to get the two variants from a heterozigous SNP
##------------------------------------------------------------------------------

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


##------------------------------------------------------------------------------
## Define a set of functions to convert RPM --> RFC and vice versa
##------------------------------------------------------------------------------

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


##----------------------------------------------------------------------------------------------
## Define a function to extract the accession id from a table hit downloaded from Blast(NCBI) 
## ---------------------------------------------------------------------------------------------

accs_tidy <- function(blast, acc_type){
    
  ### ''' Take a csv table with Blast hits (downloaded from NCBI) and return a 
  ### character vector containing the accessions id. in a tidy way '''
  
  ## blast, table hit
  ## acc_type, 'XP' for XP-type accessions; 'XM' for XM-type accessions
  
    # Initialize the cleaning by selecting the column that contains the accession 
    blast = hit_table[1:dim(hit_table)[1], 2] # column2
    blast = levels(blast)
    
    ## First, we want to get 1 hit per entry. 
    ## If >1 hit per entry, it splits by ";"
    gi = c()
    for(b in blast) {
        if(grepl(";", b)) {                       # if >1hit/entry, it contains ";"
            val  = (strsplit(b, ";"))             # split the entry by ";"
            gi = c(gi, val[[1]][1], val[[1]][2])  # modify if the entry contains >2hits
        }
        else{gi = c(gi, b)}
        
    }
    
    ## Then, we get just the ID of the hit (index: 18-29)
    accessions = c()
    for (i in gi) {
       
      s = regexec(acc_type, i)[[1]][1]     # this always returns 18
      acc = substr(i, start = s, 30)       # substring from 18-30
      accessions = c(accessions, acc)
    
      }
    
    return(accessions)
}

## Usage
## ncbi = "2MW2NVF9014-Alignment-HitTable.csv"
## hit_table  <- read.csv(ncbi, header= FALSE)
## hits = accs_tidy(hit_table, acc_type = 'XP')

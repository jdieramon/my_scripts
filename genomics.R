# genomics.R -- a set of the most common functions that I use in my genomics analyses
# Copyright (C) 2018 Jose V. Die  <jodiera@upv.es>
# Distributed under terms of the MIT license.


##------------------------------------------------------------------------------
## Get the two variants from a heterozigous SNP
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
## Convert RPM --> RFC and vice versa
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

##------------------------------------------------------------------------------
## Convert nanograms to copy number (standard curve qPCR)
##------------------------------------------------------------------------------
getCopies <- function(ng, amplicon_len) {
    (ng * 6.0221*10**23) / (amplicon_len * 660 * 1*10**9)
    }

# ng = amount of amplicon
# amplicon_len = length of amplicon
# 660 g/mole = avg mass of 1 bp dsDNA


##----------------------------------------------------------------------------------------------
## Extract the accession id from a table hit downloaded from Blast(NCBI) 
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


##--------------------------------------------------------------------------------------------
## Subset a df with unique values from a column containing repetitions
##--------------------------------------------------------------------------------------------

df_wo_duplicates <- function(dataset, ncolDuplicates, ncolTarget) {
    # Subset a data frame with unique values from a column containing repetitions
    
    ## dataset = dataset
    ## ncolDuplicates = number of column containing duplicates 
    ## ncolTarget = number of target column to extract values
    
    
    # Get target values (unique values) to be extracted 
    vals = c()
    target = c()
    for(i in 1:nrow(dataset)) {
        if(!dataset[i,ncolDuplicates] %in% vals){
            target = c(target, as.character(dataset[i, ncolTarget]))
            vals = c(vals, as.character(dataset[i,ncolDuplicates]))
        }
    }
    
    # Subset df based on the target values  
    return(subset.data.frame(dataset, dataset[,ncolTarget] %in% target))
}

## Usage
# df_silly = read.csv("df_silly.csv")
# df_wo_duplicates(df_silly, 2, 1)


##--------------------------------------------------------------------------------------------
## TSS coordinates into a GRanges object
##--------------------------------------------------------------------------------------------

getTSS <- function(inputGR, CDSseqs, bp=200) {
  
  # Output: GRanges object with the TSS coordinates (dependency on strand)
  
  # Input (inputGR) : GRanges object with the mRNA coordinates
  # Input (CDSseqs) : DNAStringSet containing the CDS sequences
  # Input (bp): integer, number of nucleotides extracted from the TSS 
  
  
  # Create an GR image from the inputGR
  gr.tss = inputGR
  
  for(i in seq_along(gr.tss)) {
    
    # If strand +
    if(gr.tss[i] %in% gr.tss[strand(gr.tss) == "+"]) {
      
      bp_cut = bp
      
      text = table4chr(gr.tss[i])
      
      if(i == 10) {bp_cut = 600}         #adapted for my current study
      if(i %in% c(41,42)) {bp_cut = 56}  #adapted for my current study
      
      pattern  = myCDS[[i]][1:bp_cut]
      start(gr.tss[i]) = start(matchPattern(pattern, text, max.mismatch = 0))
      
      bp_cut = bp # take the original bp value
      
    }
    
      
    # If strand -
    if(gr.tss[i] %in% gr.tss[strand(gr.tss) == "-"]) { #do not use 'else' bc some gene might be unmapped (strand = *)
      
      bp_cut = bp
      
      text = table4chr(gr.tss[i])
      
      if(i == 48) {bp_cut = 22}         #adapted for my current study
      
      pattern  = myCDS[[i]][1:bp_cut]
      pattern = reverseComplement(pattern)
      end(gr.tss[i]) = end(matchPattern(pattern, text, max.mismatch = 0))
      bp_cut = bp
    }
      
  }
    gr.tss
}

##--------------------------------------------------------------------------------------------
## Conversion Table for Chromosome
##--------------------------------------------------------------------------------------------

## Function to extract Chr
table4chr <- function(granges) {
  text = c()
  j = paste0("ao",as.numeric(seqnames(granges)))
  if(j == "ao1") {
    text = ao1 }
  if(j == "ao2") {
    text = ao2 }
  if(j == "ao3") {
    text = ao3 }
  if(j == "ao4") {
    text = ao4 }
  if(j == "ao5") {
    text = ao5 }
  if(j == "ao6") {
    text = ao6 }
  if(j == "ao7") {
    text = ao7 }
  if(j == "ao8") {
    text = ao8 }
  if(j == "ao9") {
    text = ao9 }
  
  text
}


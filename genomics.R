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
    pos <-  gregexpr("/", gsub("\\[|\\]", "", dna))[[1]][1]
    
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

# The file is downloaded by clicking selected boxes from the 'Sequences producing significant alignments' section.

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


##----------------------------------------------------------------------------------------------
## Extract the Hit Table ONLY for the best hit per sequence
## ---------------------------------------------------------------------------------------------

# The file is downloaded by clicking 'Download / Hit Table(csv)' on top of the page (Edit and Resubmit, Save Search, ...)

clean_hit <- function(hitFile) {
  # Read Downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  # Chnage colnames
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score", "%positives")
  head(hit)
  
  # Get 1st line per accession (=best macth)
  nline = c()
  xp = c()
  for(l in 1:nrow(hit)) {
    if(! hit$query[l] %in% xp) {
      xp = c(xp, hit$query[l])
      nline = c(nline, l) 
    } 
    
  }
  # Extract just the best hit per accession
  hit[nline,]
  }

## Usage
## Ex. multiple blastp searches (with 2 sequences)
## ncbi = "B1J2TDTZ01R-Alignment-HitTable.csv"
## clean_hit("B1J2TDTZ01R-Alignment-HitTable.csv")



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
## Get the promoter gDNA sequence
##--------------------------------------------------------------------------------------------

get_promoters <- function(inputGR) {
  
  # Input (inputGR) : GRanges object with the promoter coordinates
  # Output: DNAStringSet containing the promoter sequences (bp) 
  
  promoter.seq = c()
  for(i in seq_along(inputGR)) {
    
    text = table4chr(inputGR[i]) # load the corresponding Chr
    
    # If strand +
    if(inputGR[i] %in% inputGR[strand(inputGR) == "+"]) {
      seq = text[start(inputGR[i]):end(inputGR[i])]
      promoter.seq = c(promoter.seq,seq)
    }

    # If strand -
    if(inputGR[i] %in% inputGR[strand(inputGR) == "-"]) {
      seq = text[start(inputGR[i]):end(inputGR[i])]
      seq = reverseComplement(seq)
      promoter.seq = c(promoter.seq, seq)
    }
    
  }
  promoter.seq
  
}

##--------------------------------------------------------------------------------------------
## Conversion Table for Chromosome: for 'get_promoters' and 'getTSS'
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


##--------------------------------------------------------------------------------------------
## PCR insilico: get the coordinates from a dataset containing markers
##--------------------------------------------------------------------------------------------


mapMarkers <- function(dataset, chr, mm) {
  
  ## dataset: df with 2 colums : name, sequence
  ## chr: chromosome to map the sequences from dataset
  ## mm: max. mismatch
  
  multipos = c()
  start = c()
  end = c()
  B = 1:nrow(dataset)
  
  # loop over the markers
  for(i in seq_along(B)) {
    print(paste('SNP',i))                                            ##### print statement                                  
    pattern = DNAString(as.character(snp_variants(dataset[i,2])[1]))
    
    # loop over mismatches: strand +
    for(m in c(head(c(0:mm), n= 3), tail(c(0:mm), n= 3))) { # Check the first 3 mm, and the last 3 mm
    #for(m in 0:mm) {                                       # Check one by one
      print(m)
      mapping = matchPattern(pattern, chr, max.mismatch = m)
      
      # if 1 hit, keep coordinates and break
      if(length(mapping) == 1) {
        #print('1 hit found in strand +')                              ##### print statement
        start = c(start, start(mapping))
        end = c(end, end(mapping))
       break
      }
      # if >1 hit, keep coordinates and break
      if(length(mapping) > 1) {
        multipos <- c(multipos, rep(i, length(mapping)))
        start <- c(start, start(mapping))
        end <- c(end, end(mapping))
        break  
      }
    }
    
    # if no hits --> look in strand - 
    if(length(mapping) == 0) {
      #print(paste('Now searching strand - for',i))                    ##### print statement
      
      #loop over mismatches
      for(m in c(head(c(0:mm), n= 3), tail(c(0:mm), n= 3))) {
        #for(m in 0:mm) {
        print(m)
        pattern = reverseComplement(pattern)
        mapping = matchPattern(pattern, chr, max.mismatch = m)  
        
        # if 1 hit, keep coordinates and break
        if(length(mapping) == 1) {
          #print('1 hit found in strand -')                            ##### print statement
          start = c(start, start(mapping))
          end = c(end, end(mapping))
          break
        }
        
        # if >1 hit, keep coordinates and break
        if(length(mapping) > 1) {
          multipos <- c(multipos, rep(i, length(mapping)))
          start <- c(start, start(mapping))
          end <- c(end, end(mapping))
         break  
        }
      }
      
    }
    # if no hits -- > the sequence does not map to any strand
    if(length(mapping) == 0){
      # print('There is no hit in any strand')                         ##### print statement
      start <- c(start, 0)
      end <- c(end, 0)
    }
    
  }
  
  a = list(start=start, end=end, multiple= multipos)
  a
    
}

## Usage
# coord2 = mapMarkers(snp, ao2, 10)


##--------------------------------------------------------------------------------------------
##  Get the Conserved Protein Domain Family from a protein id
##--------------------------------------------------------------------------------------------


## Dependencies
library(rentrez)


## Define function
get_CD = function(xp) {
  
  ### ''' Take a protein id and returns a vector with the region name of the 
  ### Conserved Protein Domain Family, if applicable '''
    
  # NCBI Eutils through rentrez pckg.
  protein <-  entrez_summary(db="protein", id= xp)
  p_links = entrez_link(dbfrom='protein', id=protein$uid, db='cdd')
  cdd = p_links$links$protein_cdd_concise_2
  
  
  # initializes an empty vector
  acc = c()
  
  # if CDD
  if(! is.null(cdd)) {
    for(i in seq_along(cdd)) {
      region = entrez_summary(db="cdd", id= cdd[i])
      acc = c(acc, region$accession)
    }

  } else {acc = 0}
  
  acc
  
}

## Usage
#xp = c("XP_012572936", "XP_012572937", "YP_001949758")
#for(x in xp) {print(lapply(x, get_CD))}



## -------------------------------------------------------------------------------
## Clean up the n-end characters from a string (basic R , NO stringr)
## -------------------------------------------------------------------------------

clean_end <- function(str, x){
  
  ### str, string vector
  ### x, number of characters
    
  str = substr(str, start = 1, stop = nchar(str)-x)
  str
 }

## Usage: remove the last 2 characters
#srr = c("SRR5927133.121.1", "SRR5927133.121.2", "SRR5927133.245.2", "SRR5927133.1978.1","SRR5927133.1978.2")
# first element
srr[1]
clean_end(srr[1],2)

# for each element in the vector
lapply(srr, function(x) clean_end(x,2))
sapply(srr, function(x) clean_end(x,2), USE.NAMES = FALSE)



## -------------------------------------------------------------------------------
## Extract the end n-chars from a string (basic R , NO stringr)
## -------------------------------------------------------------------------------

get_end <- function(str, x) {
  
  ### str, string vector
  ### x, number of characters
 
  text = as.character(str)
  text = substr(text, start = nchar(text)-x+1, stop = nchar(text))
  text
  
}


## Usage: extract the last 2 characters
## first element
#get_end(srr[1],2)
## for each element in the vector
#sapply(srr, function(x) get_end(x,2), USE.NAMES = FALSE)

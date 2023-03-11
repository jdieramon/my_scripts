# genomics.R -- a set of the most common functions that I use in my genomics analyses
# Copyright (C) 2019 Jose V. Die  <jodiera@upv.es>
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

# ng = amount of amplicon
# amplicon_len = length of amplicon
# 660 g/mole = avg mass of 1 bp dsDNA


##----------------------------------------------------------------------------------------------
## Extract the accession id from a table hit downloaded from blast (NCBI) 
## ---------------------------------------------------------------------------------------------

# The file is downloaded by clicking selected boxes from the 'Sequences producing significant alignments' section.

blast_accs_tidy <- function(blast, acc_type){
    
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
## hits = blast_accs_tidy(hit_table, acc_type = 'XP')


##----------------------------------------------------------------------------------------------
## Extract the Hit Table ONLY for the best hit per sequence
## ---------------------------------------------------------------------------------------------

# This function assumes that the query is an unknown sequence. If it is an id. (GenBank) the best hit in BLAST will be itself and 
# this function is useless. In that case, use the function 'best_homolog'  
# The file is downloaded by clicking 'Download / Hit Table(csv)' on top of the page (Edit and Resubmit, Save Search, ...)

blast_best_hit <- function(hitFile) {
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
## blast_best_hit("B1J2TDTZ01R-Alignment-HitTable.csv")



##----------------------------------------------------------------------------------------------
## Extract the Hit Table ONLY for the best hit per sequence
## ---------------------------------------------------------------------------------------------

# The file is downloaded by clicking 'Download / Hit Table(csv)' on top of the page (Edit and Resubmit, Save Search, ...)

blast_best_homolog <- function(hitFile) {
  
  ### ''' Take a csv table with Blast hits and return a 
  ### subset containing only the best hit for each query '''
  
  # Read Downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # Change colnames
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score", "%positives")
  
  res <- as_tibble(hit) %>% 
    filter(query != subject) %>% 
    group_by(query) %>% 
    slice(1) %>% 
    ungroup()
  
  res
  
}

## Usage
## Ex. multiple blastp searches (with 2 sequences)
## blast_best_homolog("B1J2TDTZ01R-Alignment-HitTable.csv")



##----------------------------------------------------------------------------------------------
## Extract the Hit Table filtering by identity and coverage 
## ---------------------------------------------------------------------------------------------

# The file is downloaded by clicking 'Download / Hit Table(csv)' on top of the page (Edit and Resubmit, Save Search, ...)
# Coverage of the longer aa sequence 


get_xp_length <- function(xp_id){
  ### ''' Estimate protein length (aa)
  
  xpinfo <- entrez_summary(db = "protein", id = xp_id)
  xpinfo$slen
  
}


blast_homologs <- function(hitFile, ident, cover) {
  
  ### ''' Take a csv table with Blast hits and return a 
  ### subset filteres by %Identity & %Coverage (longer sequence)
  
  # read downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # change colnames
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score", "%positives")
  
  # make a table_length for queries & subjects 
  queries <- hit %>% 
    pull(query) %>% 
    unique()
  
  tlengths_query = tibble(query = queries, len.query = sapply(query, function(i) get_xp_length(i)))
  
  subjects <- hit %>% 
    pull(subject) %>% 
    unique()
  
  tlengths_sub = tibble(subject = subjects, len.sub = sapply(subject, function(i) get_xp_length(i)))
  
  
  # combine hit_table with queries & subjects lengths 
  # estimate coverage (longer sequence)
  # apply filters : identity * coverage
  hit %>% 
    tibble() %>% 
    select(-c(align_length, bit_score, "%positives", mismatches, gap_opens)) %>% 
    right_join(tlengths_query, by = "query") %>% 
    right_join(tlengths_sub,   by = "subject") %>% 
    mutate(cov.query = (q.end - q.start + 1)/len.query, 
           cov.sub = (s.end - s.start + 1)/len.sub) %>% 
    mutate(cov_long = case_when(len.query > len.sub ~ cov.query*100, 
                                len.query < len.sub ~ cov.sub*100, 
                                TRUE ~ cov.query*100)) %>% 
    select(-c(q.start, q.end, s.start, s.end, len.query, len.sub, cov.query, cov.sub)) %>% 
    filter(identity >= ident, cov_long >= cover) %>% 
    select(query, subject, identity, cov_long, evalue) %>% 
    arrange(identity)
  
  
}

## Usage
## blast_homologs("B1J2TDTZ01R-Alignment-HitTable.csv", ident = 55, cover = 85)


# BLAST Align two or more sequences -------------------------------------------------                                                             
# Sequence mRNA : 5UTR+CDS+3UTR (it may be important to identify the longer sequence). 
blast_homologs_mRNA <- function(hitFile, ident, cover){
  
  ### ''' Take a csv table with Blast hits "Align two or more sequences"
  ### return subset filtered by %Identity & %Coverage (longer sequence)
  # Ex. Duplication Analysis
  # mRNA = 5UTR+CDS+3UTR
  
  # read downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # change colnames
  # en blast-2 sequences contra ellas mismas no hay 13 columnas sino 12 (-"%positives")
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score")
  
  # make a table_length for queries & subjects 
  tlengths <- hit %>% 
    filter(query == subject) %>% 
    select(query, align_length) %>% 
    arrange(query, desc(align_length)) %>% 
    distinct(query, .keep_all = TRUE) %>% 
    rename(length = align_length)
  
  
  # combine hit_table with queries & subjects lengths 
  # estimate coverage (longer sequence)
  as_tibble(hit) %>% 
    select(-c(align_length, bit_score, mismatches, gap_opens)) %>% 
    filter(query != subject) %>% 
    right_join(tlengths, by = "query") %>% 
    rename(len.query = length) %>% 
    right_join(tlengths %>% rename(subject = query), by = "subject") %>% 
    rename(len.sub = length) %>% 
    mutate(cov.query = (q.end - q.start + 1)/len.query, 
           cov.sub = (s.end - s.start + 1)/len.sub) %>% 
    mutate(cov_long = case_when(len.query > len.sub ~ cov.query*100, 
                                len.query < len.sub ~ cov.sub*100, 
                                TRUE ~ cov.query*100)) %>% 
    select(-c(q.start, q.end, s.start, s.end, len.query, len.sub, cov.query, cov.sub)) %>% 
    filter(identity >= ident, cov_long >= cover) %>% 
    select(query, subject, identity, cov_long, evalue) %>% 
    arrange(identity)
  
}

## Usage
#hitFile = "dat/0SR2C6E2114-Alignment-HitTable.csv"
#blast_homologs_mRNA(hitFile, ident = 60, cover = 50)    

# Sequence CDS (it may be important to identify the longer sequence). 
extract_CDS <- function(xm) {
  # ''' extract the CDS (start/end) coordinates

efetch = rentrez::entrez_fetch(db= 'nuccore', id = xm, rettype = 'gp')
  listName = strsplit(efetch, "\n")
  for(i in seq(listName[[1]])) {
    val <- listName[[1]][i]
    #remove whitespaces from the string
    val = gsub(" ", "", val)
    
    # check for feature
    if(substr(val, 1, 3) == "CDS") {
      #remove characters "CDS" from the string
      val = gsub("CDS", "", val)
      #remove special symbols "..." from the string
      val = strsplit(val, "[..]")
      #elements 1-3 of the list contain the start/stop coordinates
      start = as.numeric(val[[1]][1])
      stop = as.numeric(val[[1]][3])
      cds = list(startCDS = start, stopCDS = stop)
      
      # return for downstream analysis
      return(cds)
    }
    }
  }

getCDS <- function(xmsIds) {
  # ''' return DNAstringobject with the CDS sequence
  mycds = c()
  mycds = sapply(xmsIds, function(i){
    
    coord = extract_CDS(i)
    cds = rentrez::entrez_fetch(db="nucleotide", id=i, rettype="fasta",
                                seq_start = coord$startCDS, seq_stop = coord$stopCDS)
    cds_tidy = strsplit(cds, "\n")
    
    cds_tidy <- as.character(paste0(cds_tidy[[1]][2:length(cds_tidy[[1]])], collapse = ""))
    mycds = c(mycds, cds_tidy)
  })
  
  # Result 
  DNAStringSet(mycds)
  
}
                                                             
##Usage
# FIRST : Write sequences to a fasta file
#writeXStringSet(getCDS(xmsIds), filepath="dat/CDS4dup.fasta", format="fasta")  

# SECOND : Load your fasta file 'Align two or more sequences' into BLAST and download the "Hit table(csv)" 
#hitFile = "dat/0SW43G15114-Alignment-HitTable.csv"

# THIRD : Analyze BLAST output
#blast_homologs_CDS(hitFile, ident = 60, cover = 50) 



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
#lapply(xp, function(x) get_CD(x))
#sapply(xp, function(x) get_CD(x))




## -------------------------------------------------------------------------------
## Clean up the n-end characters from a string (basic R , NO stringr)
## -------------------------------------------------------------------------------

clean_end <- function(str, x){
  
  ### str, string vector
  ### x, number of characters
    
  str = as.character(str)
  str = substr(str, start = 1, stop = nchar(str)-x)
  str
 }

## Usage: remove the last 2 characters
#srr = c("SRR5927133.121.1", "SRR5927133.121.2", "SRR5927133.245.2", "SRR5927133.1978.1","SRR5927133.1978.2")
# first element
#srr[1]
#clean_end(srr[1],2)


# for each element in the vector
#lapply(srr, function(x) clean_end(x,2))
#sapply(srr, function(x) clean_end(x,2), USE.NAMES = FALSE)



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


## -------------------------------------------------------------------------------
## Extract the characters after a period (stringr)
## -------------------------------------------------------------------------------

get_after_period <- function(my_vector) {
  
  # Return a string vector without the characters 
  # before a period (excluding the peiod)
  
  # my_vector, a string vector 
  
  str_sub(my_vector, str_locate(my_vector, "\\.")[,1]+1)
  
}

## Usage
#my_vector <-  c('foobar.barfoo', 'amazing.point')
#get_after_period(my_vector)



## -------------------------------------------------------------------------------
## Extract the species name from a sequence description (stringr)
## -------------------------------------------------------------------------------

# Dependency
#library(stringr)

get_spp <-
function(description) {
  
  spp <- str_sub(description, 
                 start = str_locate(description, "\\[")[1]+1,
                 end = str_locate(description, "\\]")[2]-1)
  spp
}

## Usage: 
#desc <-  c("PREDICTED: auxin response factor 19-like isoform X1 [Glycine max]")
#get_spp(desc)  


## -------------------------------------------------------------------------------
# Taxonomy identifiers <-->  scientific names (NCBI)
## -------------------------------------------------------------------------------

# Dependencies
# library(rentrez)

# Define some functions : 
get_tax_id <- function(scientific_name) {
  myterm = paste0(scientific_name,"[SCIN]")
  esearch = entrez_search(db = "taxonomy", term = myterm)
  esearch$ids
  }
  
get_scientific_name <- function(taxid) {
  esumm <- entrez_summary(db="taxonomy", id=taxid)
  esumm$scientificname
  }

## Usage
# From species name -> taxonomy id :
#mycicer <- c("Cicer macracanthum", "Cicer canariense", "Cicer chorassanicum", "Cicer cuneatum", "Cicer judaicum", 
#        "Cicer yamashitae", "Cicer bijugum", "Cicer reticulatum", "Cicer echinospermum", "Cicer pinnatifidum", 
#        "Cicer arietinum")

#mytaxids = sapply(mycicer, function(x) get_tax_id(x), USE.NAMES = FALSE)

# From taxonomy id -> scientific name :
# mynames = sapply(mytaxids, function(x) get_scientific_name(x), USE.NAMES = FALSE)


## -------------------------------------------------------------------------------
# NCBI Accessions characterization 
## -------------------------------------------------------------------------------

characterizeTable <- function(targets) {
  
  # targets, character vector with XM protein accession ids
  
  ## Initialize the vectors containing the features that we need for each sequence
  LOC    <-  vector("character") 
  Chr    <-  vector("character")
  chr_s  <-  vector("integer")
  chr_e  <-  vector("integer")
  exon   <-  vector("integer") 
  AA     <-  vector("integer")
  mol_wt <-  vector("integer")
  xm     <-  vector("integer")
  xp     <-  vector("integer")
  
  ## Extract feature info from NCBI
  for(xi in as.character(targets)) {
    xp = c(xp, xi)
    xm = c(xm, getXM(xi))
    
    xpinfo <- entrez_summary(db = "protein", id = xi)
    AA = c(AA, xpinfo$slen)
    
    protein_gb <- entrez_fetch(db = "protein", id = xi, rettype = "gp")
    mol_wt = c(mol_wt, extract_mol.wt_from_xp(protein_gb)/1000) # kDa
    
    xplink <- entrez_link(dbfrom = "protein", id = xi, db = "gene")
    genesummary = entrez_summary(db = "gene", id = xplink$links[1])
    LOC = c(LOC, paste0("LOC",genesummary$uid))
    Chr = c(Chr, genesummary$chromosome)
    
    # defensive programming
    # some XP may not have chrstart, chrstop or exon count
    if(length(genesummary$genomicinfo$chrstart) != 0) {
      chr_s = c(chr_s, genesummary$genomicinfo$chrstart)
      chr_e = c(chr_e, genesummary$genomicinfo$chrstop)
      exon = c(exon, genesummary$genomicinfo$exoncount)  
      
    } else {
      chr_s = c(chr_s, 0)
      chr_e = c(chr_e, 0)
      exon = c(exon, 0)  
      
    }
  }
  
  ## Build the dataset 
  t1 <- data.frame(LOC, XP = xp, Chr, chr_s, chr_e, AA, mol_wt, XM = xm, exon)
  
  ## Add strand
  t1 <- t1 %>% mutate(Strand = ifelse(chr_e > chr_s, "+", "-"))
  
  ## Sort data set by Chr + start coordinate
  t1 = t1 %>% arrange(Chr, chr_s)
  
  t1
  
}
## -------------------------------------------------------------------------------
# Re-order a vector based on values shown in another vector
## -------------------------------------------------------------------------------

make_index <- function(vector1, vector2) {
  # re-order vector2 based on the order of vector1
  # returns an index 
  
  sapply(vector1, function(x) which(x == as.character(vector2)), USE.NAMES = F)
  
}

#Usage: 
#Equivalent to base function : `match`
#pets <- c("dog", "cat", "bird", "deer", "chicken")
#wild <- c("chicken", "cat",  "dog", "deer", "bird") 
#pets[make_index(wild, pets)]

#pets[match(wild, pets)]

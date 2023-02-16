# MADS-box.R -- a set of common functions for analysis of MADS-box family in Ca
# Copyright (C) 2023 Jose V. Die  <jose.die@uco.es>
# Distributed under terms of the MIT license.


best_homolog <- function(hitFile) {
  
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
## best_homolog("B1J2TDTZ01R-Alignment-HitTable.csv"



#library(dplyr)
#library(refseqR)
#library(rentrez)

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

## Usage
##targets <- c("XP_004498126.1")
##targets <- c("XP_004513719.1", "XP_004485955.1", "XP_004498126.1", "XP_004498126.1")
##tdat <- characterizeTable(targets)

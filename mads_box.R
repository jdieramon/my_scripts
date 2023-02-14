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

## Dependencies
library(rentrez)
library(stringr)

## ----------------------------------------------------------------------------------------------
## Read the prot. ids as data frame
ids <- read.csv("prot_ids", header = FALSE)

## Convert into vector
ids = ids[,1]


## Loop over protein ids

genome = c()

  
for(i in seq_along(ids)) {
  
  target = ids[i]
  
  # NCBI E-utilities
  foo = entrez_link(dbfrom='protein', id=target, db='all')
  
  # Link to Gene database
  faa = entrez_summary(db = "gene", id = foo$links$protein_gene)
  
  # From genomic info, extract the genome
  genome[i] = faa$genomicinfo[,2]
  
  }


## Build df containing columns: prot.ids, genome
typhi = data.frame(ids, genome) 

## Export to csv file
write.csv(typhi, file = "summary")



## ----------------------------------------------------------------------------------------------
## 'getSize' takes a vector containing NCBI SRA accessions and returns and vector 
## with the run size (Mb) 


getSize <- function(ids) {
  
  sizes <- c()
  
  for(i in ids) {
    
    term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    
    ## extract range defined by "total bases" and "load_done"
    size = str_sub(x, start =str_locate(x, "total_bases=")[2]+2, 
                   end =str_locate(x, "load_done")[1]-3)
    
    size_mega = as.numeric(size)/1e6
    
    sizes <- c(sizes, size_mega)
    
  }
  
  sizes
  
}

## Usage 
# ids <- c("DRR071071", "ERR1912953")
# getSize(ids)



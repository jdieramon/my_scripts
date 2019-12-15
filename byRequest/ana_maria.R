# Copyright (c) 2019 Jose V. Die

# Dependencies
library(dplyr)
library(tidyr)
library(tibble)
library(rentrez)
library(XML)

# ------------------------------------------------------------------------
##  GO from Uniprot ids 
# ------------------------------------------------------------------------

# Retrieve/Id mapping : https://www.uniprot.org/uploadlists/
# Add columns GO : BP, MF, CC
# Download all / Format : Tab-separated / uncompressed 

# Change file name 
system("mv ~/Downloads/*.tab ~/Downloads/uniprot.tab")

# Read file 
uniprot <- read.delim("data/uniprot.tab")

# Check file
head(uniprot)
names(uniprot)

uniprot <- as_tibble(uniprot %>% 
                       select(Entry, Protein.names, Gene.names, 
                              Organism, Gene.ontology..molecular.function., 
                              Gene.ontology..cellular.component., Gene.ontology..biological.process.) %>% 
                       rename(MF = Gene.ontology..molecular.function., 
                              BP = Gene.ontology..biological.process., 
                              CC = Gene.ontology..cellular.component.))

uniprot  

# ------------------------------------------------------------------------
##  TAXONOMY 
# ------------------------------------------------------------------------
# Read csv file with protein ids
prot_ids <- read.csv("data/prot_ids.txt", header = F, stringsAsFactors = F)

# Check out the df
prot_ids %>% count(V1)

# Make a df with unique prot ids. 
prot_ids %>% distinct(V1)

# Extract the unique ids into a vector
ids <- prot_ids %>% 
  distinct(V1) %>% 
  pull(V1)


## Programatic access to UniProt
# Mapping database identifiers
uri <- 'http://www.uniprot.org/uniprot/?query='
idStr <- paste(ids, collapse="+or+")
format <- '&format=tab'
fullUri <- paste0(uri,idStr,format)
dat <- read.delim(fullUri)

# Check UniProt results
as_tibble(dat)
glimpse(dat)


# Programatic access to ENTREZ : Taxonomy database
# -----------------------------------------------------
# Make a vector with taxonomy ids 
orgs <- levels(dat$Organism)

taxids <- c()
for(org in orgs) {
  foo = entrez_search(db = "taxonomy", term = paste0(org,"[SCIN]"))
  taxids <- c(taxids, foo$ids)
}


# Check equality in vector length 
length(orgs) == length(taxids)


# Programatic access to ENTREZ : xml file
source("https://raw.githubusercontent.com/jdieramon/my_scripts/master/byRequest/functions.R")
vals <- get_taxonomy(taxids)

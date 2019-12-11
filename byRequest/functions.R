# Dependencies
library(dplyr)
library(tidyr)
library(tibble)
library(rentrez)
library(XML)


# Programatic access to ENTREZ : xml file

get_taxonomy <- function(taxids){
  # Take a vector with protein ids (uniprot)
  # Return a tibble with taxonomic data
  
  # make an empty dataframe (outcome)
  outcome <- data.frame()
  
  # loop over the taxids
  for(id in taxids) {
    result = entrez_fetch(db = "taxonomy", id = id,
                          rettype="xml", parsed=TRUE)  
    
    # make an empty df (tmp)
    mydf = tibble()
    
    # extract the root node from the xml file
    rootnode = xmlRoot(result)
    
    mydf = cbind(setNames(xmlToDataFrame(getNodeSet(result, "//Rank")), "Rank"),
                 setNames(xmlToDataFrame(getNodeSet(result, "//ScientificName")), "Taxon"))
    
    outcome = rbind(outcome, mydf %>% 
                      filter(Rank %in% c("superkingdom", "phylum", "class", 
                                         "order", "family")) %>% 
                      mutate(Spp = as.character(mydf$Taxon[1])) %>% 
                      select(Spp, Rank, Taxon)  %>% 
                      spread(Rank, Taxon) %>% 
                      select(Spp, superkingdom, phylum, class, order, family)) 
    
  }
  
  
  as_tibble(outcome)
  
  
}
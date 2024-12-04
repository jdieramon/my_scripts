# Dependencies 
library(dplyr)
library(stringr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Analize the output file from OmicsBox 
# -----------------------------------------------------------------------------

# OMICSBOX / EXPORT 



# --- Analysis GENERIC EXPORT --------------------------------------------------

genexport <- read.delim("dat/generic_export")

# Sequence.Name Annotation.GO.ID            Annotation.GO.Term Annotation.GO.Category
# 1   XP_004495769.1       GO:0000166       nucleotide binding     Molecular Function
# 2   XP_004495769.1       GO:0005515       protein binding     Molecular Function
# 3   XP_004495769.1       GO:0005634       nucleus     Cellular Component


# Quick view
genexport %>% glimpse()

# Filter out sequences without GO terms 
genexport <- genexport %>% 
  filter(Annotation.GO.ID %in% str_subset(Annotation.GO.ID, "")) %>% 
  as_tibble() 


# Number of  GO terms 
genexport %>% count(Annotation.GO.Category)

# Sequences with GO terms from 3, 2, or 1 Functional Category  

FC_by_seq <- function(i) {
    #''' i, number of GOs
  
  genexport %>% 
    count(Sequence.Name, Annotation.GO.Category) %>% 
    arrange(Sequence.Name) %>% 
    count(Sequence.Name) %>% 
    filter(n == i)
  
}

lapply(1:3, function(i) FC_by_seq(i))



# Number of BP terms per sequence 
genexport %>% 
  filter(Annotation.GO.Category == "Biological Process") %>% 
  count(Annotation.GO.Category, Sequence.Name, sort = T) %>% 
  select(Sequence.Name, n)

# Number of MF terms per sequence 
genexport %>% 
  filter(Annotation.GO.Category == "Molecular Function") %>% 
  count(Annotation.GO.Category, Sequence.Name, sort = T) %>% 
  select(Sequence.Name, n)

# Number of CC terms per sequence 
genexport %>% 
  filter(Annotation.GO.Category == "Cellular Component") %>% 
  count(Annotation.GO.Category, Sequence.Name, sort = T) %>% 
  select(Sequence.Name, n)




# --- Analysis SLIM ----------------------------------------------------------

# Check if any BP, CC, MF term is not available in the 'slim.rda' object 

# Load slim classification
load("dat/slim.rda")

# check available categories per Ontology
bp %>% count(slim_class, sort = T)
mf %>% count(slim_class, sort = T)
cc %>% count(slim_class, sort = T)

genexport_cc <- genexport %>% 
  rename("slim_GO_ID" = Annotation.GO.ID) %>% 
  inner_join(cc, by = "slim_GO_ID")

genexport_mf <- genexport %>% 
  rename("slim_GO_ID" = Annotation.GO.ID) %>% 
  inner_join(mf, by = "slim_GO_ID")

genexport_bp <- genexport %>% 
  rename("slim_GO_ID" = Annotation.GO.ID) %>% 
  inner_join(bp, by = "slim_GO_ID")


# Check if you need to update the `slim´ object.
add_slim <- function() {
  
  slims <- genexport_bp %>% 
    bind_rows(genexport_mf, genexport_cc)
  
  slims_uniq <- slims %>% pull(slim_GO_ID) %>% unique()
  
  omics <- genexport %>% pull(Annotation.GO.ID) %>% unique()
  
  
  if(length(slims_uniq) < length(omics)) {
    cat(paste0("The OmicsBox annotation has ", length(omics), 
" GO terms but the `slim` object has  ", length(slims_uniq), " GO terms.
        \nYou need to update the `slim` object."))
    
    cat("\n")
    omics[!omics %in% slims_uniq]
    
    
  } else {cat("No need to update the `slim´ object.")}
}

add_slim()

# --- how to update the slim object : 
# Example : "GO:0009507"
## goto : https://www.ebi.ac.uk/QuickGO/search/GO:0042221
# get the class : 'cc'
# get the name : chloroplast 
# get the class (from the ancestor chart; hierarchy icon)  : plastid / chloroplast
# 
# cc <- cc %>% bind_rows(data.frame(GO_type = c("CC"), 
#              slim_GO_ID = c("GO:0009507"), 
#              slim_GO_Name = c("chloroplast"), 
#              slim_class = c("chloroplast")))
# 
# bp <- bp %>% bind_rows(data.frame(GO_type = c("BP"), 
#            slim_GO_ID = c("GO:0065009"), 
#            slim_GO_Name = c("regulation of molecular function"), 
#            slim_class = c("other biological processes")))
# 
# mf <- mf %>% bind_rows(data.frame(GO_type = c("MF"), 
#            slim_GO_ID = c("GO:0038023"), 
#            slim_GO_Name = c("signaling receptor activity"), 
#            slim_class = c("other molecular functions")))
# 
# # Save updated slim classification object
# save(bp, mf, cc, file = "slim.rda")




#Annotationanalysis 
genexport_mf  %>% 
  filter(slim_class == "transcription factor activity") %>% 
  select(-c(GO_type, slim_class))

genexport_bp %>% 
  filter(slim_class %in% str_subset(slim_class, "response"))

genexport_bp %>% 
  count(slim_class, sort = T) %>% 
  ggplot(aes(reorder(slim_class, n), n)) + 
  geom_col(fill = "steelblue")  + 
  coord_flip() + xlab("") + ylab("")



## analisis for Agbase output file
# used in : NILs


# Load Ag file
nils <- read.delim("dat/AgBase.txt")

as_tibble(nils)

# Make ontology datasets
nils_cc <- as_tibble(nils) %>% 
  filter(GO_Type == "C")

nils_cc %>% 
count(Slim_GO_Name, sort = T)


nils_bp <- as_tibble(nils) %>% 
  filter(GO_Type == "P")

nils_mf <- as_tibble(nils) %>% 
  filter(GO_Type == "F")



# Make GO slim plots
nils_bp %>% 
  count(Slim_GO_Name, sort = T) %>% 
  ggplot(aes(reorder(Slim_GO_Name, n), n)) + 
  geom_col()  + 
  coord_flip()

nils_mf %>% 
  count(Slim_GO_Name, sort = T) %>% 
  ggplot(aes(reorder(Slim_GO_Name, n), n)) + 
  geom_col()  + 
  coord_flip()


# check specific categories
nils_bp %>% filter(Slim_GO_Name == "response to light stimulus") 
nils_bp %>% filter(Slim_GO_Name == "reproduction") 
nils_bp %>% filter(Slim_GO_Name == "circadian rhythm") 


# general descriptions
nils_bp %>% 
  filter(Slim_GO_Name == "response to light stimulus") %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "reproduction")) %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "circadian rhythm")) %>% 
  count(Input_Accession)

nils_bp %>% 
  filter(Slim_GO_Name == "response to light stimulus") %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "reproduction")) %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "circadian rhythm")) %>% 
  count(Input_GO_Name, sort = T)



# match to genome annotation
nils_ds <- nils_bp %>% 
  filter(Slim_GO_Name == "response to light stimulus") %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "reproduction")) %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "circadian rhythm"))
  
xms = nils_ds %>% pull(Input_Accession) %>% unique()


getXMdescription <- function(xm) {
  transcript_elink = rentrez::entrez_link(dbfrom = "nuccore", 
                                          id = xm, 
                                          db= "gene")
  gene_id = transcript_elink$links$nuccore_gene
  gene = rentrez::entrez_summary(db = "gene", id = gene_id)
  gene$description
  
}


desc = sapply(xms, function(i) getXMdescription(i), USE.NAMES = F)

ds <- data_frame(Input_Accession = xms, desc = desc) %>% tibble()

nils_ds %>% inner_join(ds, by = "Input_Accession") %>% 
  transmute(Input_Accession, desc, Slim_GO_Name, GO_Name = Input_GO_Name)


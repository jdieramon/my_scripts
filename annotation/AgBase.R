## analisis for Agbase output file
# used in : NILs


# Load Ag file
nils <- read.delim("data/AgBase.txt")

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
nils_bp %>% 
  filter(Slim_GO_Name == "response to light stimulus") %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "reproduction")) %>% 
  bind_rows(nils_bp %>% filter(Slim_GO_Name == "circadian rhythm"))
  




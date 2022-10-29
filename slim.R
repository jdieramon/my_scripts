## Read the ontologies 'bp' , 'mf', 'cc' and update it with new GO slims

# Load slim class object
load("data/slim.rda")


# New GO slims to add into 'bp' , 'mf', 'cc' 
my_go <- "GO:0036211" "GO:0009416" "GO:0140110" "GO:0042221" "GO:0065009"


# Example : "GO:0042221"
## goto : https://www.ebi.ac.uk/QuickGO/search/GO:0042221
# get the class : 'bp'
# get the name : response to chemical 
# get the class (from the ancestor chart; hierarchy icon)  : response to stimulus


# update the objects : 'bp' , 'mf', 'cc' 
bp <- bp %>% 
  bind_rows(
  data.frame(GO_type = c("BP", "BP", "BP"), 
             slim_GO_ID = c("GO:0036211", "GO:0009416", "GO:0042221"), 
             slim_GO_Name = c("protein modification process", 
                              "response to light stimulus", "response to chemical"), 
             slim_class = c("protein metabolism", "response to abiotic or biotic stimulus", 
                            "response to abiotic or biotic stimulus")))



mf <- mf %>% 
  bind_rows(c(GO_type = "MF", 
              slim_GO_ID = "GO:0140110", 
              slim_GO_Name = "transcription regulator activity", 
              slim_class = "other enzyme activity "))

# Save updated slim classification object
save(bp, mf, cc, file = "slim.rda")

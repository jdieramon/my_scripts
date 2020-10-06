# Dependencies
library(dplyr)
library(tibble)
library(stringr)


# Define a function to change "." -> ","            
make_comma <- function(string) {
  str_replace(string, pattern = "\\.", replacement = ",")
  
}

# Load data 
myfile <- read.table("Downloads/scan_data2", sep = " ", stringsAsFactors=FALSE)

# Check data
glimpse(myfile)

myfile %>% 
  as_tibble() 
  

# Change "." to ","
myfile %>% 
  as_tibble() %>% 
  mutate_if(is.double, make_comma)
            

# Apply changes
myfile <- myfile %>% 
  as_tibble() %>% 
  mutate_if(is.double, make_comma)

# Write csv into a file 
write.csv(myfile, file = "Desktop/mary.csv", row.names = FALSE)

            


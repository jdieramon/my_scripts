# Dependencies
library(dplyr)
library(stringr)

# Read csv file
prots <- read.csv("byRequest/data/prot.csv")

# Make an index for rows with >1 ID
ids <- prots %>% pull(ID)
idx <- str_which(string = ids, pattern = ";")

# Subset the dataframe to tidy 
df <- prots %>% slice(idx)

# Make some tidy (in an object list)
mydf = list()
for(i in seq(nrow(df))) {
  
  p <- df %>% slice(i) %>% pull(ID)
  df_tmp <- df %>% filter(ID == p) %>% select(-ID)
  ids_split = str_split(p, ";")
  n = length(ids_split[[1]])
  list_tmp = vector("list", n)
  for(c in seq_along(df_tmp)) {
    list_tmp[[c]] = rep(df_tmp[[c]], n)  
  }
  
  mydf[[i]] = tibble(tmp = rep(NA, n))
  for(c in seq_along(list_tmp)) {
    my_col = tibble(list_tmp[[c]])  
    mydf[[i]] = mydf[[i]] %>% bind_cols(my_col)
  }
  mydf[[i]] = mydf[[i]] %>% select(-tmp)
  names(mydf[[i]]) = names(df_tmp)
  mydf[[i]] <- mydf[[i]] %>% 
    bind_cols(ID = ids_split[[1]]) %>% 
    select(ID, everything())
  
}

# Extract the list elements into a dataframe 
res = tibble()
for(i in seq_along(mydf)) {
  res = res %>% bind_rows(mydf[[i]])  
}

# Combine with raws == 1ID from original data 
res <- res %>% 
  bind_rows(prots %>% slice(-idx))

res

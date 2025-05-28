---
title: "structure"
author: "Jose V. Die"
date: "5/23/2025"
output: 
  html_document: 
    keep_md: yes
---



## Rationale

The code transforms a typical matrix with fragment analysis data into an input format for the program [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html), the free software package for using multi-locus genotype data to investigate population structure. 



```r
library(polysat)
library(dplyr)
library(tidyr)
library(stringr)
```

## Data pre-processing
Load data and check content.   

```r
dat = read.csv2("data/Alnus_madresvs brinzales.csv")
#dat = read.csv2("data/Alnus_madres vs brinzales_2x_3x_4x_homocigotos completos.csv")
dat %>% tibble()
```

```
## # A tibble: 276 × 67
##    SAMPLE POPULATION Ploidía  AG01     X   X.1   X.2  AG05   X.3   X.4   X.5
##    <chr>  <chr>      <chr>   <int> <int> <int> <int> <int> <int> <int> <int>
##  1 T-481  Madres     2x        134   146    NA    NA   148   159    NA    NA
##  2 T-482  Madres     4x        132    NA    NA    NA   141   155   157    NA
##  3 T-483  Madres     4x        132    NA    NA    NA   141   146   157    NA
##  4 T-484  Madres     4x        132    NA    NA    NA   146   151   155    NA
##  5 T-485  Madres     4x        132    NA    NA    NA   141   146    NA    NA
##  6 T-486  Madres     4x        132   134    NA    NA   141   155    NA    NA
##  7 T-487  Madres     4x        128   132   166    NA   141   148   155    NA
##  8 T-488  Madres     4x        132    NA    NA    NA   141   148   155    NA
##  9 T-489  Madres     2x        134   135    NA    NA   151   159    NA    NA
## 10 T-490  Madres     4x        126   132    NA    NA   148   155   171    NA
## # ℹ 266 more rows
## # ℹ 56 more variables: AG09 <int>, X.6 <int>, X.7 <int>, X.8 <int>, AG10 <int>,
## #   X.9 <int>, X.10 <int>, X.11 <lgl>, AG13 <int>, X.12 <int>, X.13 <int>,
## #   X.14 <int>, AG23 <int>, X.15 <int>, X.16 <int>, X.17 <lgl>, AG25 <int>,
## #   X.18 <int>, X.19 <int>, X.20 <int>, AG27 <int>, X.21 <int>, X.22 <int>,
## #   X.23 <lgl>, AG30 <int>, X.24 <int>, X.25 <int>, X.26 <int>, AG35 <int>,
## #   X.27 <int>, X.28 <int>, X.29 <int>, A2 <int>, X.30 <int>, X.31 <lgl>, …
```

In the input file (.csv), two things can happen: either the last column contains marker data or it does not. Sometimes, an extra column is added at the end of the .csv file that does not contain marker data. It is necessary to check whether the last column contains any information (in which case nothing is done) or if it is empty (in which case it should be removed).


```r
dat %>% glimpse()
```

```
## Rows: 276
## Columns: 67
## $ SAMPLE     <chr> "T-481", "T-482", "T-483", "T-484", "T-485", "T-486", "T-48…
## $ POPULATION <chr> "Madres", "Madres", "Madres", "Madres", "Madres", "Madres",…
## $ Ploidía    <chr> "2x", "4x", "4x", "4x", "4x", "4x", "4x", "4x", "2x", "4x",…
## $ AG01       <int> 134, 132, 132, 132, 132, 132, 128, 132, 134, 126, 126, 132,…
## $ X          <int> 146, NA, NA, NA, NA, 134, 132, NA, 135, 132, 128, 166, 142,…
## $ X.1        <int> NA, NA, NA, NA, NA, NA, 166, NA, NA, NA, 132, NA, 166, 132,…
## $ X.2        <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 166, NA, NA, NA, NA…
## $ AG05       <int> 148, 141, 141, 146, 141, 141, 141, 141, 151, 148, 141, 141,…
## $ X.3        <int> 159, 155, 146, 151, 146, 155, 148, 148, 159, 155, 144, 146,…
## $ X.4        <int> NA, 157, 157, 155, NA, NA, 155, 155, NA, 171, 155, 155, NA,…
## $ X.5        <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ AG09       <int> 250, 246, 244, 244, 246, 244, 244, 246, 250, 246, 246, 242,…
## $ X.6        <int> 254, 248, 246, 248, NA, 246, 246, 248, 250, NA, 248, 248, 2…
## $ X.7        <int> NA, 250, 248, 250, NA, 262, 250, NA, NA, NA, 250, NA, NA, N…
## $ X.8        <int> NA, NA, NA, 256, NA, NA, NA, NA, NA, NA, 258, NA, NA, NA, N…
## $ AG10       <int> 216, 218, 218, 216, 218, 216, 218, 216, 216, 218, 216, 216,…
## $ X.9        <int> 218, NA, NA, 218, 222, 218, NA, 218, 218, NA, 218, 218, NA,…
## $ X.10       <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 220, NA, NA, NA, 22…
## $ X.11       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ AG13       <int> 251, 251, 262, 262, 255, 255, 262, 268, 251, 255, 262, 264,…
## $ X.12       <int> 268, 266, 266, 266, 270, 264, 266, 272, 268, 262, 270, 268,…
## $ X.13       <int> NA, 270, 268, 268, 272, 268, 272, 282, NA, 268, 272, 272, 2…
## $ X.14       <int> NA, 272, NA, 272, 286, 272, NA, NA, NA, NA, NA, NA, 280, NA…
## $ AG23       <int> 356, 356, 356, 356, 356, 356, 356, 356, 358, 356, 356, 356,…
## $ X.15       <int> 360, NA, 358, 358, 358, 358, 358, 358, 358, NA, 358, 358, N…
## $ X.16       <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ X.17       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ AG25       <int> 93, 93, 93, 93, 93, 99, 99, 93, 93, 99, 93, 93, 91, 93, 93,…
## $ X.18       <int> 93, 105, 99, 99, NA, 101, 105, 99, 103, 101, NA, 101, 93, 1…
## $ X.19       <int> NA, NA, 103, NA, NA, 113, NA, 109, NA, 111, NA, NA, 103, 11…
## $ X.20       <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, 113, NA, NA, NA, NA, NA…
## $ AG27       <int> 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101,…
## $ X.21       <int> 101, 103, 103, NA, 103, NA, NA, NA, 101, NA, NA, NA, NA, 10…
## $ X.22       <int> NA, 105, 105, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 103, …
## $ X.23       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ AG30       <int> 99, 97, 97, 97, 99, 99, 97, 97, 99, 97, 97, 97, 97, 97, 97,…
## $ X.24       <int> 107, 99, 99, 105, NA, 101, 101, 99, 101, 99, 107, 99, 105, …
## $ X.25       <int> NA, 103, 103, NA, NA, 109, 111, 109, NA, NA, 109, NA, NA, 1…
## $ X.26       <int> NA, NA, NA, NA, NA, NA, 113, NA, NA, NA, 113, NA, NA, 109, …
## $ AG35       <int> 196, 184, 181, 179, 179, 179, 183, 186, 181, 183, 188, 186,…
## $ X.27       <int> 202, 200, 186, 184, 181, 192, 186, 192, 184, 186, 192, 196,…
## $ X.28       <int> NA, NA, 190, 186, 192, 198, 198, 194, NA, 192, 194, NA, 194…
## $ X.29       <int> NA, NA, 200, 194, NA, NA, 202, NA, NA, 194, NA, NA, NA, 200…
## $ A2         <int> 142, 142, 142, 140, 142, 142, 142, 142, 142, 142, 142, 142,…
## $ X.30       <int> NA, NA, NA, 142, NA, NA, NA, NA, NA, NA, NA, 146, NA, NA, N…
## $ X.31       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ X.32       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ A22        <int> 161, 159, 159, 159, 161, 159, 159, 159, 161, 159, 161, 155,…
## $ X.33       <int> NA, 163, 163, 163, 177, 161, 161, 163, NA, 165, 163, 159, 1…
## $ X.34       <int> NA, 167, NA, 175, NA, NA, 165, 167, NA, 167, 169, 173, 163,…
## $ X.35       <int> NA, 169, NA, NA, NA, NA, NA, NA, NA, NA, NA, 175, NA, NA, 1…
## $ A35        <int> 228, 224, 224, 224, 224, 218, 218, 224, 224, 226, 224, 224,…
## $ X.36       <int> NA, 228, 226, 226, 226, 226, 222, 230, 228, 238, 226, 228, …
## $ X.37       <int> NA, 230, 230, 228, 232, 232, 226, 232, NA, NA, 228, NA, 232…
## $ X.38       <int> NA, NA, NA, 230, 238, NA, 230, 238, NA, NA, NA, NA, NA, 232…
## $ A37        <int> 249, 239, 239, 239, 239, 239, 239, 239, 249, 239, 239, 243,…
## $ X.39       <int> 255, 247, 249, 243, 243, 243, 247, 251, 255, 243, 243, 245,…
## $ X.40       <int> NA, NA, 251, 251, 245, 251, 249, NA, NA, 251, 265, 247, 249…
## $ X.41       <int> NA, NA, 265, NA, 249, NA, NA, NA, NA, NA, NA, 251, NA, NA, …
## $ A38        <int> 119, 110, 119, 113, 113, 119, 108, 113, 108, 113, 119, 119,…
## $ X.42       <int> 136, 119, 122, 122, 119, 122, 119, 119, 110, 119, 122, 125,…
## $ X.43       <int> NA, 122, 148, 125, 122, 125, NA, 122, NA, 122, 125, 148, 12…
## $ X.44       <int> NA, NA, NA, 148, NA, 133, NA, NA, NA, NA, NA, NA, NA, 136, …
## $ AG20       <int> 304, 302, 302, 304, 302, 302, 302, 302, 304, 302, 304, 302,…
## $ X.45       <int> NA, 304, 306, 306, 308, 304, 306, 306, 306, 306, 306, 318, …
## $ X.46       <int> NA, 306, 308, 308, NA, 308, NA, 308, NA, 308, 308, NA, 318,…
## $ X.47       <int> NA, NA, 318, NA, NA, NA, NA, NA, NA, 318, NA, NA, NA, NA, N…
```



```r
if(dat %>%
  summarise(tiene_datos = any(!is.na(pull(., last_col())))) %>% pull()) {
  dat <-  dat
  
  } else {dat <- dat %>% select(-last_col())}
```

Finally, an object containing only the loci data (in columns) will be created.  

```r
# Ask the user how many columns to remove from the beginning
#n_col_remove <- as.numeric(readline(prompt = "How many columns do you want to #remove from the start? "))

n_col_remove = 3

# Optional: validate the input
if (is.na(n_col_remove) || n_col_remove < 1 || n_col_remove > ncol(dat)) {
  stop("Invalid number of columns.")
}

# Keep those col names for later use (in the proper format) 
cols_recover = str_to_lower(names(dat)[1:n_col_remove])
cols_recover = stringi::stri_trans_general(cols_recover, "Latin-ASCII")

# Remove the specified number of columns from the start of the data frame
dat_pr <- dat[, -(1:n_col_remove)]

# Show the result
dat_pr %>% tibble()
```

```
## # A tibble: 276 × 64
##     AG01     X   X.1   X.2  AG05   X.3   X.4   X.5  AG09   X.6   X.7   X.8  AG10
##    <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
##  1   134   146    NA    NA   148   159    NA    NA   250   254    NA    NA   216
##  2   132    NA    NA    NA   141   155   157    NA   246   248   250    NA   218
##  3   132    NA    NA    NA   141   146   157    NA   244   246   248    NA   218
##  4   132    NA    NA    NA   146   151   155    NA   244   248   250   256   216
##  5   132    NA    NA    NA   141   146    NA    NA   246    NA    NA    NA   218
##  6   132   134    NA    NA   141   155    NA    NA   244   246   262    NA   216
##  7   128   132   166    NA   141   148   155    NA   244   246   250    NA   218
##  8   132    NA    NA    NA   141   148   155    NA   246   248    NA    NA   216
##  9   134   135    NA    NA   151   159    NA    NA   250   250    NA    NA   216
## 10   126   132    NA    NA   148   155   171    NA   246    NA    NA    NA   218
## # ℹ 266 more rows
## # ℹ 51 more variables: X.9 <int>, X.10 <int>, X.11 <lgl>, AG13 <int>,
## #   X.12 <int>, X.13 <int>, X.14 <int>, AG23 <int>, X.15 <int>, X.16 <int>,
## #   X.17 <lgl>, AG25 <int>, X.18 <int>, X.19 <int>, X.20 <int>, AG27 <int>,
## #   X.21 <int>, X.22 <int>, X.23 <lgl>, AG30 <int>, X.24 <int>, X.25 <int>,
## #   X.26 <int>, AG35 <int>, X.27 <int>, X.28 <int>, X.29 <int>, A2 <int>,
## #   X.30 <int>, X.31 <lgl>, X.32 <lgl>, A22 <int>, X.33 <int>, X.34 <int>, …
```

## Build the matrix for polyploids
The code assumes that some samples are 4x. It is important to check that no marker name contains "X"

```r
# names of the markers
markers = dat_pr %>% names()

# IMPORTANT : make sure that no name contains "X"
markers = markers[!str_detect(markers, "X")]

# Optional: check the markers name input
if (sum(str_detect(markers, "X")) > 0) {
  stop("Invalid marker name : it contains 'X'")
}

start = which(str_detect(dat_pr %>% names(), "A"))
end = c(start - 1, length(dat_pr))

# internal coding
tmp <- sapply(1:length(markers), function(i) str_c("A", i), USE.NAMES = FALSE)
```



```r
nAG <- seq_along(tmp)

for(m in nAG) {
  marker = dat_pr[,start[m]:end[m+1]]
  
  one = c()
  for(i in 1:nrow(marker)) {
    # combine into vector
    one = c(one, c(t(marker[i,])))
    
  }
  if(m == 1) {
    mydat = data.frame(one)
  } else {
    mydat = bind_cols(mydat, one)
  }
  
}
```

```
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## New names:
## • `` -> `...2`
```

```r
# Recover original marker names 
names(mydat) = markers

mydat %>% as_tibble()
```

```
## # A tibble: 1,104 × 16
##     AG01  AG05  AG09  AG10  AG13  AG23  AG25  AG27  AG30  AG35    A2   A22   A35
##    <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
##  1   134   148   250   216   251   356    93   101    99   196   142   161   228
##  2   146   159   254   218   268   360    93   101   107   202    NA    NA    NA
##  3    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
##  4    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
##  5   132   141   246   218   251   356    93   101    97   184   142   159   224
##  6    NA   155   248    NA   266    NA   105   103    99   200    NA   163   228
##  7    NA   157   250    NA   270    NA    NA   105   103    NA    NA   167   230
##  8    NA    NA    NA    NA   272    NA    NA    NA    NA    NA    NA   169    NA
##  9   132   141   244   218   262   356    93   101    97   181   142   159   224
## 10    NA   146   246    NA   266   358    99   103    99   186    NA   163   226
## # ℹ 1,094 more rows
## # ℹ 3 more variables: A37 <int>, A38 <int>, AG20 <int>
```

Combine dataframe with the columns other than loci information.   

```r
new_cols <- purrr::imap_dfc(cols_recover, function(colname, i) {
  # i es el índice de la columna en dat
  tibble(!!colname := rep(dat[[i]], each = 4))
})

# Combinamos con mydat
mydat <- bind_cols(new_cols, mydat)

mydat %>% tibble()
```

```
## # A tibble: 1,104 × 19
##    sample population ploidia  AG01  AG05  AG09  AG10  AG13  AG23  AG25  AG27
##    <chr>  <chr>      <chr>   <int> <int> <int> <int> <int> <int> <int> <int>
##  1 T-481  Madres     2x        134   148   250   216   251   356    93   101
##  2 T-481  Madres     2x        146   159   254   218   268   360    93   101
##  3 T-481  Madres     2x         NA    NA    NA    NA    NA    NA    NA    NA
##  4 T-481  Madres     2x         NA    NA    NA    NA    NA    NA    NA    NA
##  5 T-482  Madres     4x        132   141   246   218   251   356    93   101
##  6 T-482  Madres     4x         NA   155   248    NA   266    NA   105   103
##  7 T-482  Madres     4x         NA   157   250    NA   270    NA    NA   105
##  8 T-482  Madres     4x         NA    NA    NA    NA   272    NA    NA    NA
##  9 T-483  Madres     4x        132   141   244   218   262   356    93   101
## 10 T-483  Madres     4x         NA   146   246    NA   266   358    99   103
## # ℹ 1,094 more rows
## # ℹ 8 more variables: AG30 <int>, AG35 <int>, A2 <int>, A22 <int>, A35 <int>,
## #   A37 <int>, A38 <int>, AG20 <int>
```

Replace missing data with special value "-5" (STRUCTURE syntax)

```r
mydat <-  mydat %>% 
  group_by(sample) %>% 
  mutate(across(all_of(markers), 
                ~ if (all(is.na(.))) -5 else .,
                .names = "{.col}"))
```


Replace NA with special value "-9" (STRUCTURE syntax)

```r
mydat[is.na(mydat)] = -9
```

Export data 

```r
write.csv(mydat, file = "res/dataset.csv", row.names = F)
```

**Note**: This dataset should be accessible to the user for downloading.  
  

Finally, build the STRUCTURE input file. STRUCTURE requires an specific format that is built in the next code : 

```r
res <-  mydat %>% select(-3)
n_cols <- ncol(res)

# Define header manually (without modifying 'res')
col_names <- c("", "", rep("-9", ncol(res) - 2))

# Open connection for writing
#file_conn <- file("res/structure_input2025.txt", open = "wt")
file_conn <- file("res/resto_no.txt", open = "wt")

# Write header manually
writeLines(paste(col_names, collapse = "\t"), con = file_conn)

# Write data without header and without quotes
write.table(res, file = file_conn,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Close conexion
close(file_conn)
```


---
output: 
  html_document: 
    keep_md: yes
---
title: "binary matrix conversion"
output: html_notebook
editor_options: 
  chunk_output_type: inline
---

## Rationale
This code processes a dataset containing allele data per locus for multiple individuals and converts it into a binary matrix. Each column in the binary matrix represents the presence (1) or absence (0) of a specific allele at a given locus.



```r
library(polysat)
library(dplyr)
library(tidyr)
library(stringr)
```

Load data 

```r
dat <- read.csv(file = "res/dataset.csv")
dat %>% tibble()
```

```
## # A tibble: 1,104 × 19
##    sample population ploidia  AG01  AG05  AG09  AG10  AG13  AG23  AG25  AG27
##    <chr>  <chr>      <chr>   <int> <int> <int> <int> <int> <int> <int> <int>
##  1 T-481  Madres     2x        134   148   250   216   251   356    93   101
##  2 T-481  Madres     2x        146   159   254   218   268   360    93   101
##  3 T-481  Madres     2x         -9    -9    -9    -9    -9    -9    -9    -9
##  4 T-481  Madres     2x         -9    -9    -9    -9    -9    -9    -9    -9
##  5 T-482  Madres     4x        132   141   246   218   251   356    93   101
##  6 T-482  Madres     4x         -9   155   248    -9   266    -9   105   103
##  7 T-482  Madres     4x         -9   157   250    -9   270    -9    -9   105
##  8 T-482  Madres     4x         -9    -9    -9    -9   272    -9    -9    -9
##  9 T-483  Madres     4x        132   141   244   218   262   356    93   101
## 10 T-483  Madres     4x         -9   146   246    -9   266   358    99   103
## # ℹ 1,094 more rows
## # ℹ 8 more variables: AG30 <int>, AG35 <int>, A2 <int>, A22 <int>, A35 <int>,
## #   A37 <int>, A38 <int>, AG20 <int>
```

Pre-processing   
If missing data are coded as -9, convert these values to NA for proper handling in R. 

```r
if(any(dat==-9)) {
  dat[dat == -9] <- NA
}
```


### Step 1: Make long format
Here, we reshape the dataset from wide to long format. We assume the first `n_col_remove` columns are metadata and exclude them from reshaping.


```r
# Ask the user how many columns to remove from the beginning
#n_col_remove <- as.numeric(readline(prompt = "How many columns do you want to #remove from the start? "))

# Number of columns to exclude from reshaping (usually metadata columns)
n_col_remove = 3

long <- dat %>%
  pivot_longer(cols = -c(1:n_col_remove), names_to = "locus", 
               values_to = "allele") %>%
  filter(!is.na(allele))

long
```

```
## # A tibble: 9,534 × 5
##    sample population ploidia locus allele
##    <chr>  <chr>      <chr>   <chr>  <int>
##  1 T-481  Madres     2x      AG01     134
##  2 T-481  Madres     2x      AG05     148
##  3 T-481  Madres     2x      AG09     250
##  4 T-481  Madres     2x      AG10     216
##  5 T-481  Madres     2x      AG13     251
##  6 T-481  Madres     2x      AG23     356
##  7 T-481  Madres     2x      AG25      93
##  8 T-481  Madres     2x      AG27     101
##  9 T-481  Madres     2x      AG30      99
## 10 T-481  Madres     2x      AG35     196
## # ℹ 9,524 more rows
```

### Step 2: make column "locus_allele"
We create a new column, `marker`, which uniquely identifies each allele at each locus by concatenating locus and allele information. We also add a `value` column set to 1 to indicate presence.

```r
long <- long  %>% 
  mutate(marker = paste0(locus, "_", allele),
         value = 1L)

long
```

```
## # A tibble: 9,534 × 7
##    sample population ploidia locus allele marker   value
##    <chr>  <chr>      <chr>   <chr>  <int> <chr>    <int>
##  1 T-481  Madres     2x      AG01     134 AG01_134     1
##  2 T-481  Madres     2x      AG05     148 AG05_148     1
##  3 T-481  Madres     2x      AG09     250 AG09_250     1
##  4 T-481  Madres     2x      AG10     216 AG10_216     1
##  5 T-481  Madres     2x      AG13     251 AG13_251     1
##  6 T-481  Madres     2x      AG23     356 AG23_356     1
##  7 T-481  Madres     2x      AG25      93 AG25_93      1
##  8 T-481  Madres     2x      AG27     101 AG27_101     1
##  9 T-481  Madres     2x      AG30      99 AG30_99      1
## 10 T-481  Madres     2x      AG35     196 AG35_196     1
## # ℹ 9,524 more rows
```

### Step 3: remove homocygotes

Homozygous genotypes (identical alleles at the same locus) can cause problems in binary pivoting. Here is an example filtering step for a specific sample and locus.


```r
# Example inspection for potential homozygotes
long %>% 
  filter(sample %in% c("T-489"),
         locus %in% c("AG09"))
```

```
## # A tibble: 2 × 7
##   sample population ploidia locus allele marker   value
##   <chr>  <chr>      <chr>   <chr>  <int> <chr>    <int>
## 1 T-489  Madres     2x      AG09     250 AG09_250     1
## 2 T-489  Madres     2x      AG09     250 AG09_250     1
```

```r
# Remove duplicate rows to avoid issues during pivoting
long <- long %>% distinct()
```



### Step 4: Pivotear a formato binario
Finally, we convert the long-format data into a wide binary matrix. Each column corresponds to a unique marker, and cells contain 1 (presence) or 0 (absence). The matrix keeps sample and population as identifier columns.

```r
binary <- pivot_wider(
  long, id_cols = c(sample, population),
  names_from = marker,
  values_from = value,
  values_fill = 0) %>% 
  select(sample, population, sort(names(.))) 
```


### Result
Display the resulting binary matrix:

```r
print(binary)
```

```
## # A tibble: 276 × 170
##    sample population A2_140 A2_142 A2_146 A22_155 A22_157 A22_159 A22_161
##    <chr>  <chr>       <int>  <int>  <int>   <int>   <int>   <int>   <int>
##  1 T-481  Madres          0      1      0       0       0       0       1
##  2 T-482  Madres          0      1      0       0       0       1       0
##  3 T-483  Madres          0      1      0       0       0       1       0
##  4 T-484  Madres          1      1      0       0       0       1       0
##  5 T-485  Madres          0      1      0       0       0       0       1
##  6 T-486  Madres          0      1      0       0       0       1       1
##  7 T-487  Madres          0      1      0       0       0       1       1
##  8 T-488  Madres          0      1      0       0       0       1       0
##  9 T-489  Madres          0      1      0       0       0       0       1
## 10 T-490  Madres          0      1      0       0       0       1       0
## # ℹ 266 more rows
## # ℹ 161 more variables: A22_163 <int>, A22_165 <int>, A22_167 <int>,
## #   A22_169 <int>, A22_171 <int>, A22_173 <int>, A22_175 <int>, A22_177 <int>,
## #   A22_179 <int>, A35_210 <int>, A35_214 <int>, A35_218 <int>, A35_220 <int>,
## #   A35_222 <int>, A35_224 <int>, A35_226 <int>, A35_228 <int>, A35_230 <int>,
## #   A35_232 <int>, A35_234 <int>, A35_236 <int>, A35_238 <int>, A35_240 <int>,
## #   A35_242 <int>, A37_239 <int>, A37_243 <int>, A37_245 <int>, …
```



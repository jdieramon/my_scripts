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
From the binary matrix, we can further apply Principal Component Analysis (PCA) to visualize the similarity or dissimilarity among the samples.     
  


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
if(any(dat == -9)) {
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
## # A tibble: 9,554 × 5
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
## # ℹ 9,544 more rows
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
## # A tibble: 9,554 × 7
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
## # ℹ 9,544 more rows
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
## # A tibble: 276 × 174
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
## # ℹ 165 more variables: A22_163 <int>, A22_165 <int>, A22_167 <int>,
## #   A22_169 <int>, A22_171 <int>, A22_173 <int>, A22_175 <int>, A22_177 <int>,
## #   A22_179 <int>, A35_210 <int>, A35_214 <int>, A35_218 <int>, A35_220 <int>,
## #   A35_222 <int>, A35_224 <int>, A35_226 <int>, A35_228 <int>, A35_230 <int>,
## #   A35_232 <int>, A35_234 <int>, A35_236 <int>, A35_238 <int>, A35_240 <int>,
## #   A35_242 <int>, A37_239 <int>, A37_243 <int>, A37_245 <int>, …
```



### PCA

```r
# Ask the user how many columns to remove from the beginning
#n_col_remove <- as.numeric(readline(prompt = "How many columns do you want to #remove from the start? "))

n_col_remove = 2

# Convert the dataframe to a matrix
bin <- as.matrix(binary[,-c(1:n_col_remove)])

# Column names should be of form "locus.allele"
colnames(bin) <- str_trim(str_replace(colnames(bin), "_", "."))

# check potentital typos and fix, if any
colnames(bin)
```

```
##   [1] "A2.140"   "A2.142"   "A2.146"   "A22.155"  "A22.157"  "A22.159" 
##   [7] "A22.161"  "A22.163"  "A22.165"  "A22.167"  "A22.169"  "A22.171" 
##  [13] "A22.173"  "A22.175"  "A22.177"  "A22.179"  "A35.210"  "A35.214" 
##  [19] "A35.218"  "A35.220"  "A35.222"  "A35.224"  "A35.226"  "A35.228" 
##  [25] "A35.230"  "A35.232"  "A35.234"  "A35.236"  "A35.238"  "A35.240" 
##  [31] "A35.242"  "A37.239"  "A37.243"  "A37.245"  "A37.247"  "A37.249" 
##  [37] "A37.251"  "A37.253"  "A37.255"  "A37.257"  "A37.259"  "A37.261" 
##  [43] "A37.265"  "A38.108"  "A38.110"  "A38.113"  "A38.116"  "A38.119" 
##  [49] "A38.122"  "A38.125"  "A38.133"  "A38.136"  "A38.145"  "A38.148" 
##  [55] "A38.150"  "AG01.-5"  "AG01.126" "AG01.128" "AG01.130" "AG01.132"
##  [61] "AG01.134" "AG01.135" "AG01.139" "AG01.142" "AG01.144" "AG01.146"
##  [67] "AG01.162" "AG01.164" "AG01.166" "AG01.168" "AG05.141" "AG05.144"
##  [73] "AG05.146" "AG05.148" "AG05.151" "AG05.153" "AG05.155" "AG05.157"
##  [79] "AG05.159" "AG05.161" "AG05.163" "AG05.169" "AG05.171" "AG09.242"
##  [85] "AG09.244" "AG09.246" "AG09.248" "AG09.250" "AG09.252" "AG09.254"
##  [91] "AG09.256" "AG09.258" "AG09.262" "AG10.216" "AG10.218" "AG10.220"
##  [97] "AG10.222" "AG10.224" "AG10.232" "AG13.-5"  "AG13.251" "AG13.253"
## [103] "AG13.255" "AG13.257" "AG13.262" "AG13.264" "AG13.266" "AG13.268"
## [109] "AG13.270" "AG13.272" "AG13.274" "AG13.280" "AG13.282" "AG13.286"
## [115] "AG20.302" "AG20.304" "AG20.306" "AG20.308" "AG20.312" "AG20.314"
## [121] "AG20.316" "AG20.318" "AG20.320" "AG23.-5"  "AG23.356" "AG23.358"
## [127] "AG23.360" "AG25.101" "AG25.103" "AG25.105" "AG25.107" "AG25.109"
## [133] "AG25.111" "AG25.113" "AG25.91"  "AG25.93"  "AG25.95"  "AG25.97" 
## [139] "AG25.99"  "AG27.101" "AG27.103" "AG27.105" "AG27.108" "AG27.116"
## [145] "AG27.99"  "AG30.101" "AG30.103" "AG30.105" "AG30.107" "AG30.109"
## [151] "AG30.111" "AG30.113" "AG30.97"  "AG30.99"  "AG35.-5"  "AG35.173"
## [157] "AG35.175" "AG35.177" "AG35.179" "AG35.181" "AG35.183" "AG35.184"
## [163] "AG35.186" "AG35.188" "AG35.190" "AG35.192" "AG35.194" "AG35.196"
## [169] "AG35.198" "AG35.200" "AG35.202" "AG35.204"
```

```r
# Calculate the distance
d <- dist(bin)

# Apply classical scaling (Classical Multidimensional Scaling, MDS)
coords <- cmdscale(d, k = 2)  # k = número de dimensiones

# Plot 'coords' : a matrix with 2D coordinates
plot(coords, main = "Principal Coordinates Analysis")  # PCoA
```

![](binary_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Another way : 

```r
coords <- cmdscale(d, eig = T)
variance = coords$eig / sum(coords$eig)
variance1 <- 100 * signif(variance[1], 2)
variance2 <- 100 * signif(variance[2], 2)
```
Represented by ploidy : 

```r
#Extract the ploidy level
dat <- dat %>% select(1:3) %>% distinct()

# PCoA
coords <- cmdscale(d, k = 2)  # k = número de dimensiones
mycol = as.numeric(as.factor(dat$ploidia))
plot(coords, main = "Principal Coordinates Analysis", # PCoA
     pch = 19, 
     col = scales::alpha(mycol, 0.6), 
     xlab = paste("PCA1", paste0("(", variance1, "%)")), 
     ylab = paste("PCA2", paste0("(", variance2, "%)")))  
```

![](binary_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Represented by population : 

```r
mycol = as.numeric(as.factor(dat$population))
plot(coords, main = "Principal Coordinates Analysis", # PCoA
     pch = 19, 
     col = scales::alpha(mycol, 0.6), 
     xlab = paste("PCA1", paste0("(", variance1, "%)")), 
     ylab = paste("PCA2", paste0("(", variance2, "%)")))  
```

![](binary_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

PCA

```r
pca <- prcomp(bin, scale. = F)
biplot(pca)

pca <- prcomp(bin, scale. = T)
biplot(pca)
```


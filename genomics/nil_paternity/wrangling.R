# Dependencies
library(stringr)
library(tibble)
library(dplyr)

# Download datasets with SNP
#https://cegresources.icrisat.org/cicerseq/?page_id=3444

# Download files with the SNP data for the cultivated chickpea accessions 
url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca1.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca1.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca2.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca2.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca3.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca3.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca4.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca4.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca5.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca5.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca6.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca6.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca7.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca7.hmp.txt.gz")

url = "https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca8.hmp.txt.gz"
download.file(url = url, destfile = "data/Cultivated_Ca8.hmp.txt.gz")

ls ~/Downloads/*.gz
mv ~/Downloads/*.gz data


# Uncompress data
junk <- list.files(path = "data/", pattern = "*.hmp.txt.gz", full.names = TRUE)
sapply(junk, function(x) R.utils::gunzip(filename = x, remove = TRUE))
rm(junk)


# Load the passport data 
passport <- read.csv("data/passport_results.csv")


# Extract passport data for a given genotype -----------------------------------

# Use the pangenome_app.R for passport data 
# https://genehummus.shinyapps.io/pangenome/

# and extract the column data containing the SNPs for that paricular genotype. 
# that colmn number will be needed in the next section 

  
  
# SNP pangenome -------------------------------------------------------
# Build the data frame for selected genotypes (Linux)
# use the column numbers obtained in the previous section (pangenome_app.R)

# cut -f1,-5,3134,3181,103 data/Cultivated_Ca1.hmp.txt > data/ca1_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca2.hmp.txt > data/ca2_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca3.hmp.txt > data/ca3_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca4.hmp.txt > data/ca4_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca5.hmp.txt > data/ca5_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca6.hmp.txt > data/ca6_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca7.hmp.txt > data/ca7_genotypes.txt #~2Gb
# cut -f1,-5,3134,3181,103 data/Cultivated_Ca8.hmp.txt > data/ca8_genotypes.txt #~1Gb

#mkdir cultivated
#mv data/*.hmp.txt cultivated/
#mv cultivated/ ~/Desktop/


# NILs ------------------------------------------------------------------------
# Uncompress *. vcf file
vcf_file = "data/TORRESANA_01.genotyped.snpEff.p.GATKfilters.SAL.SAL10_1.diffGT.vcf.gz"
R.utils::gunzip(filename = vcf_file, remove = TRUE)
# Rename
file.rename(from = "data/TORRESANA_01.genotyped.snpEff.p.GATKfilters.SAL.SAL10_1.diffGT.vcf", 
            to = "data/nil.vcf")

# Remove header from vcf file 
egrep -v "^#" data/nil.cvf

# Get the line for the columns header (remove the header #)
system("grep -n CHROM data/nil.vcf")

# Get first data (from that number of line)
system("sed -n 7196,7200p data/nil.vcf | less -S")

# See SNPs for chromosome 8 (f. ex)
grep -E "CHROM|NC_021167.1" data/nil.vcf | head | less -S
grep -E "CHROM|NC_021167.1" data/nil.vcf > data/nils8.txt

# Create files for each chromosome 
grep -E "CHROM|NC_021160.1" data/nil.vcf > data/nils1.txt
grep -E "CHROM|NC_021161.1" data/nil.vcf > data/nils2.txt
etc ...
#edit with nano : remove ##contig=<ID=NC_021160.1,length=48359943>



# Read NIL files 
nils_ds = read.delim2("data/nils1.txt")   ## -> select chr1 file (NIL)
# nils_ds = read.delim2("data/nils2.txt") ## -> select chr2 file (NIL)
# nils_ds = read.delim2("data/nils3.txt") ## -> select chr3 file (NIL)
dim(nils_ds)


# Pre-processing
##filter by PASS 
## clean columns NILS
nils_ds <- nils_ds %>% 
  as_tibble() %>% 
  filter(FILTER == "PASS") %>%
  select(-c(QUAL, FILTER, INFO, FORMAT)) %>% 
  mutate(AN0205 = str_sub(AN0205, 1, 3), 
         AN0207 = str_sub(AN0207, 1, 3))

dim(nils_ds)

# enconde 0/1 by nucleotide : ver posibles combinaciones a codificar 
table(nils_ds$AN0205)
table(nils_ds$AN0207)

nils_ds <- nils_ds %>% 
  mutate(AN0205 = case_match(AN0205, 
                             "0/0" ~ REF, 
                             "0/1" ~ paste0(REF, "/", ALT), 
                             "0/2" ~ paste0(REF, "/", str_split(ALT, ",")[[1]][2]),
                             "1/1" ~ ALT, 
                             "1/2" ~ paste0(str_split(ALT, ",")[[1]][1], "/", str_split(ALT, ",")[[1]][2]),
                             "2/2" ~ str_replace(ALT, pattern = ",", replacement = "/")), 
         AN0207 = case_match(AN0207, 
                             "0/0" ~ REF, 
                             "0/1" ~ paste0(REF, "/", ALT), 
                             "1/1" ~ ALT,
                             "0/2" ~ paste0(REF, "/", str_split(ALT, ",")[[1]][2]), 
                             "2/2" ~ str_replace(ALT, pattern = ",", replacement = "/"), 
                             "1/2" ~ paste0(str_split(ALT, ",")[[1]][1], "/", str_split(ALT, ",")[[1]][2])) ) %>% 
  janitor::clean_names() %>% 
  rename(chr = x_chrom, 
         nil_tardia = an0205, 
         nil_precoz = an0207)

dim(nils_ds)


# Pangenome ---------------------------------------

ca = read.delim("data/ca1_genotypes.txt")   ## -> select chr1 file (pangenome)
# ca = read.delim("data/ca2_genotypes.txt") ## -> select chr2 file (pangenome)
# ca = read.delim("data/ca3_genotypes.txt") ## -> select chr3 file (pangenome)

# Pre-processing 
ca = ca %>% 
  as_tibble() %>% 
  mutate(CPAM.88 = case_match(CPAM.88, 
                              "R" ~ "A/G", 
                              "S" ~ "G/C", 
                              "Y" ~ "T/C", 
                              "K" ~"G/T", 
                              "M" ~"A/C", 
                              "W" ~"A/T", 
                              .default = CPAM.88), 
         JG.62 = case_match(JG.62, 
                            "R" ~ "A/G", 
                            "S" ~ "G/C", 
                            "Y" ~ "T/C", 
                            "K" ~"G/T", 
                            "M" ~"A/C", 
                            "W" ~"A/T", 
                            .default = JG.62), 
         WR315 = case_match(WR315, 
                            "R" ~ "A/G", 
                            "S" ~ "G/C", 
                            "Y" ~ "T/C", 
                            "K" ~"G/T", 
                            "M" ~"A/C", 
                            "W" ~"A/T", 
                            .default = WR315))


# posiciones con el mismo SNP de cada parental no aportan información
# simplifico dataset eliminando esas coordenadas 
non_info <- ca %>% 
  filter(CPAM.88 == JG.62) %>% 
  filter(CPAM.88 == WR315) %>% 
  pull(pos)

# check all genotypes with same allele
ca %>% filter(pos %in% sample(non_info, 5))


ca_filtered <- ca %>% 
  filter(!pos %in% non_info)

ca_filtered

# Combine caCHROMOSOME_filtered and nils_snps
ds <- nils_ds %>% 
  inner_join(ca_filtered, by = "pos") %>% 
  select(-c(id, rs., alt, alleles))

# Number of SNPs shared between NILs and parents (pangenome)
nrow(ds)

# heterozygous SNPS
idx = which(str_detect(ds$nil_tardia, "/"))
pos_heteroz = ds %>% slice(idx) %>% pull(pos)

idx = which(str_detect(ds$nil_precoz, "/"))
pos_heteroz = c(pos_heteroz, ds %>% slice(idx) %>% pull(pos))


## homozygous SNP ------
ds_homoz <- ds %>% 
  filter(! pos %in% pos_heteroz) %>% 
  mutate(ilc72_tardia = ifelse(nil_tardia == CPAM.88, 1, 0), 
         ilc72_precoz = ifelse(nil_precoz == CPAM.88, 1, 0), 
         jg62_tardia  = ifelse(nil_tardia == JG.62, 1, 0), 
         jg62_precoz  = ifelse(nil_precoz == JG.62, 1, 0),
         wr315_tardia = ifelse(nil_tardia == WR315, 1, 0), 
         wr315_precoz = ifelse(nil_precoz == WR315, 1, 0)) %>% 
  # fix by "N" : all bases available
  mutate(ilc72_tardia = case_when(CPAM.88 == "N" ~ 1, 
                                  TRUE ~ ilc72_tardia), 
         ilc72_precoz = case_when(CPAM.88 == "N" ~ 1, 
                                  TRUE ~ ilc72_precoz), 
         jg62_tardia = case_when(JG.62 == "N" ~ 1, 
                                 TRUE ~ jg62_tardia), 
         jg62_precoz = case_when(JG.62 == "N" ~ 1, 
                                 TRUE ~ jg62_precoz), 
         wr315_tardia = case_when(WR315 == "N" ~ 1, 
                                  TRUE ~ wr315_tardia), 
         wr315_precoz = case_when(WR315 == "N" ~ 1, 
                                  TRUE ~ wr315_precoz))


dim(ds_homoz)

## heteroztgous SNP ------
ds_res_heteroz <- ds %>% 
  filter(pos %in% pos_heteroz) %>% 
  
  separate(nil_tardia, into = c("nil_t1", "nil_t2"), sep = "/") %>% 
  separate(nil_precoz, into = c("nil_p1", "nil_p2"), sep = "/") %>% 
  
  #crea 2 alelos : por nil y parental
  mutate(ilc72_t1 = str_detect(CPAM.88, nil_t1), 
         ilc72_t2 = str_detect(CPAM.88, nil_t2), 
         ilc72_p1 = str_detect(CPAM.88, nil_p1), 
         ilc72_p2 = str_detect(CPAM.88, nil_p2), 
         
         jg62_t1 = str_detect(JG.62, nil_t1), 
         jg62_t2 = str_detect(JG.62, nil_t2), 
         jg62_p1 = str_detect(JG.62, nil_p1), 
         jg62_p2 = str_detect(JG.62, nil_p2),
         
         wr315_t1 = str_detect(WR315, nil_t1), 
         wr315_t2 = str_detect(WR315, nil_t2), 
         wr315_p1 = str_detect(WR315, nil_p1), 
         wr315_p2 = str_detect(WR315, nil_p2) ) %>% 
  
  
  
  
  # fix by "N" : all bases available
  mutate(ilc72_t1 = case_when(CPAM.88 == "N" ~ TRUE,
                              TRUE ~ ilc72_t1), 
         ilc72_t2 = ifelse(!is.na(nil_t2), 
                           case_when(CPAM.88 == "N" ~ TRUE,
                                     TRUE ~ ilc72_t2), NA),
         ilc72_p1 = case_when(CPAM.88 == "N" ~ TRUE,
                              TRUE ~ ilc72_p1), 
         ilc72_p2 = ifelse(!is.na(nil_p2), 
                           case_when(CPAM.88 == "N" ~ TRUE,
                                     TRUE ~ ilc72_p2), NA), 
         
         
         jg62_t1 = case_when(JG.62 == "N" ~ TRUE,
                             TRUE ~ jg62_t1), 
         jg62_t2 = ifelse(!is.na(nil_t2), 
                          case_when(JG.62 == "N" ~ TRUE,
                                    TRUE ~ jg62_t2), NA),
         jg62_p1 = case_when(JG.62 == "N" ~ TRUE,
                             TRUE ~ jg62_p1), 
         jg62_p2 = ifelse(!is.na(nil_p2), 
                          case_when(JG.62 == "N" ~ TRUE,
                                    TRUE ~ jg62_p2), NA), 
         
         
         
         wr315_t1 = case_when(WR315 == "N" ~ TRUE,
                              TRUE ~ wr315_t1), 
         wr315_t2 = ifelse(!is.na(nil_t2), 
                           case_when(WR315 == "N" ~ TRUE,
                                     TRUE ~ wr315_t2), NA),
         wr315_p1 = case_when(WR315 == "N" ~ TRUE,
                              TRUE ~ wr315_p1), 
         wr315_p2 = ifelse(!is.na(nil_p2), 
                           case_when(WR315 == "N" ~ TRUE,
                                     TRUE ~ wr315_p2), NA), 
         
         
         
  ) %>%   
  
  
  # unifica alelos en parentales 
  rowwise() %>% 
  
  mutate(ilc72_tardia = sum(ilc72_t1, ilc72_t2, na.rm = T),
         ilc72_precoz = sum(ilc72_p1, ilc72_p2, na.rm = T), 
         jg62_tardia = sum(jg62_t1, jg62_t2, na.rm = T),
         jg62_precoz = sum(jg62_p1, jg62_p2, na.rm = T), 
         wr315_tardia = sum(wr315_t1, wr315_t2, na.rm = T),
         wr315_precoz = sum(wr315_p1, wr315_p2, na.rm = T)
         
  ) %>% 
  
  mutate(ilc72_tardia = ifelse(ilc72_tardia > 0, 1, ilc72_tardia), 
         ilc72_precoz = ifelse(ilc72_precoz > 0, 1, ilc72_precoz), 
         jg62_tardia = ifelse(jg62_tardia > 0, 1, jg62_tardia), 
         jg62_precoz = ifelse(jg62_precoz > 0, 1, jg62_precoz), 
         wr315_tardia = ifelse(wr315_tardia > 0, 1, wr315_tardia), 
         wr315_precoz = ifelse(wr315_precoz > 0, 1, wr315_precoz))

dim(ds_res_heteroz)

# combine data from homozygous and heterozygous SNPs
ds_res <- ds_homoz %>% 
  select(chr, pos, ref, strand, ilc72_tardia:wr315_precoz ) %>% 
  bind_rows(ds_res_heteroz %>% 
              select(chr, pos, ref, strand, ilc72_tardia:wr315_precoz)) %>% 
  arrange(chr, pos)

dim(ds_res)


# Resultado
# ver documento final : (Drive) SNP_nils_pangenome.gdoc

cat("NIL precoz comparte", 
    mean(ds_res$jg62_precoz) *100, "% de los SNPS con JG62")
cat("NIL precoz comparte", 
    mean(ds_res$ilc72_precoz)*100, "% de los SNPS con ILC72")
cat("NIL precoz comparte", 
    mean(ds_res$wr315_precoz)*100, "% de los SNPS con WR315")


cat("NIL tardía comparte", 
    mean(ds_res$jg62_tardia) *100, "% de los SNPS con JG62")
cat("NIL tardía comparte", 
    mean(ds_res$ilc72_tardia)*100, "% de los SNPS con ILC72")
cat("NIL tardía comparte", 
    mean(ds_res$wr315_tardia)*100, "% de los SNPS con WR315")



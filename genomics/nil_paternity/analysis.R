# general analysis combining all chromosomes 

# NILs -------------------------------------------

# Build dataset 
nils_ds <- nils_ds <- read.delim2("data/nils1.txt") %>% 
  bind_rows(read.delim("data/nils2.txt")) %>% 
  bind_rows(read.delim("data/nils3.txt")) %>% 
  bind_rows(read.delim("data/nils4.txt")) %>% 
  bind_rows(read.delim("data/nils5.txt")) %>% 
  bind_rows(read.delim("data/nils6.txt")) %>% 
  bind_rows(read.delim("data/nils7.txt")) %>% 
  bind_rows(read.delim("data/nils8.txt")) %>% 
  as_tibble()

dim(nils_ds)



# Pre-processing
nils_ds <- nils_ds %>% 
  filter(FILTER == "PASS") %>%
  select(-c(QUAL, FILTER, INFO, FORMAT)) %>% 
  mutate(AN0205 = str_sub(AN0205, 1, 3), 
         AN0207 = str_sub(AN0207, 1, 3)) %>% 
  # unify names : GenBank [.vcf file] -> CicerSeq
  mutate(X.CHROM = case_when(X.CHROM == "NC_021160.1" ~ "CA1", 
                             X.CHROM == "NC_021161.1" ~ "CA2", 
                             X.CHROM == "NC_021162.1" ~ "CA3",
                             X.CHROM == "NC_021163.1" ~ "CA4",
                             X.CHROM == "NC_021164.1" ~ "CA5",
                             X.CHROM == "NC_021165.1" ~ "CA6",
                             X.CHROM == "NC_021166.1" ~ "CA7",
                             X.CHROM == "NC_021167.1" ~ "CA8")) %>% 
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
ca = read.delim("data/ca1_genotypes.txt") %>% 
  bind_rows(read.delim("data/ca2_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca3_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca4_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca5_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca6_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca7_genotypes.txt")) %>% 
  bind_rows(read.delim("data/ca8_genotypes.txt"))




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
                            .default = WR315)) %>% 
  rename(chr = chrom)

dim(ca)

non_info <- ca %>% 
  filter(CPAM.88 == JG.62) %>% 
  filter(CPAM.88 == WR315) %>% 
  pull(pos)

# check all genotypes with same allele
ca %>% filter(pos %in% sample(non_info, 5))


ca_filtered <- ca %>% 
  filter(!pos %in% non_info)

dim(ca_filtered)


#COMBINED datasets
# Combine ca8_filtered and nils_snps
ds <- nils_ds %>% 
  inner_join(ca_filtered, by = c("chr", "pos")) %>% 
  select(-c(id, rs., alt, alleles))

# Number of SNPs shared between NILs and parents (pangenome)
nrow(ds)


idx = which(str_detect(ds$nil_tardia, "/"))
pos_heteroz = ds %>% slice(idx) %>% pull(pos)

idx = which(str_detect(ds$nil_precoz, "/"))
pos_heteroz = c(pos_heteroz, ds %>% slice(idx) %>% pull(pos))


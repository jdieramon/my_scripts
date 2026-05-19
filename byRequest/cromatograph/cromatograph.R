# Dependencies 
library(rvest)
library(stringr)
library(dplyr)
library(purrr)

# Read dataset
volatil <- read.csv("~/Desktop/volatil.csv", sep=";")
volatil <- read.csv("~/Desktop/nuevos_cas.csv", sep=";")
volatil[1:3, ]

# Extract CAS numbers as vector 
cas <- volatil$Hit.1..CAS
cas <- unique(cas)


# Define a function 

# Encuentra la sección correcta (aunque cambie ligeramente el texto)
# Evita errores si falta la tabla
# Devuelve NA si no hay datos
# Funciona con distintas variantes HTML del sitio



get_kovats_RI <- function(cas, phase = "DB-Wax") {
  
  url <- paste0(
    "https://webbook.nist.gov/cgi/cbook.cgi?ID=",
    cas,
    "&Units=SI&cGC=on"
  )
  
  # Intentar leer la página
  page <- tryCatch(read_html(url), error = function(e) return(NULL))
  if (is.null(page)) return(NA)
  
  # Extraer nodos (títulos + tablas)
  nodes <- page %>% html_elements(xpath = "//h3 | //table")
  texts <- nodes %>% html_text(trim = TRUE)
  
  # Buscar sección correcta
  idx <- which(str_detect(
    texts,
    regex("Kovats.*polar column.*temperature ramp", ignore_case = TRUE)
  ))
  
  if (length(idx) == 0) return(NA)
  
  # Buscar la siguiente tabla válida
  for (i in idx) {
    if ((i + 1) <= length(nodes)) {
      
      tbl <- tryCatch(
        nodes[[i + 1]] %>% html_table(fill = TRUE),
        error = function(e) NULL
      )
      
      if (!is.null(tbl)) {
        
        # Normalizar nombres de columnas
        colnames(tbl) <- str_trim(colnames(tbl))
        
        if (!("Active phase" %in% colnames(tbl)) || !("I" %in% colnames(tbl))) {
          next
        }
        
        # Filtrar fase
        row <- tbl %>%
          filter(str_detect(`Active phase`, fixed(phase, ignore_case = TRUE)))
        
        if (nrow(row) > 0) {
          return(as.numeric(row$I[1]))
        }
      }
    }
  }
  
  return(NA)
}



# Usage 
get_kovats_RI("66-25-1")
get_kovats_RI("112-40-3")

sapply(c("66-25-1","112-40-3"), function(i) get_kovats_RI(i))


# Uso en lote 
#cas_list <- c("66-25-1", "64-17-5", "67-56-1")
cas_list <- cas
results <- data.frame(
  CAS = cas_list,
  RI = map_dbl(cas_list, get_kovats_RI)
)

results

#vals <- sapply(cas, function(i) get_kovats_RI(i))


# Save the results into a .csv file
write.csv(results, file = "~/Desktop/results.csv", row.names = FALSE)





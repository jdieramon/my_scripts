## Build a binary matrix (presence/absence)


library(dplyr)
library(tidyr)

# Datos iniciales
df <- data.frame(
  sample = c("T-481", "T-481", "T-481", "T-482", "T-482", "T-482"),
  AG01 = c(134, 146, NA, NA, NA, 132),
  AG5  = c(148, 159, NA, 159, NA, NA)
)

# Paso 1: Convertir a formato largo
long <- pivot_longer(df, cols = -sample, names_to = "locus", values_to = "allele") %>%
  filter(!is.na(allele))

# Paso 2: Crear columna de alelo completo
long <- long %>%
  mutate(marker = paste0(locus, "_", allele),
         value = 1)

# Paso 3: Pivotear a formato binario
binary <- pivot_wider(long, names_from = marker, values_from = value, values_fill = 0)

# Resultado
print(binary)

# ------------------- Cargar librerías -------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(FSA)
library(rcompanion)

# ------------------- Establecer directorio de trabajo -------------------
setwd("PATH_LOCAL_AQUÍ")

# ------------------- Leer y procesar neuroblastoma_df (E-GEOD-73518) -------------------
probes_file <- "PATH/pruebas.txt" #Archivo TXT con las sondas pertenecientes a las ilsla CpG de RASSF1A y RASSF1C
sdrf_file <- "PATH/E-GEOD-73518.sdrf.txt"
methylation_dir <- "PATH_DIRECTORIO_ARCHIVOS_METILACIÓN"

# Leer probes de interés
probes_info <- read.table(probes_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
probes_of_interest <- probes_info$V5

# Leer metadata SDRF
sdrf <- read.table(sdrf_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
sdrf$Archivo <- tools::file_path_sans_ext(sdrf$`Derived Array Data File`)
sdrf_metadata <- sdrf[, c("Archivo", 
                          "Characteristics [age at diagnosis]", 
                          "Characteristics [11q status]", 
                          "Characteristics [mycn status]",
                          "Characteristics [organism part]",
                          "FactorValue [current risk category]",
                          "FactorValue [age at diagnosis]")]

# Leer archivos de metilación
setwd(methylation_dir)
methylation_files <- list.files(pattern = "*.txt")
methylation_data_list <- lapply(methylation_files, function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  required_cols <- c("Reporter.Identifier", "VALUE")
  if (!all(required_cols %in% colnames(data))) return(NULL)
  probe_data <- data[data$Reporter.Identifier %in% probes_of_interest, ]
  if (nrow(probe_data) > 0) {
    probe_data$Archivo <- tools::file_path_sans_ext(file)
    return(probe_data[, c("Reporter.Identifier", "VALUE", "Archivo")])
  } else {
    return(NULL)
  }
})
methylation_data <- do.call(rbind, methylation_data_list[!sapply(methylation_data_list, is.null)])
methylation_data <- merge(methylation_data, sdrf_metadata, by = "Archivo", all.x = TRUE)
neuroblastoma_df <- methylation_data %>%
  filter(Reporter.Identifier %in% probes_of_interest) %>%
  dplyr::select(
    probe = Reporter.Identifier,
    beta = VALUE,
    grupo_muestra = `FactorValue [current risk category]`,
    status_11q = `Characteristics [11q status]`,
    status_mycn = `Characteristics [mycn status]`
  )

# ------------------- Leer y procesar tejido_adrenal_df (TARGET) -------------------
tejido_dir <- "PATH_ARCHIVOS_COHORTE_TARGET"
setwd(tejido_dir)
archivos_txt <- list.files(path = tejido_dir, pattern = "\\.txt$", full.names = TRUE)
lista_datos <- lapply(archivos_txt, function(archivo) {
  datos <- read.table(archivo, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(datos) <- c("probe", "beta")
  datos$archivo <- tools::file_path_sans_ext(basename(archivo))
  return(datos)
})
datos_combinados_control <- do.call(rbind, lista_datos)
datos_filtrados <- subset(datos_combinados_control, probe %in% probes_of_interest)
datos_filtrados$beta <- as.numeric(datos_filtrados$beta)
tejido_adrenal_df <- datos_filtrados %>%
  dplyr::select(probe, beta) %>%
  dplyr::mutate(
    grupo_muestra = "tejido_adrenal",
    status_11q = NA_character_,
    status_mycn = NA_character_
  )

# ------------------- Unificar y trabajar en tu directorio de análisis -------------------
setwd("PATH_DIRECTORIO_TRABAJO")
datos_combinados <- bind_rows(neuroblastoma_df, tejido_adrenal_df)
write.csv(datos_combinados, "datos_combinados.csv", row.names = FALSE)

# ------------------- Crear variables de grupo -------------------
datos_combinados <- datos_combinados %>%
  mutate(
    grupo_plot = case_when(
      grupo_muestra == "tejido_adrenal" ~ "Control",
      status_mycn == "MYCN-amplified" ~ "MYCN amplified",
      status_mycn == "MYCN-nonamplified" ~ "MYCN nonamplified",
      TRUE ~ NA_character_
    ),
    grupo_interaccion = case_when(
      status_11q == "11q-deletion" & status_mycn == "MYCN-amplified" ~ "11q del + MYCN amplified",
      status_11q == "11q-deletion" & status_mycn == "MYCN-nonamplified" ~ "11q del + MYCN nonamplified",
      status_11q == "11q-normal" & status_mycn == "MYCN-amplified" ~ "11q normal + MYCN amplified",
      status_11q == "11q-normal" & status_mycn == "MYCN-nonamplified" ~ "11q normal + MYCN nonamplified",
      status_11q == "11q-whole loss" & status_mycn == "MYCN-amplified" ~ "11q whole loss + MYCN amplified",
      status_11q == "11q-whole loss" & status_mycn == "MYCN-nonamplified" ~ "11q whole loss + MYCN nonamplified",
      status_11q == "11q-whole gain" & status_mycn == "MYCN-amplified" ~ "11q whole gain + MYCN amplified",
      status_11q == "11q-whole gain" & status_mycn == "MYCN-nonamplified" ~ "11q whole gain + MYCN nonamplified",
      grupo_muestra == "tejido_adrenal" ~ "control",
      TRUE ~ "Otros"
    ),
    grupo_simple = ifelse(grupo_muestra == "tejido_adrenal", "Control", "Tumor")
  )



# ------------------- Comparación Control vs Tumor (Wilcoxon) -------------------
probes <- unique(datos_combinados$probe)

test_control_vs_tumor <- function(probe_id) {
  df <- datos_combinados %>% filter(probe == probe_id)
  if(length(unique(df$grupo_simple)) < 2) return(data.frame(probe = probe_id, p_val = NA))
  test_res <- tryCatch({ wilcox.test(beta ~ grupo_simple, data = df) }, error = function(e) NULL)
  if (is.null(test_res)) return(data.frame(probe = probe_id, p_val = NA))
  data.frame(probe = probe_id, p_val = test_res$p.value)
}

result_control_vs_tumor <- map_df(probes, test_control_vs_tumor) %>%
  filter(!is.na(p_val) & p_val < 0.05)

print(result_control_vs_tumor)

# ------------------- Clasificación por metilación (Control vs Tumor) -------------------
dir.create("boxplots_agrupados_pval2", showWarnings = FALSE)

significant_probes_df <- result_control_vs_tumor %>%
  filter(p_val < 0.05)

control_higher_pval <- list()
tumor_higher_pval <- list()

for (row in 1:nrow(significant_probes_df)) {
  p <- significant_probes_df$probe[row]
  pval <- significant_probes_df$p_val[row]
  
  if (is.na(p)) next  # <- salta si el valor es NA
  
  if (!(p %in% datos_combinados$probe)) next
  
  df_probe <- datos_combinados %>%
    filter(probe == p, !is.na(grupo_simple))
  
  if (nrow(df_probe) == 0) next
  
  if (n_distinct(df_probe$grupo_simple) == 2) {
    mean_beta <- df_probe %>%
      group_by(grupo_simple) %>%
      summarise(mean_beta = mean(beta, na.rm = TRUE)) %>%
      pivot_wider(names_from = grupo_simple, values_from = mean_beta)
    
    if (!is.na(mean_beta$Control) && !is.na(mean_beta$Tumor)) {
      if (mean_beta$Control > mean_beta$Tumor) {
        control_higher_pval[[length(control_higher_pval) + 1]] <- list(probe = p, pval = pval)
      } else if (mean_beta$Tumor > mean_beta$Control) {
        tumor_higher_pval[[length(tumor_higher_pval) + 1]] <- list(probe = p, pval = pval)
      }
    }
  }
}



order_by_pval <- function(probe_list) {
  if (length(probe_list) > 0) {
    print(str(probe_list))
    return(probe_list[order(sapply(probe_list, function(x) x$pval))])
  } else {
    return(list())
  }
}

control_higher_pval_ordered <- order_by_pval(control_higher_pval)
tumor_higher_pval_ordered <- order_by_pval(tumor_higher_pval)

generate_multi_boxplot_ordered <- function(probe_pval_list, filename, main_title, num_cols = 3) {
  if (length(probe_pval_list) > 0) {
    png(filename = paste0("boxplots_agrupados_pval2/", filename, ".png"),
        width = 1800, height = 1000 * ceiling(length(probe_pval_list) / num_cols), res = 300)
    par(mfrow = c(ceiling(length(probe_pval_list) / num_cols), num_cols), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 3, 0))
    
    for (item in probe_pval_list) {
      p <- item$probe
      pval <- signif(item$pval, 3)
      stars <- if (item$pval < 0.001) "***" else if (item$pval < 0.01) "**" else if (item$pval < 0.05) "*" else ""
      
      df_probe <- datos_combinados %>% filter(probe == p, !is.na(grupo_simple))
      
      boxplot(beta ~ grupo_simple, data = df_probe,
              main = paste("Probe:", p, "\np =", pval, stars),
              ylab = "Valor Beta", xlab = "Grupo",
              col = c("skyblue", "tomato"), ylim = c(0, 1),
              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1)
      abline(h = 0.3, lty = 2, col = "red")
      abline(h = 0.7, lty = 2, col = "blue")
    }
    
    mtext(main_title, outer = TRUE, cex = 1.2, line = 1)
    dev.off()
  } else {
    print(paste("No probes to plot for:", main_title))
  }
}

generate_multi_boxplot_ordered(control_higher_pval_ordered, "control_mas_metilado_pval1", "Sondas con mayor metilación en Control (ordenado por p-valor)", 3)
generate_multi_boxplot_ordered(tumor_higher_pval_ordered, "tumor_mas_metilado_pval2", "Sondas con mayor metilación en Tumor (ordenado por p-valor)", 3)



# ---- Guarda la tabla resumen ----
write.csv(summary_table, "boxplots_agrupados_pval2/tabla_sondas_graficadas.csv", row.names = FALSE)


# ------------------- Análisis por grupo de riesgo (Kruskal-Wallis) -------------------
datos_combinados <- datos_combinados %>%
  mutate(grupo_muestra = factor(grupo_muestra, levels = c("low-risk", "intermediate-risk", "high-risk")))

test_riesgo <- function(probe_id) {
  df <- datos_combinados %>% filter(probe == probe_id)
  if (length(unique(df$grupo_muestra)) < 2) return(data.frame(probe = probe_id, p_val = NA))
  test_res <- tryCatch({ kruskal.test(beta ~ grupo_muestra, data = df) }, error = function(e) NULL)
  if (is.null(test_res)) return(data.frame(probe = probe_id, p_val = NA))
  data.frame(probe = probe_id, p_val = test_res$p.value)
}

result_kruskal <- map_df(probes, test_riesgo) %>%
  filter(!is.na(p_val)) %>%
  arrange(p_val)

significant_probes_risk <- result_kruskal %>% filter(p_val < 0.05)
print(significant_probes_risk)

# ------------------- Probes que separan los 3 grupos (Dunn + Kruskal) -------------------
test_separacion_total <- function(probe_id) {
  df <- datos_combinados %>% filter(probe == probe_id, !is.na(grupo_muestra))
  if (length(unique(df$grupo_muestra)) < 3) return(data.frame(probe = probe_id, separa_3 = FALSE, p_val = NA))
  kw <- tryCatch({ kruskal.test(beta ~ grupo_muestra, data = df) }, error = function(e) NULL)
  if (is.null(kw) || kw$p.value >= 0.05) return(data.frame(probe = probe_id, separa_3 = FALSE, p_val = NA))
  dunn <- tryCatch({ dunnTest(beta ~ grupo_muestra, data = df, method = "bonferroni") }, error = function(e) NULL)
  if (is.null(dunn)) return(data.frame(probe = probe_id, separa_3 = FALSE, p_val = kw$p.value))
  all_sig <- all(dunn$res$P.adj < 0.05)
  data.frame(probe = probe_id, separa_3 = all_sig, p_val = kw$p.value)
}

resultados_dunn <- map_df(probes, test_separacion_total)

probes_separan_tres <- resultados_dunn %>% filter(separa_3 == TRUE) %>% arrange(p_val)

dir.create("boxplots_3riesgos", showWarnings = FALSE)

library(ggpubr)  # para las anotaciones estadísticas

library(dplyr)
library(ggplot2)
library(ggpubr)


comparaciones <- list(
  c("low-risk", "intermediate-risk"),
  c("low-risk", "high-risk"),
  c("intermediate-risk", "high-risk")
)

dir.create("boxplots_3riesgos", showWarnings = FALSE)

for (p in probes_separan_tres$probe) {
  df <- datos_combinados %>% filter(probe == p, !is.na(grupo_muestra))
  p_val <- signif(probes_separan_tres %>% filter(probe == p) %>% pull(p_val), 3)
  
  g <- ggplot(df, aes(x = grupo_muestra, y = beta, fill = grupo_muestra)) +
    geom_boxplot() +
    labs(title = paste0("Probe: ", p, "\nKruskal-Wallis p = ", p_val),
         x = "Grupo de riesgo", y = "Valor Beta") +
    scale_fill_manual(values = c("skyblue", "gold", "tomato")) +
    theme(plot.title = element_text(size = 10), legend.position = "none") +
    stat_compare_means(comparisons = comparaciones, method = "wilcox.test", label = "p.signif")
  
  ggsave(filename = paste0("boxplots_3riesgos/boxplot_", p, ".png"),
         plot = g, width = 5, height = 4, dpi = 300)
}




# Análisis de metilación y expresión de RASSF1 en neuroblastoma

Repositorio asociado al Trabajo Fin de Máster: “RASSF1A en neuroblastoma: caracterización transcripcional, epigenética y funcional”.

## Contenido

### 1. `Script análisis metilación neuroblastomas y cohorte TARGET.R`

- **Función**:  
  Analiza los datos de arrays de metilación (Illumina) de neuroblastoma (E-GEOD-73518) y tejido adrenal normal (TARGET), específicamente para sondas de las islas CpG de RASSF1A y RASSF1C.
- **Incluye**:
  - Lectura y filtrado de datos de metilación y metadatos clínicos.
  - Unión de cohortes tumor/control.
  - Comparaciones estadísticas (Wilcoxon, Kruskal-Wallis, Dunn).
  - Generación de boxplots y resúmenes por grupo de riesgo.
- **Notas**:  
  Es necesario ajustar las rutas (`PATH_*`) a los archivos y carpetas correspondientes de tu entorno local antes de ejecutar.

### 2. `Script scRNA seq NBATLAS.R`

- **Función**:  
  Análisis de la expresión de RASSF1 en datos de scRNA-seq de neuroblastoma (NBAtlas) y estudio de comunicación celular diferencial mediante CellChat.
- **Incluye**:
  - Visualización de expresión de RASSF1 a nivel de célula única.
  - Comparación entre subgrupos de células tumorales según expresión de RASSF1 ("high"/"low").
  - Análisis de ciclo celular y distribución de expresión por grupos de riesgo.
  - Análisis expresión diferencial de genes según expresión de RASSF1 ("high"/"low").
  - Análisis de comunicación celular (CellChat) entre grupos con distinta expresión de RASSF1.
- **Notas**:  
  Requiere instalar los paquetes indicados. Las rutas de los datos deben ser modificadas según tu ubicación local.

## Requisitos

- R >= 4.0
- Paquetes CRAN: `dplyr`, `tidyr`, `ggplot2`, `purrr`, `FSA`, `rcompanion`
- Paquetes Bioconductor/otros: `Seurat`, `CellChatV2`, `biomaRt`, `scran`, `SingleCellExperiment`, `Matrix`, `ggridges`, `patchwork`, `ggpubr`
- Archivos de datos originales (no incluidos en este repositorio, por motivos de tamaño y privacidad).

## Instrucciones de uso

1. Clona o descarga este repositorio.
2. Descarga los datos originales de los repositorios públicos correspondientes o solicita acceso a los mismos.
3. Modifica los valores de las rutas (`PATH_*`) al inicio de cada script para adaptarlos a tu sistema local.
4. Instala los paquetes necesarios siguiendo las instrucciones del propio script.
5. Ejecuta cada script en el orden adecuado para reproducir los análisis y figuras.

## Contacto

Si tienes dudas sobre el análisis o la ejecución de los scripts, puedes abrir una issue en el repositorio o contactar con el autor/a.

---

Trabajo Fin de Máster  
Raúl Golfe  
Máster biotecnología Biomédica  
Año 2025

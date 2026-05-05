# TFM-Bioinformatica
Este repositorio contiene el código completo y reproducible desarrollado para el Trabajo de Fin de Máster en Bioinformática, cuyo objetivo es la **identificación y validación de biomarcadores génicos en pacientes con ictus isquémico** mediante análisis de expresión diferencial y modelos predictivos de machine learning.

El análisis compara perfiles de expresión génica entre pacientes con ictus isquémico y controles sanos, con el fin de identificar genes candidatos con valor diagnóstico o pronóstico.

---

## Estructura del repositorio

```
TFM-Bioinformatica/
│
├── modelos_predictivos.R     # Script principal con todo el pipeline de análisis
└── README.md                 # Documentación del proyecto
```
## Pipeline de análisis

El script `modelos_predictivos.R` implementa el siguiente flujo de trabajo:

1. **Descarga y preprocesamiento de datos**
   - Obtención de datos de expresión génica desde bases de datos públicas (GEO / ArrayExpress)
   - Control de calidad y normalización de las muestras

2. **Análisis de expresión diferencial (DEGs)**
   - Identificación de genes diferencialmente expresados entre ictus isquémico y controles
   - Herramientas: `limma`, `DESeq2` / `edgeR`

3. **Análisis de enriquecimiento funcional**
   - Enriquecimiento en términos GO (Gene Ontology) y rutas KEGG
   - Interpretación biológica de los genes significativos

4. **Construcción y comparativa de modelos predictivos**
   - Entrenamiento de modelos de machine learning para clasificación de muestras
   - Herramientas: `caret`, `randomForest`, `glmnet` (regresión regularizada / LASSO)
   - Comparativa del rendimiento entre modelos (AUC, accuracy, sensibilidad, especificidad)

5. **Visualización de resultados**
   - Volcano plots, heatmaps, curvas ROC y gráficas de comparativa entre modelos

---

## Requisitos

### Versión de R
R ≥ 4.0.0 recomendado

### Paquetes necesarios

```r
install.packages(c("caret", "randomForest", "glmnet", "clusterProfiler" "EnrichR))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "DESeq2", "edgeR", "clusterProfiler", "org.Hs.eg.db"))
```

---

## Datos

Los datos de expresión génica utilizados provienen de repositorios públicos:

- **GEO** (Gene Expression Omnibus): [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/)


>  ID de acceso para  la base de datos publicas (ej. `GSE22255`)

---

## Cómo reproducir el análisis

1. Clona el repositorio:
   ```bash
   git clone https://github.com/esthernom/TFM-Bioinformatica.git
   cd TFM-Bioinformatica
   ```

2. Instala los paquetes necesarios (ver sección anterior)

3. Ejecuta el script principal en R o RStudio:
   ```r
   source("modelos_predictivos.R")
   ```

---

## Autora

**Esther Nombela ** — Máster en Bioinformática  
Repositorio disponible en: [https://github.com/esthernom/TFM-Bioinformatica](https://github.com/esthernom/TFM-Bioinformatica)

---

## Licencia

Este proyecto se comparte con fines académicos.

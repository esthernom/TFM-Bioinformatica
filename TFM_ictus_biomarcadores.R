#TFM PARA BUSCAR BIOMACARCADORES EN EL IS

rm(list = ls())

library(GEOquery) # buscar la base de datos
library(tidyverse)
library(limma) #buscar los Genes diferencialmente expresados
library(glmnet) #LASSO
library(pROC)
library(ComplexHeatmap) #graficos(heatmap) 
library(circlize)
library(stringr) #para hacer las listas 
library(randomForest) #Random Forest (RF)
library(e1071)      # SVM
library(xgboost)    # Gradient Boosting
library(caret)      # Matrices de confusión

set.seed(123) #SEMILLA PARA REPRODUCTIBILIDAD


# 1. CARGA DE DATOS Y PREPROCESADO


gse <- getGEO("GSE22255", GSEMatrix = TRUE)
if (is.list(gse)) gse <- gse[[1]]

exp_matrix   <- exprs(gse)   # matriz expresión (sondas x muestras)
feature.Data <- fData(gse)   # anotación sondas
pdata        <- pData(gse)   # para ver los datos que tenemos y uqe corresponde a cada columna 

# Mapa ID sonda --> para pasarlo a simbolos de los genes 
feature.Data <- feature.Data[, c(1, 11)]
colnames(feature.Data) <- c("ID", "Gene.Symbol")
gene_symbols_map <- feature.Data

# Variables fenotípicas
condicion <- factor(ifelse(grepl("IS", pdata$characteristics_ch1),
                                  "IS", "Control"))

sexo <- factor(ifelse(grepl("female", pdata$characteristics_ch1.1,
                                   ignore.case = TRUE),
                             "female", "male"))

edad <- as.numeric(gsub("age-at-examination: ", "",
                                 pdata$characteristics_ch1.2))

hipertension <- factor(ifelse(grepl("Hypertension", pdata$characteristics_ch1.5,
                                      ignore.case = TRUE),
                                "Yes", "No"))

############################################################
# 2.LIMMA SIMPLE: CONTROL vs IS (sin los cofactores de riesgo )
limma_simple <- model.matrix(~ 0 + condicion)
colnames(limma_simple) <- levels(condicion)

fit_simple <- lmFit(exp_matrix, limma_simple)
cont.matrix <- makeContrasts(IS_vs_Control = IS - Control,
                             levels = limma_simple)
fit2_simple <- contrasts.fit(fit_simple, cont.matrix)
fit2_simple <- eBayes(fit2_simple)

results_simple <- topTable(fit2_simple,
                           coef = "IS_vs_Control",
                           number = Inf) %>%
  rownames_to_column("ID") %>%
  inner_join(gene_symbols_map, by = "ID") %>%
  filter(!is.na(Gene.Symbol) & Gene.Symbol != "" & Gene.Symbol != "---")

genes_ids <- results_simple %>%
  filter(P.Value < 0.05) %>%
  pull(ID)

cat("DEGs  (IS vs Control):", length(genes_ids), "\n")

###############
# 3.MATRIZ PARA LASSO Y MODELOS (CONTROLvsIS)
#######

X_full <- t(exp_matrix[genes_ids, ])   # muestras x genes
Y_full <- condicion          # factor: Control / IS

# División 70/30
set.seed(123)
idx_train <- sample(1:nrow(X_full), 0.7 * nrow(X_full))

X_train <- X_full[idx_train, ]
X_test  <- X_full[-idx_train, ]

Y_train <- Y_full[idx_train]
Y_test  <- Y_full[-idx_train]

######
# 4. LASSO (FIRMA DE 15 GENES: CONTROL vs IS)
# Firma principal de biomarcadores para modelos
###

set.seed(123)
cv_lasso <- cv.glmnet(
  X_train, Y_train,
  family = "binomial",
  alpha = 1
)

coefs <- coef(cv_lasso, s = "lambda.min")
coefs_mat <- as.matrix(coefs)

genes_firma_ids <- rownames(coefs_mat)[coefs_mat != 0 &
                                         rownames(coefs_mat) != "(Intercept)"]

# Si hay más de 15 genes, nos quedamos con los 15 de mayor |coef|
if (length(genes_firma_ids) > 15) {
  coefs_sub <- coefs_mat[genes_firma_ids, , drop = FALSE]
  ord <- order(abs(coefs_sub[, 1]), decreasing = TRUE)
  genes_firma_ids <- genes_firma_ids[ord][1:15]
}

cat("Genes seleccionados por LASSO (firma):", length(genes_firma_ids), "\n")
print(genes_firma_ids)

# Pesos de la firma (coeficientes LASSO)
pesos_lasso <- coefs_mat[genes_firma_ids, , drop = FALSE]

###########
# 5. HEATMAP LASSO (CONTROL vs IS)
############

mat_firma <- exp_matrix[genes_firma_ids, ]

rownames(mat_firma) <- str_split_i(
  gene_symbols_map$Gene.Symbol[match(genes_firma_ids, gene_symbols_map$ID)],
  " /// ", 1
)

mat_firma_scaled <- t(scale(t(mat_firma)))

orden_muestras <- order(condicion)
mat_firma_ord  <- mat_firma_scaled[, orden_muestras]

estado_ordenado <- condicion[orden_muestras]

ha_firma <- HeatmapAnnotation(
  Estado = estado_ordenado,
  col = list(
    Estado = c("Control" = "#4DA3FF", "IS" = "#FF8C42")
  )
)

Heatmap(
  mat_firma_ord,
  name = "Z-score",
  top_annotation = ha_firma,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_title = "HEATMAP Control vs IS/ Genes seleccionados con LASSO",
  border = TRUE
)

# 6.DATA.FRAMES PARA LOS 4 MODELOS PREDICTIVOS
#    (Todos usan ControlvsIS y cpn los 15 genes LASSO, con todos colapsarían )
#######

X_train_firma <- X_train[, genes_firma_ids]
X_test_firma  <- X_test[, genes_firma_ids]

df_train <- data.frame(X_train_firma, Estado = Y_train)
df_test  <- data.frame(X_test_firma,  Estado = Y_test)

#####
# 7.MODELO 1: LASSO SCORE Y PROBABILIDAD EN TEST
####

score_lasso_test <- as.matrix(X_test_firma) %*% pesos_lasso
prob_lasso_test  <- plogis(score_lasso_test)  # transformar a probabilidad

df_score_lasso <- data.frame(
  Score  = as.numeric(score_lasso_test),
  Prob   = as.numeric(prob_lasso_test),
  Estado = Y_test
)

roc_lasso <- roc(df_score_lasso$Estado, df_score_lasso$Prob,
                 levels = c("Control", "IS"),
                 direction = ">")

#########
# 8.MODELO 2: RANDOM FOREST
######

set.seed(123)
rf_model <- randomForest(Estado ~ ., data = df_train,
                         ntree = 500, importance = TRUE)

rf_probs <- predict(rf_model, newdata = df_test, type = "prob")[, "IS"]

roc_rf <- roc(df_test$Estado, rf_probs,
              levels = c("Control", "IS"),
              direction = ">")

###########
# 9. MODELO 3: SVM (KERNEL RADIAL)
##########

set.seed(123)
svm_model <- svm(
  Estado ~ .,
  data = df_train,
  kernel = "radial",
  probability = TRUE
)

svm_pred <- predict(svm_model, newdata = df_test, probability = TRUE)
svm_probs <- attr(svm_pred, "probabilities")[, "IS"]

roc_svm <- roc(df_test$Estado, svm_probs,
               levels = c("Control", "IS"),
               direction = ">")

#####
# 10. MODELO 4: XGBOOST (GRADIENT BOOSTING)
########

y_train_num <- ifelse(df_train$Estado == "IS", 1, 0)
y_test_num  <- ifelse(df_test$Estado == "IS", 1, 0)

dtrain <- xgb.DMatrix(data = as.matrix(X_train_firma), label = y_train_num)
dtest  <- xgb.DMatrix(data = as.matrix(X_test_firma),  label = y_test_num)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 200,
  verbose = 0
)

xgb_probs <- predict(xgb_model, dtest)

roc_xgb <- roc(df_test$Estado, xgb_probs,
               levels = c("Control", "IS"),
               direction = ">")

#########
# 11. CURVAS ROC INDIVIDUALES (por cada modelos)
#########

# LASSO
plot(roc_lasso,
     col = "#2c7bb6",
     lwd = 3,
     main = "Curva ROC – LASSO",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "#f0f0f0",
     grid = TRUE)

# Random Forest
plot(roc_rf,
     col = "forestgreen",
     lwd = 3,
     main = "Curva ROC – Random Forest",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "#e0ffe0",
     grid = TRUE)

# SVM
plot(roc_svm,
     col = "darkorange",
     lwd = 3,
     main = "Curva ROC – SVM radial",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "#ffe0c0",
     grid = TRUE)

# XGBoost
plot(roc_xgb,
     col = "firebrick",
     lwd = 3,
     main = "Curva ROC – XGBoost",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "#ffd6d6",
     grid = TRUE)

#######
# 12. CURVA ROC COMPARATIVA
#######

plot(roc_lasso,
     col = "#2c7bb6",
     lwd = 3,
     main = "Comparación de modelos (15 genes LASSO)",
     print.auc = FALSE)

plot(roc_rf,
     col = "forestgreen",
     lwd = 3,
     add = TRUE)

plot(roc_svm,
     col = "darkorange",
     lwd = 3,
     add = TRUE)

plot(roc_xgb,
     col = "firebrick",
     lwd = 3,
     add = TRUE)

legend("bottomright",
       legend = c(
         paste0("LASSO (AUC = ", round(auc(roc_lasso), 3), ")"),
         paste0("Random Forest (AUC = ", round(auc(roc_rf), 3), ")"),
         paste0("SVM (AUC = ", round(auc(roc_svm), 3), ")"),
         paste0("XGBoost (AUC = ", round(auc(roc_xgb), 3), ")")
       ),
       col = c("#2c7bb6", "forestgreen", "darkorange", "firebrick"),
       lty = 1, lwd = 3)

#########
# 13. MATRICES DE CONFUSIÓN PARA LOS 4 MODELOS
#######

# 1) LASSO
lasso_pred_class <- ifelse(df_score_lasso$Prob > 0.5, "IS", "Control")
lasso_pred_class <- factor(lasso_pred_class, levels = c("Control", "IS"))

cat("\nMATRIZ DE CONFUSIÓN: LASSO\n")
print(confusionMatrix(lasso_pred_class, Y_test))

# 2) Random Forest
rf_pred_class <- ifelse(rf_probs > 0.5, "IS", "Control")
rf_pred_class <- factor(rf_pred_class, levels = c("Control", "IS"))

cat("\nMATRIZ DE CONFUSIÓN: RANDOM FOREST\n")
print(confusionMatrix(rf_pred_class, Y_test))

# 3) SVM
svm_pred_class <- ifelse(svm_probs > 0.5, "IS", "Control")
svm_pred_class <- factor(svm_pred_class, levels = c("Control", "IS"))

cat("\nMATRIZ DE CONFUSIÓN: SVM\n")
print(confusionMatrix(svm_pred_class, Y_test))

# 4) XGBoost
xgb_pred_class <- ifelse(xgb_probs > 0.5, "IS", "Control")
xgb_pred_class <- factor(xgb_pred_class, levels = c("Control", "IS"))

cat("\n--- MATRIZ DE CONFUSIÓN: XGBOOST ---\n")
print(confusionMatrix(xgb_pred_class, Y_test))

######### 14.GENES DIFERENCIALES AJUSTADOS(LIMMA)
limma_cofactores <- model.matrix(~ condicion + sexo + edad + hipertension)
colnames(limma_cofactores)

fit_cof <- lmFit(exp_matrix, limma_cofactores)
fit2_cof <- eBayes(fit_cof)

results_cof <- topTable(fit2_cof,
                        coef = "condicionIS",   # nombre correcto
                        number = Inf) %>%
  rownames_to_column("ID") %>%
  inner_join(gene_symbols_map, by="ID") %>%
  filter(Gene.Symbol != "" & Gene.Symbol != "---")



######
#SELECCIÓN DE GENES AJUSTADOS (UMBRAL ROBUSTO) ##COFACTORES
###

simbolos_ids_cof <- results_cof %>%
  filter(P.Value < 0.01) %>%     # antes 0.005 → demasiado estricto
  pull(ID)

cat("Genes ajustados significativos (P<0.01):", length(simbolos_ids_cof), "\n")

# Si hay pocos genes, ampliamos automáticamente
if (length(simbolos_ids_cof) < 20) {
  simbolos_ids_cof <- results_cof %>%
    filter(P.Value < 0.05) %>%
    pull(ID)
  cat("Umbral ampliado a P < 0.05. Genes:", length(sig_ids_adj), "\n")
}

#########
#  LASSO AJUSTADO 
#############

X_cof <- t(exp_matrix[simbolos_ids_cof, ])
Y_cof <- condicion

if (ncol(X_cof) < 5) {
  stop("No hay suficientes genes significativos tras el ajuste para ejecutar LASSO.")
}

cv_lasso_cof <- cv.glmnet(X_cof, Y_cof, family="binomial", alpha=1)

coefs_cof <- coef(cv_lasso_cof, s="lambda.min")
coefs_cof_mat <- as.matrix(coefs_cof)

genes_cof <- rownames(coefs_cof_mat)[coefs_cof_mat != 0 &
                                           rownames(coefs_cof_mat) != "(Intercept)"]

cat("Genes seleccionados por LASSO ajustado:", length(genes_cof), "\n")
print(genes_cof)

##########
# HEATMAP LASSO AJUSTADO + COFACTORES
##########

mat_cof <- exp_matrix[genes_cof, ]
rownames(mat_cof) <- str_split_i(
  gene_symbols_map$Gene.Symbol[match(genes_cof, gene_symbols_map$ID)],
  " /// ", 1
)

mat_adj_scaled <- t(scale(t(mat_cof)))
orden <- order(condicion)

Heatmap(
  mat_adj_scaled[, orden],
  name="Z-score",
  top_annotation = HeatmapAnnotation(
    Estado = condicion[orden],
    Sexo   = sexo[orden],
    Edad   = edad[orden],
    HTA    = hipertension[orden],
    col = list(
      Estado = c("Control"="#4DA3FF","IS"="#FF8C42"),
      Sexo   = c("male"="#6A5ACD","female"="#FF69B4"),
      HTA    = c("Yes"="#E41A1C","No"="#4DAF4A"),
      Edad   = colorRamp2(
        c(min(edad), median(edad), max(edad)),
        c("green3","yellow","firebrick3")
      )
    )
  ),
  col=colorRamp2(c(-2,0,2),c("navy","white","firebrick3")),
  cluster_rows=TRUE,
  cluster_columns=FALSE,
  show_column_names=FALSE,
  column_title="LASSO ajustado (sexo, edad, HTA)"
)


#########
# PREPARAR LISTAS DE GENES PARA ENRIQUECIMIENTO
#########

# 1. DEGs ajustados 
genes_deg_cof <- results_cof %>%
  filter(P.Value < 0.05) %>%
  pull(Gene.Symbol) %>%
  str_split_i(" /// ", 1) %>%
  unique()
#con esto lo que vemos son los genes diferenciales ajustados por los cofactores
#tener una visión global, y hacer los analisis funcioanles amplios
# 2. Firma LASSO simple
genes_lasso_simple <- str_split_i(
  gene_symbols_map$Gene.Symbol[match(genes_firma_ids, gene_symbols_map$ID)],
  " /// ", 1
)
#los 15 genes seleccionados con lasso para ver el valos discriminante 
# 3. Firma LASSO ajustada
genes_lasso_adj <- str_split_i(
  gene_symbols_map$Gene.Symbol[match(genes_cof, gene_symbols_map$ID)],
  " /// ", 1
)
## análisis funcional con enrich 
library(enrichR)

dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Reactome_2022"
)

enrich_deg_adj <- enrichr(genes_deg_cof, dbs)
enrich_lasso_adj <- enrichr(genes_lasso_adj, dbs)
##visualizamos GO
plotEnrich(
  enrich_deg_adj[["GO_Biological_Process_2023"]],
  showTerms = 15,
  numChar = 45,
  y = "Count",
  orderBy = "P.value",
  title = "Procesos biológicos (GO)"
)
##VISUALIZAR CON KEGGS
plotEnrich(
  enrich_deg_adj[["KEGG_2021_Human"]],
  showTerms = 15,
  numChar = 45,
  y = "Count",
  orderBy = "P.value",
  title = "Procesos biológicos (KEGG)")


#analisis con clusterprofile 
library(clusterProfiler)
library(org.Hs.eg.db)

genes_entrez <- bitr(genes_deg_cof,
                     fromType="SYMBOL",
                     toType="ENTREZID",
                     OrgDb="org.Hs.eg.db")
## GO
ego <- enrichGO(
  gene = genes_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
)
dotplot(ego, showCategory = 20) + ggtitle("Procesos biológicos(GO)")

##kegg
ekegg <- enrichKEGG(
  gene = genes_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

dotplot(ekegg, showCategory = 20) + ggtitle("Procesos biológicos enriquecidos (KEGG)")
###Hasta aquí se hae el analisis funcional con cluster profiler con los DEGS ajustados
#procesos globales estamos observando 

###LASSO
#aqui vemos tambien los degs ajustados, pero tras ser seleccionados con LASSO
#que funciones comparten los genes que mejor discriminan entre ISvs Control 
genes_lasso_entrez <- bitr(genes_lasso_adj,
                           fromType="SYMBOL",
                           toType="ENTREZID",
                           OrgDb="org.Hs.eg.db")

ego_lasso <- enrichGO(
  gene = genes_lasso_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  readable = TRUE
)

dotplot(ego_lasso, showCategory = 15) +  ggtitle("Procesos biológicos enriquecidos (LASSO)")
##string por si queremos tener una lsta con todos los genes 
write.table(
  genes_deg_cof,
  "DEGs_ajustados_para_STRING.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

####### COMO SVM es el que mejor elige hacemos aqui para ver como salen un heatmap con los genes seleccionados 
#  ELECCIÓN DE GENES CON SVM LINEAL 
########

# Usamos solo los DEGs simples (IS vs Control)
X_svm <- t(exp_matrix[genes_ids, ])
Y_svm <- condicion

# Entrenamos SVM lineal (permite extraer pesos)
set.seed(123)
svm_linear <- svm(
  X_svm, Y_svm,
  kernel = "linear",
  scale = TRUE
)

# Extraer pesos del modelo SVM
w <- t(svm_linear$coefs) %*% svm_linear$SV
w <- as.numeric(w)

names(w) <- genes_ids

# Ordenar por importancia absoluta
w_ord <- sort(abs(w), decreasing = TRUE)

# Seleccionar los 15 genes más importantes
genes_svm_ids <- names(w_ord)[1:15]

cat("Genes seleccionados por SVM (firma):", length(genes_svm_ids), "\n")
print(genes_svm_ids)

# Obtener símbolos de genes
genes_svm_symbols <- str_split_i(
  gene_symbols_map$Gene.Symbol[match(genes_svm_ids, gene_symbols_map$ID)],
  " /// ", 1
)

#######
# 14.HEATMAP SVM (solo Control vs IS)
######

mat_svm <- exp_matrix[genes_svm_ids, ]

rownames(mat_svm) <- genes_svm_symbols

# Escalado por gen
mat_svm_scaled <- t(scale(t(mat_svm)))

# Ordenar muestras por condición
orden_svm <- order(condicion)
mat_svm_ord <- mat_svm_scaled[, orden_svm]

estado_ordenado_svm <- condicion[orden_svm]

ha_svm <- HeatmapAnnotation(
  Estado = estado_ordenado_svm,
  col = list(
    Estado = c("Control" = "#4DA3FF", "IS" = "#FF8C42")
  )
)

Heatmap(
  mat_svm_ord,
  name = "Z-score",
  top_annotation = ha_svm,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_title = "HEATMAP Control vs IS / Genes seleccionados con SVM",
  border = TRUE
)

##### MODELOS CON LOS 4070 genes

genes_ids <- results_simple %>%
  filter(P.Value < 0.05) %>%
  pull(ID)
############################################################
# MATRICES PARA MODELOS CON TODOS LOS DEGs (≈4707 genes)
############################################################

X_full_all <- t(exp_matrix[genes_ids, ])   # muestras x genes
Y_full_all <- condicion

set.seed(123)
idx_train_all <- sample(1:nrow(X_full_all), 0.7 * nrow(X_full_all))

X_train_all <- X_full_all[idx_train_all, ]
X_test_all  <- X_full_all[-idx_train_all, ]

Y_train_all <- Y_full_all[idx_train_all]
Y_test_all  <- Y_full_all[-idx_train_all]


############################################################
# MODELO 1: LASSO con todos los DEGs
############################################################

set.seed(123)
cv_lasso_all <- cv.glmnet(
  X_train_all, Y_train_all,
  family = "binomial",
  alpha = 1
)

lasso_probs_all <- predict(cv_lasso_all, newx = X_test_all,
                           s = "lambda.min", type = "response")

roc_lasso_all <- roc(Y_test_all, as.numeric(lasso_probs_all),
                     levels = c("Control", "IS"),
                     direction = ">")
############################################################
# MODELO 2: Random Forest con todos los DEGs
############################################################

set.seed(123)
rf_model_all <- randomForest(
  x = X_train_all,
  y = Y_train_all,
  ntree = 500,
  importance = TRUE
)

rf_probs_all <- predict(rf_model_all, newdata = X_test_all, type = "prob")[, "IS"]

roc_rf_all <- roc(Y_test_all, rf_probs_all,
                  levels = c("Control", "IS"),
                  direction = ">")
###
# MODELO 3: SVM radial con todos los DEGs
###

set.seed(123)
svm_model_all <- svm(
  x = X_train_all,
  y = Y_train_all,
  kernel = "radial",
  probability = TRUE
)

svm_pred_all <- predict(svm_model_all, newdata = X_test_all, probability = TRUE)
svm_probs_all <- attr(svm_pred_all, "probabilities")[, "IS"]

roc_svm_all <- roc(Y_test_all, svm_probs_all,
                   levels = c("Control", "IS"),
                   direction = ">")
############################################################
# MODELO 4: XGBOOST con todos los DEGs
############################################################

y_train_num_all <- ifelse(Y_train_all == "IS", 1, 0)
y_test_num_all  <- ifelse(Y_test_all == "IS", 1, 0)

dtrain_all <- xgb.DMatrix(data = as.matrix(X_train_all), label = y_train_num_all)
dtest_all  <- xgb.DMatrix(data = as.matrix(X_test_all),  label = y_test_num_all)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(123)
xgb_model_all <- xgb.train(
  params = params,
  data = dtrain_all,
  nrounds = 200,
  verbose = 0
)

xgb_probs_all <- predict(xgb_model_all, dtest_all)

roc_xgb_all <- roc(Y_test_all, xgb_probs_all,
                   levels = c("Control", "IS"),
                   direction = ">")
#######
# COMPARACIÓN ROC (todos los DEGs)
###########

plot(roc_lasso_all, col = "#2c7bb6", lwd = 3,
     main = "Comparación de modelos (todos los DEGs)",
     print.auc = FALSE)

plot(roc_rf_all, col = "forestgreen", lwd = 3, add = TRUE)
plot(roc_svm_all, col = "darkorange", lwd = 3, add = TRUE)
plot(roc_xgb_all, col = "firebrick", lwd = 3, add = TRUE)

legend("bottomright",
       legend = c(
         paste0("LASSO (AUC = ", round(auc(roc_lasso_all), 3), ")"),
         paste0("Random Forest (AUC = ", round(auc(roc_rf_all), 3), ")"),
         paste0("SVM (AUC = ", round(auc(roc_svm_all), 3), ")"),
         paste0("XGBoost (AUC = ", round(auc(roc_xgb_all), 3), ")")
       ),
       col = c("#2c7bb6", "forestgreen", "darkorange", "firebrick"),
       lty = 1, lwd = 3)


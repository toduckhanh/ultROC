require(ggplot2)
require(gridExtra)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### ---- Obtaining data ----
# install.packages("BiocManager")
# BiocManager::install("GEOquery")

require(GEOquery)

gset <- getGEO("GSE19188", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
rownames(ex)
colnames(ex)
dim(ex)

table(gset@phenoData@data$`cell type:ch1`)

probe_set <- c("203914_x_at", "211548_s_at", "218087_s_at", "201578_at")

data_genes <- data.frame(
  Group = factor(gset@phenoData@data$`cell type:ch1`,
                  levels = c("healthy", "SCC", "ADC", "LCC"))
)

for(i in 1:length(probe_set)) {
  data_genes[, i + 1] <- ex[probe_set[i], ]
}

colnames(data_genes) <- c("Group", "HPGD_1", "HPGD_2", "SORBS1", "PODXL")

###---- FINAL ANALYSIS ----
n <- table(data_genes$Group)
n

### ---- plot densities ----
plot_dens_1 <- ggplot(data = data_genes, aes(x = HPGD_1, fill = Group)) +
  geom_density(color = "#e9ecef", alpha = 0.4, position = "identity") +
  theme_bw() + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                     legend.position = "bottom")

plot_dens_2 <- ggplot(data = data_genes, aes(x = HPGD_2, fill = Group)) +
  geom_density(color = "#e9ecef", alpha = 0.4, position = "identity") +
  theme_bw() + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

plot_dens_3 <- ggplot(data = data_genes, aes(x = SORBS1, fill = Group)) +
  geom_density(color = "#e9ecef", alpha = 0.4, position = "identity") +
  theme_bw() + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

plot_dens_4 <- ggplot(data = data_genes, aes(x = PODXL, fill = Group)) +
  geom_density(color = "#e9ecef", alpha = 0.4, position = "identity") +
  theme_bw() + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

mylegend <- g_legend(plot_dens_1)

grid.arrange(arrangeGrob(plot_dens_1 + theme(legend.position = "none"),
                         plot_dens_2 + theme(legend.position = "none"),
                         plot_dens_3 + theme(legend.position = "none"),
                         plot_dens_4 + theme(legend.position = "none"),
                         nrow = 2), mylegend, nrow = 2, heights = c(10, 1))

### ---- define the data for biomarkers ----
HPGD_1_1 <- (-1) * data_genes$HPGD_1[data_genes$Group == "healthy"]
HPGD_1_2 <- (-1) * data_genes$HPGD_1[data_genes$Group == "SCC"]
HPGD_1_3 <- (-1) * data_genes$HPGD_1[data_genes$Group == "ADC"]
HPGD_1_4 <- (-1) * data_genes$HPGD_1[data_genes$Group == "LCC"]
HPGD_1_ls <- list(HPGD_1_1, HPGD_1_2, HPGD_1_3, HPGD_1_4)

HPGD_2_1 <- (-1) * data_genes$HPGD_2[data_genes$Group == "healthy"]
HPGD_2_2 <- (-1) * data_genes$HPGD_2[data_genes$Group == "SCC"]
HPGD_2_3 <- (-1) * data_genes$HPGD_2[data_genes$Group == "ADC"]
HPGD_2_4 <- (-1) * data_genes$HPGD_2[data_genes$Group == "LCC"]
HPGD_2_ls <- list(HPGD_2_1, HPGD_2_2, HPGD_2_3, HPGD_2_4)

SORBS1_1 <- (-1) * data_genes$SORBS1[data_genes$Group == "healthy"]
SORBS1_2 <- (-1) * data_genes$SORBS1[data_genes$Group == "SCC"]
SORBS1_3 <- (-1) * data_genes$SORBS1[data_genes$Group == "ADC"]
SORBS1_4 <- (-1) * data_genes$SORBS1[data_genes$Group == "LCC"]
SORBS1_ls <- list(SORBS1_1, SORBS1_2, SORBS1_3, SORBS1_4)

PODXL_1 <- (-1) * data_genes$PODXL[data_genes$Group == "healthy"]
PODXL_2 <- (-1) * data_genes$PODXL[data_genes$Group == "SCC"]
PODXL_3 <- (-1) * data_genes$PODXL[data_genes$Group == "ADC"]
PODXL_4 <- (-1) * data_genes$PODXL[data_genes$Group == "LCC"]
PODXL_ls <- list(PODXL_1, PODXL_2, PODXL_3, PODXL_4)

### ---- Empirical estimation for indexes ----
## ---- HPGD_1 ----
LTAUC_HPGD_1 <- LTAUC_emp(HPGD_1_ls)
UTAUC_HPGD_1 <- UTAUC_emp(HPGD_1_ls)

res_bts_HPGD_1 <- bts_func(HPGD_1_ls, n, LTAUC_HPGD_1, UTAUC_HPGD_1, B = 200)

w1_HPGD_1 <- 1 / mean(res_bts_HPGD_1[1,])
w2_HPGD_1 <- (LTAUC_HPGD_1 * (1 - LTAUC_HPGD_1) / n[1]) /
  var(res_bts_HPGD_1[3,])

w3_HPGD_1 <- 1 / mean(res_bts_HPGD_1[2,])
w4_HPGD_1 <- (UTAUC_HPGD_1 * (1 - UTAUC_HPGD_1) / n[1]) /
  var(res_bts_HPGD_1[4,])

ci_LTAUC_HPGD_1 <- EL_ci_LUTAUC(LTAUC_HPGD_1, n = n[1], w_est = w1_HPGD_1,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(LTAUC_HPGD_1, n = n[1], w_est = w2_HPGD_1,
             ci_level = c(0.9, 0.95, 0.99))

ci_UTAUC_HPGD_1 <- EL_ci_LUTAUC(UTAUC_HPGD_1, n = n[1], w_est = w3_HPGD_1,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(UTAUC_HPGD_1, n = n[1], w_est = w4_HPGD_1,
             ci_level = c(0.9, 0.95, 0.99))


c(LTAUC_HPGD_1, ci_LTAUC_HPGD_1[[2]], UTAUC_HPGD_1, ci_UTAUC_HPGD_1[[2]])

## ---- HPGD_2 ----
LTAUC_HPGD_2 <- LTAUC_emp(HPGD_2_ls)
UTAUC_HPGD_2 <- UTAUC_emp(HPGD_2_ls)

res_bts_HPGD_2 <- bts_func(HPGD_2_ls, n, LTAUC_HPGD_2, UTAUC_HPGD_2, B = 200)

w1_HPGD_2 <- 1 / mean(res_bts_HPGD_2[1,])
w2_HPGD_2 <- (LTAUC_HPGD_2 * (1 - LTAUC_HPGD_2) / n[1]) /
  var(res_bts_HPGD_2[3,])

w3_HPGD_2 <- 1 / mean(res_bts_HPGD_2[2,])
w4_HPGD_2 <- (UTAUC_HPGD_2 * (1 - UTAUC_HPGD_2) / n[1]) /
  var(res_bts_HPGD_2[4,])

ci_LTAUC_HPGD_2 <- EL_ci_LUTAUC(LTAUC_HPGD_2, n = n[1], w_est = w1_HPGD_2,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(LTAUC_HPGD_2, n = n[1], w_est = w2_HPGD_2,
             ci_level = c(0.9, 0.95, 0.99))

ci_UTAUC_HPGD_2 <- EL_ci_LUTAUC(UTAUC_HPGD_2, n = n[1], w_est = w3_HPGD_2,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(UTAUC_HPGD_2, n = n[1], w_est = w4_HPGD_2,
             ci_level = c(0.9, 0.95, 0.99))

c(LTAUC_HPGD_2, ci_LTAUC_HPGD_2[[2]], UTAUC_HPGD_2, ci_UTAUC_HPGD_2[[2]])

## ---- SORBS1 ----
LTAUC_SORBS1 <- LTAUC_emp(SORBS1_ls)
UTAUC_SORBS1 <- UTAUC_emp(SORBS1_ls)

res_bts_SORBS1 <- bts_func(SORBS1_ls, n, LTAUC_SORBS1, UTAUC_SORBS1, B = 200)

w1_SORBS1 <- 1 / mean(res_bts_SORBS1[1,])
w2_SORBS1 <- (LTAUC_SORBS1 * (1 - LTAUC_SORBS1) / n[1]) /
  var(res_bts_SORBS1[3,])

w3_SORBS1 <- 1 / mean(res_bts_SORBS1[2,])
w4_SORBS1 <- (UTAUC_SORBS1 * (1 - UTAUC_SORBS1) / n[1]) /
  var(res_bts_SORBS1[4,])

ci_LTAUC_SORBS1 <- EL_ci_LUTAUC(LTAUC_SORBS1, n = n[1], w_est = w1_SORBS1,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(LTAUC_SORBS1, n = n[1], w_est = w2_SORBS1,
             ci_level = c(0.9, 0.95, 0.99))

ci_UTAUC_SORBS1 <- EL_ci_LUTAUC(UTAUC_SORBS1, n = n[1], w_est = w3_SORBS1,
                                ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(UTAUC_SORBS1, n = n[1], w_est = w4_SORBS1,
             ci_level = c(0.9, 0.95, 0.99))

c(LTAUC_SORBS1, ci_LTAUC_SORBS1[[2]], UTAUC_SORBS1, ci_UTAUC_SORBS1[[2]])

## ---- PODXL ----
LTAUC_PODXL <- LTAUC_emp(PODXL_ls)
UTAUC_PODXL <- UTAUC_emp(PODXL_ls)

res_bts_PODXL <- bts_func(PODXL_ls, n, LTAUC_PODXL, UTAUC_PODXL, B = 200)

w1_PODXL <- 1 / mean(res_bts_PODXL[1,])
w2_PODXL <- (LTAUC_PODXL * (1 - LTAUC_PODXL) / n[1]) /
  var(res_bts_PODXL[3,])

w3_PODXL <- 1 / mean(res_bts_PODXL[2,])
w4_PODXL <- (UTAUC_PODXL * (1 - UTAUC_PODXL) / n[1]) /
  var(res_bts_PODXL[4,])

ci_LTAUC_PODXL <- EL_ci_LUTAUC(LTAUC_PODXL, n = n[1], w_est = w1_PODXL,
                               ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(LTAUC_PODXL, n = n[1], w_est = w2_PODXL,
             ci_level = c(0.9, 0.95, 0.99))

ci_UTAUC_PODXL <- EL_ci_LUTAUC(UTAUC_PODXL, n = n[1], w_est = w3_PODXL,
                               ci_level = c(0.9, 0.95, 0.99))

EL_ci_LUTAUC(UTAUC_PODXL, n = n[1], w_est = w4_PODXL,
             ci_level = c(0.9, 0.95, 0.99))

c(LTAUC_PODXL, ci_LTAUC_PODXL[[2]], UTAUC_PODXL, ci_UTAUC_PODXL[[2]])

round(
  rbind("HPGD-1" = c(LTAUC_HPGD_1, ci_LTAUC_HPGD_1[[2]], UTAUC_HPGD_1,
                     ci_UTAUC_HPGD_1[[2]]),
        "HPGD-2" = c(LTAUC_HPGD_2, ci_LTAUC_HPGD_2[[2]], UTAUC_HPGD_2,
                     ci_UTAUC_HPGD_2[[2]]),
        "SORBS1" = c(LTAUC_SORBS1, ci_LTAUC_SORBS1[[2]], UTAUC_SORBS1,
                     ci_UTAUC_SORBS1[[2]]),
        "PODXL" = c(LTAUC_PODXL, ci_LTAUC_PODXL[[2]], UTAUC_PODXL,
                    ci_UTAUC_PODXL[[2]])), 4
)

### ---- LTROC and UTROC curves ----
p1 <- seq(0, 1, by = 0.001)
LTROC_emp_HPGD_1 <- LTROC_emp(p1, Xlist = HPGD_1_ls)
UTROC_emp_HPGD_1 <- UTROC_emp(p1, Xlist = HPGD_1_ls)

LTROC_emp_HPGD_1[p1 == 0] <- 0
LTROC_emp_HPGD_1[p1 == 1] <- 1
UTROC_emp_HPGD_1[p1 == 0] <- 0
UTROC_emp_HPGD_1[p1 == 1] <- 1

t1_HPGD_1 <- -1.1

Sp_emp_HPGD_1 <- mean(HPGD_1_ls[[1]] <= t1_HPGD_1)
all_Se_emp_HPGD_1 <- sapply(HPGD_1_ls[-1], function(y) {
  mean(y > t1_HPGD_1)
})

LSe_emp_HPGD_1 <- min(all_Se_emp_HPGD_1)
id_min_Se_HPGD_1 <- which.min(all_Se_emp_HPGD_1)

x1_Sp_HPGD_1 <- seq(0.01, 0.3, length.out = 300)
y1_Se_HPGD_1 <- seq(0.5, 0.95, length.out = 300)

ll_Sp_LSe_estHPGD_1 <- sapply(x1_Sp_HPGD_1, function(x){
  sapply(y1_Se_HPGD_1, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_HPGD_1, LSe_emp_HPGD_1),
              id_min_max = id_min_Se_HPGD_1)
  })
})

t2_HPGD_1 <- 0.8

Sp_emp_HPGD_1_2 <- mean(HPGD_1_ls[[1]] <= t2_HPGD_1)
all_Se_emp_HPGD_1_2 <- sapply(HPGD_1_ls[-1], function(y){
  mean(y > t2_HPGD_1)
})
USe_emp_HPGD_1 <- max(all_Se_emp_HPGD_1_2)
id_min_Se_HPGD_1_2 <- which.max(all_Se_emp_HPGD_1_2)

x2_Sp_HPGD_1 <- seq(0, 0.3, length.out = 300)
y2_Se_HPGD_1 <- seq(0.6, 1, length.out = 300)

ll_Sp_USe_estHPGD_1 <- sapply(x2_Sp_HPGD_1, function(x){
  sapply(y2_Se_HPGD_1, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_HPGD_1_2, USe_emp_HPGD_1),
              id_min_max = id_min_Se_HPGD_1_2)
  })
})

plot(x = p1, y = LTROC_emp_HPGD_1, type = "l", col = "forestgreen",
     xlim = c(0, 1), ylim = c(0, 1), # xaxs = "i", yaxs = "i",
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid()) #, axes = FALSE)
contour(x1_Sp_HPGD_1, y1_Se_HPGD_1, t(ll_Sp_LSe_estHPGD_1),
        levels = c(4.60517, 5.99, 9.21034), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_HPGD_1, y = LSe_emp_HPGD_1, pch = 16,
       col = "forestgreen")

lines(x = p1, y = UTROC_emp_HPGD_1, col = "blue")
contour(x2_Sp_HPGD_1, y2_Se_HPGD_1, t(ll_Sp_USe_estHPGD_1),
        levels = c(4.60517, 5.99, 9.21034), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_HPGD_1_2, y = USe_emp_HPGD_1, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves ----
p1 <- seq(0, 1, by = 0.001)
LTROC_emp_HPGD_2 <- LTROC_emp(p1, Xlist = HPGD_2_ls)
UTROC_emp_HPGD_2 <- UTROC_emp(p1, Xlist = HPGD_2_ls)

LTROC_emp_HPGD_2[p1 == 0] <- 0
LTROC_emp_HPGD_2[p1 == 1] <- 1
UTROC_emp_HPGD_2[p1 == 0] <- 0
UTROC_emp_HPGD_2[p1 == 1] <- 1

t1_HPGD_2 <- -1.1

Sp_emp_HPGD_2 <- mean(HPGD_2_ls[[1]] <= t1_HPGD_2)
all_Se_emp_HPGD_2 <- sapply(HPGD_2_ls[-1], function(y){
  mean(y > t1_HPGD_2)
})
LSe_emp_HPGD_2 <- min(all_Se_emp_HPGD_2)
id_min_Se_HPGD_2 <- which.min(all_Se_emp_HPGD_2)

x1_Sp_HPGD_2 <- seq(0.01, 0.3, length.out = 300)
y1_Se_HPGD_2 <- seq(0.5, 0.95, length.out = 300)

ll_Sp_LSe_estHPGD_2 <- sapply(x1_Sp_HPGD_2, function(x){
  sapply(y1_Se_HPGD_2, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_HPGD_2, LSe_emp_HPGD_2),
              id_min_max = id_min_Se_HPGD_2)
  })
})

t2_HPGD_2 <- 0.8

Sp_emp_HPGD_2_2 <- mean(HPGD_2_ls[[1]] <= t2_HPGD_2)
all_Se_emp_HPGD_2_2 <- sapply(HPGD_2_ls[-1], function(y){
  mean(y > t2_HPGD_2)
})
USe_emp_HPGD_2 <- max(all_Se_emp_HPGD_2_2)
id_min_Se_HPGD_2_2 <- which.max(all_Se_emp_HPGD_2_2)

x2_Sp_HPGD_2 <- seq(0, 0.3, length.out = 300)
y2_Se_HPGD_2 <- seq(0.6, 1, length.out = 300)

ll_Sp_USe_estHPGD_2 <- sapply(x2_Sp_HPGD_2, function(x){
  sapply(y2_Se_HPGD_2, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_HPGD_2_2, USe_emp_HPGD_2),
              id_min_max = id_min_Se_HPGD_2_2)
  })
})

plot(x = p1, y = LTROC_emp_HPGD_2, type = "l", col = "forestgreen",
     xlim = c(0, 1), ylim = c(0, 1), # xaxs = "i", yaxs = "i",
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid()) #, axes = FALSE)
contour(x1_Sp_HPGD_2, y1_Se_HPGD_2, t(ll_Sp_LSe_estHPGD_2),
        levels = c(4.60517, 5.99, 9.21034), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_HPGD_2, y = LSe_emp_HPGD_2, pch = 16,
       col = "forestgreen")

lines(x = p1, y = UTROC_emp_HPGD_2, col = "blue")
contour(x2_Sp_HPGD_2, y2_Se_HPGD_2, t(ll_Sp_USe_estHPGD_2),
        levels = c(4.60517, 5.99, 9.21034), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_HPGD_2_2, y = USe_emp_HPGD_2, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves: SORBS1 ----
p1 <- seq(0, 1, by = 0.001)
LTROC_emp_SORBS1 <- LTROC_emp(p1, Xlist = SORBS1_ls)
UTROC_emp_SORBS1 <- UTROC_emp(p1, Xlist = SORBS1_ls)

LTROC_emp_SORBS1[p1 == 0] <- 0
LTROC_emp_SORBS1[p1 == 1] <- 1
UTROC_emp_SORBS1[p1 == 0] <- 0
UTROC_emp_SORBS1[p1 == 1] <- 1

t1_SORBS1 <- -1

Sp_emp_SORBS1 <- mean(SORBS1_ls[[1]] <= t1_SORBS1)
all_Se_emp_SORBS1 <- sapply(SORBS1_ls[-1], function(y){
  mean(y > t1_SORBS1)
})
LSe_emp_SORBS1 <- min(all_Se_emp_SORBS1)
id_min_Se_SORBS1 <- which.min(all_Se_emp_SORBS1)

x1_Sp_SORBS1 <- seq(0.01, 0.3, length.out = 300)
y1_Se_SORBS1 <- seq(0.4, 1, length.out = 300)

ll_Sp_LSe_estSORBS1 <- sapply(x1_Sp_SORBS1, function(x){
  sapply(y1_Se_SORBS1, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_SORBS1, LSe_emp_SORBS1),
              id_min_max = id_min_Se_SORBS1)
  })
})

t2_SORBS1 <- 0.5

Sp_emp_SORBS1_2 <- mean(SORBS1_ls[[1]] <= t2_SORBS1)
all_Se_emp_SORBS1_2 <- sapply(SORBS1_ls[-1], function(y){
  mean(y > t2_SORBS1)
})
USe_emp_SORBS1 <- max(all_Se_emp_SORBS1_2)
id_min_Se_SORBS1_2 <- which.max(all_Se_emp_SORBS1_2)

x2_Sp_SORBS1 <- seq(0, 0.3, length.out = 300)
y2_Se_SORBS1 <- seq(0.6, 1, length.out = 300)

ll_Sp_USe_estSORBS1 <- sapply(x2_Sp_SORBS1, function(x){
  sapply(y2_Se_SORBS1, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_SORBS1_2, USe_emp_SORBS1),
              id_min_max = id_min_Se_SORBS1_2)
  })
})

plot(x = p1, y = LTROC_emp_SORBS1, type = "l", col = "forestgreen",
     xlim = c(0, 1), ylim = c(0, 1), # xaxs = "i", yaxs = "i",
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid()) #, axes = FALSE)
contour(x1_Sp_SORBS1, y1_Se_SORBS1, t(ll_Sp_LSe_estSORBS1),
        levels = c(4.60517, 5.99, 9.21034), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_SORBS1, y = LSe_emp_SORBS1, pch = 16,
       col = "forestgreen")

lines(x = p1, y = UTROC_emp_SORBS1, col = "blue")
contour(x2_Sp_SORBS1, y2_Se_SORBS1, t(ll_Sp_USe_estSORBS1),
        levels = c(4.60517, 5.99, 9.21034), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_SORBS1_2, y = USe_emp_SORBS1, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves ----
p1 <- seq(0, 1, by = 0.001)
LTROC_emp_PODXL <- LTROC_emp(p1, Xlist = PODXL_ls)
UTROC_emp_PODXL <- UTROC_emp(p1, Xlist = PODXL_ls)

LTROC_emp_PODXL[p1 == 0] <- 0
LTROC_emp_PODXL[p1 == 1] <- 1
UTROC_emp_PODXL[p1 == 0] <- 0
UTROC_emp_PODXL[p1 == 1] <- 1

t1_PODXL <- -0.4

Sp_emp_PODXL <- mean(PODXL_ls[[1]] <= t1_PODXL)
all_Se_emp_PODXL <- sapply(PODXL_ls[-1], function(y){
  mean(y > t1_PODXL)
})
LSe_emp_PODXL <- min(all_Se_emp_PODXL)
id_min_Se_PODXL <- which.min(all_Se_emp_PODXL)

x1_Sp_PODXL <- seq(0.01, 0.4, length.out = 300)
y1_Se_PODXL <- seq(0.4, 1, length.out = 300)

ll_Sp_LSe_estPODXL <- sapply(x1_Sp_PODXL, function(x){
  sapply(y1_Se_PODXL, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_PODXL, LSe_emp_PODXL),
              id_min_max = id_min_Se_PODXL)
  })
})

t2_PODXL <- 0.2

Sp_emp_PODXL_2 <- mean(PODXL_ls[[1]] <= t2_PODXL)
all_Se_emp_PODXL_2 <- sapply(PODXL_ls[-1], function(y){
  mean(y > t2_PODXL)
})
USe_emp_PODXL <- max(all_Se_emp_PODXL_2)
id_min_Se_PODXL_2 <- which.max(all_Se_emp_PODXL_2)

x2_Sp_PODXL <- seq(0, 0.4, length.out = 300)
y2_Se_PODXL <- seq(0.6, 1, length.out = 300)

ll_Sp_USe_estPODXL <- sapply(x2_Sp_PODXL, function(x){
  sapply(y2_Se_PODXL, function(y){
    ll_emp_2s(n = n, p0 = c(x, y),
              p_est = c(1 - Sp_emp_PODXL_2, USe_emp_PODXL),
              id_min_max = id_min_Se_PODXL_2)
  })
})

plot(x = p1, y = LTROC_emp_PODXL, type = "l", col = "forestgreen",
     xlim = c(0, 1), ylim = c(0, 1), # xaxs = "i", yaxs = "i",
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid()) #, axes = FALSE)
contour(x1_Sp_PODXL, y1_Se_PODXL, t(ll_Sp_LSe_estPODXL),
        levels = c(4.60517, 5.99, 9.21034), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_PODXL, y = LSe_emp_PODXL, pch = 16,
       col = "forestgreen")

lines(x = p1, y = UTROC_emp_PODXL, col = "blue")
contour(x2_Sp_PODXL, y2_Se_PODXL, t(ll_Sp_USe_estPODXL),
        levels = c(4.60517, 5.99, 9.21034), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - Sp_emp_PODXL_2, y = USe_emp_PODXL, pch = 16,
       col = "blue")

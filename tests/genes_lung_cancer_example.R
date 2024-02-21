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

### ---- Empirical estimation for LTAUC/UTAUC ----
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

### ---- LTROC and UTROC curves: HPGD_1 ----
LTROC_emp_HPGD_1 <- LTROC_emp(Xlist = HPGD_1_ls)
UTROC_emp_HPGD_1 <- UTROC_emp(Xlist = HPGD_1_ls)

t1_HPGD_1 <- -1.1

ll_Sp_LSe_est_HPGD_1 <- EL_cr_LTROC(cpt = t1_HPGD_1, Xlist = HPGD_1_ls,
                                    xlim = c(0.01, 0.3), ylim = c(0.5, 0.95))

t2_HPGD_1 <- 0.8

ll_Sp_USe_est_HPGD_1 <- EL_cr_UTROC(cpt = t2_HPGD_1, Xlist = HPGD_1_ls,
                                    xlim = c(0, 0.3), ylim = c(0.6, 1))

plot(x = LTROC_emp_HPGD_1[, 1], y = LTROC_emp_HPGD_1[, 2], type = "l",
     col = "forestgreen", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid())
contour(ll_Sp_LSe_est_HPGD_1$x_Sp, ll_Sp_LSe_est_HPGD_1$y_Se,
        ll_Sp_LSe_est_HPGD_1$ll_Sp_LSe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_LSe_est_HPGD_1$Sp, y = ll_Sp_LSe_est_HPGD_1$LSe,
       pch = 16, col = "forestgreen")

lines(x = UTROC_emp_HPGD_1[, 1], y = UTROC_emp_HPGD_1[, 2], col = "blue")
contour(ll_Sp_USe_est_HPGD_1$x_Sp, ll_Sp_USe_est_HPGD_1$y_Se,
        ll_Sp_USe_est_HPGD_1$ll_Sp_USe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_USe_est_HPGD_1$Sp, y = ll_Sp_USe_est_HPGD_1$USe, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves: HPGD_2 ----
LTROC_emp_HPGD_2 <- LTROC_emp(Xlist = HPGD_2_ls)
UTROC_emp_HPGD_2 <- UTROC_emp(Xlist = HPGD_2_ls)

t1_HPGD_2 <- -1.1
ll_Sp_LSe_est_HPGD_2 <- EL_cr_LTROC(cpt = t1_HPGD_2, Xlist = HPGD_2_ls,
                                    xlim = c(0.01, 0.3), ylim = c(0.5, 0.95))

t2_HPGD_2 <- 0.8
ll_Sp_USe_est_HPGD_2 <- EL_cr_UTROC(cpt = t2_HPGD_2, Xlist = HPGD_2_ls,
                                    xlim = c(0, 0.3), ylim = c(0.6, 1))

plot(x = LTROC_emp_HPGD_2[, 1], y = LTROC_emp_HPGD_2[, 2], type = "l",
     col = "forestgreen", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid())
contour(ll_Sp_LSe_est_HPGD_2$x_Sp, ll_Sp_LSe_est_HPGD_2$y_Se,
        ll_Sp_LSe_est_HPGD_2$ll_Sp_LSe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_LSe_est_HPGD_2$Sp, y = ll_Sp_LSe_est_HPGD_2$LSe, pch = 16,
       col = "forestgreen")

lines(x = UTROC_emp_HPGD_2[, 1], y = UTROC_emp_HPGD_2[, 2], col = "blue")
contour(ll_Sp_USe_est_HPGD_2$x_Sp, ll_Sp_USe_est_HPGD_2$y_Se,
        ll_Sp_USe_est_HPGD_2$ll_Sp_USe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_USe_est_HPGD_2$Sp, y = ll_Sp_USe_est_HPGD_2$USe, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves: SORBS1 ----
LTROC_emp_SORBS1 <- LTROC_emp(Xlist = SORBS1_ls)
UTROC_emp_SORBS1 <- UTROC_emp(Xlist = SORBS1_ls)

t1_SORBS1 <- -1
ll_Sp_LSe_est_SORBS1 <- EL_cr_LTROC(cpt = t1_SORBS1, Xlist = SORBS1_ls,
                                    xlim = c(0.01, 0.3), ylim = c(0.4, 1))

t2_SORBS1 <- 0.5
ll_Sp_USe_est_SORBS1 <- EL_cr_UTROC(cpt = t2_SORBS1, Xlist = SORBS1_ls,
                                    xlim = c(0, 0.3), ylim = c(0.6, 1))

plot(x = LTROC_emp_SORBS1[, 1], y = LTROC_emp_SORBS1[, 2],
     type = "l", col = "forestgreen", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid())
contour(ll_Sp_LSe_est_SORBS1$x_Sp, ll_Sp_LSe_est_SORBS1$y_Se,
        ll_Sp_LSe_est_SORBS1$ll_Sp_LSe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_LSe_est_SORBS1$Sp, y = ll_Sp_LSe_est_SORBS1$LSe, pch = 16,
       col = "forestgreen")

lines(x = UTROC_emp_SORBS1[, 1], y = UTROC_emp_SORBS1[, 2], col = "blue")
contour(ll_Sp_USe_est_SORBS1$x_Sp, ll_Sp_USe_est_SORBS1$y_Se,
        ll_Sp_USe_est_SORBS1$ll_Sp_USe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_USe_est_SORBS1$Sp, y = ll_Sp_USe_est_SORBS1$USe, pch = 16,
       col = "blue")

### ---- LTROC and UTROC curves: PODXL ----
LTROC_emp_PODXL <- LTROC_emp(Xlist = PODXL_ls)
UTROC_emp_PODXL <- UTROC_emp(Xlist = PODXL_ls)

t1_PODXL <- -0.4
ll_Sp_LSe_est_SORBS1 <- EL_cr_LTROC(cpt = t1_PODXL, Xlist = PODXL_ls,
                                    xlim = c(0.01, 0.4), ylim = c(0.4, 1))

t2_PODXL <- 0.2
ll_Sp_USe_est_SORBS1 <- EL_cr_UTROC(cpt = t2_PODXL, Xlist = PODXL_ls,
                                    xlim = c(0, 0.4), ylim = c(0.6, 1))

plot(x = LTROC_emp_PODXL[, 1], y = LTROC_emp_PODXL[, 2], type = "l",
     col = "forestgreen", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "LSe / USe",
     firstpanel = grid())
contour(ll_Sp_LSe_est_SORBS1$x_Sp, ll_Sp_LSe_est_SORBS1$y_Se,
        ll_Sp_LSe_est_SORBS1$ll_Sp_LSe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "forestgreen",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_LSe_est_SORBS1$Sp, y = ll_Sp_LSe_est_SORBS1$LSe, pch = 16,
       col = "forestgreen")

lines(x = UTROC_emp_PODXL[, 1], y = UTROC_emp_PODXL[, 2], col = "blue")
contour(ll_Sp_USe_est_SORBS1$x_Sp, ll_Sp_USe_est_SORBS1$y_Se,
        ll_Sp_USe_est_SORBS1$ll_Sp_USe,
        levels = qchisq(c(0.9, 0.95, 0.99), df = 2), col = "blue",
        labels = c("0.90", "0.95", "0.99"), add = TRUE)
points(x = 1 - ll_Sp_USe_est_SORBS1$Sp, y = ll_Sp_USe_est_SORBS1$USe, pch = 16,
       col = "blue")

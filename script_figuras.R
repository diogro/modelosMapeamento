library(ggplot2)
library(cowplot)
library(plyr)

set.seed(42)
genotypes = c("AA", "AB", "BB")
phenotypes = c(0.1, 0.05, -0.1)
p = 0.5; q = 1 - p
N = 300
pop_index = c(rep(1, floor(N*p^2)), rep(2, floor(N*2*p*q)), rep(3, floor(N*q^2)))
pop_index_2 = sample(pop_index)
pop = phenotypes[pop_index] + 3 + rnorm(N, 0, 0.2)
pop_df = data.frame(Genotypes = factor(genotypes[pop_index], levels = c("BB", "AB", "AA")), 
                    Phenotype = pop)
medias = ddply(pop_df, .(Genotypes), numcolwise(mean))


set.seed(42)
anova_plot = ggplot(pop_df, aes(Genotypes, Phenotype)) + geom_jitter(size = 2, width = 0.2) + 
  geom_point(data = medias, size = 5, color = "red") + 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
save_plot("~/anova_plot.png", anova_plot, base_height = 4, base_aspect_ratio = 1.5, bg = "transparent")

set.seed(42)
centered_anova_plot = ggplot(pop_df, aes(Genotypes, Phenotype)) + geom_jitter(size = 2, width = 0.2) + geom_hline(yintercept = 3) +
  scale_y_continuous(breaks = c(1, 0.5, 0, -1)+3, labels = c("a", "d", "0", "-a"))+ 
  geom_point(data = medias, size = 5, color = "red")+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
save_plot("~/centered_anova_plot.png", centered_anova_plot, base_height = 4, base_aspect_ratio = 1.5, bg = "transparent")

set.seed(42)
additive_regression = ggplot(pop_df, aes(Genotypes, Phenotype)) + geom_jitter(size = 2, width = 0.2) +
  scale_y_continuous(breaks = c(0.1, 0.0, 0.05, -0.1)+3, labels = c("a", "d", "0", "-a"))+ 
  scale_x_discrete(labels = c(-1, 0, 1)) + labs (x = "Genotypes", y = "Phenotype") +
  geom_point(data = medias, size = 5, color = "red") + ggtitle("Additive effects") + 
  geom_abline(slope = 0.1, intercept = 2.8)+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
save_plot("~/non_sig_additive_regression.png", additive_regression, base_height = 4, base_aspect_ratio = 1.5, bg = "transparent")


pop_df$additive = ifelse(pop_df$Genotypes == "LS", 0, ifelse(pop_df$Genotypes=="LL", 1, -1)) 
pop_df$dominance = ifelse(pop_df$Genotypes == "LS",1,0)
medias_dm = ddply(pop_df, .(dominance), numcolwise(mean))
dominance_regression = ggplot(pop_df, aes(dominance, Phenotype)) + 
  geom_jitter(size = 2, width = 0.1) +
  scale_y_continuous(breaks = c(1, 0.5, 0, -1)+3, labels = c("a", "d", "0", "-a"))+ 
  scale_x_continuous(breaks = c(0, 1))+ 
  labs (x = "Genotypes", y = "Phenotype") +
  geom_point(data = medias_dm, size = 5, color = "red") + 
  ggtitle("Regressão dos efeitos de dominância") + 
  geom_abline(slope = 0.5, intercept = 3)+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

ortogonal_regression = plot_grid(additive_regression, dominance_regression)
save_plot("apresentacao/images/ortogonal_regression_plot.png", ortogonal_regression, base_height = 4, base_aspect_ratio = 1, ncol = 2, bg = "transparent")

summary(lm(Phenotype ~ additive + dominance, data = pop_df))

list_pkgs <- c("plyr", "dplyr", "tidyr", "readr", "car", "ggplot2", "cowplot")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(car)
library(ggplot2)
library(cowplot)

f2_data = read_csv("~/projects/QGcourse/tutorials/F2 geno pheno with QTL effect.csv")
n_markers = 31
marker_fits_Trait2a = list()
for(i in seq(n_markers)){
  current_marker = paste0('G', i)
  f2_data[[paste0(current_marker, '_D')]] = ifelse(f2_data[[current_marker]], 0, 1)
  model_formula = paste0("Trait2a ~ Sex + LSB + LSW + ", current_marker, "+", current_marker, "_D")
  marker_fits_Trait2a[[i]] = lm(as.formula(model_formula), data = f2_data)
}

LPR_plot = ldply(marker_fits_Trait2a, function(x) summary(x)$coefficients[5, 'Pr(>|t|)']) %>%
  ggplot(aes(1:n_markers, -log10(V1))) + geom_line() + labs(x = "marker", y = "LPR") +
  scale_x_continuous(breaks = seq(1, n_markers, by = 2))+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
save_plot("apresentacao/images/LPR_f2.png", LPR_plot, base_height = 4, base_aspect_ratio = 1.5, bg = "transparent")

n_markers = 31
marker_fits_Trait3a = list()
for(i in seq(n_markers)){
  current_marker = paste0('G', i)
  model_formula = paste0("Trait3a ~ Sex + LSB + LSW + ", current_marker, "+", current_marker, "_D")
  marker_fits_Trait3a[[i]] = lm(as.formula(model_formula), data = f2_data)
}

LPR_2 = ldply(marker_fits_Trait3a, function(x) summary(x)$coefficients[c(5, 6), 'Pr(>|t|)']) %>% 
  ggplot(aes(1:n_markers, -log10(G1))) + geom_line() + labs(x = "marker", y = "LPR") +
  scale_x_continuous(breaks = seq(1, n_markers, by = 2))+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
multipleQTL = plot_grid(LPR_plot + ggtitle("Um QTL"), LPR_2 + ggtitle("Dois QTLs"))
save_plot("apresentacao/images/multiQTL_reg.png", multipleQTL, base_height = 4, base_aspect_ratio = 1, ncol = 2, bg = "transparent")


n_markers = 31
interval = 6 # Play around with this value. 
fl_marker_fits_Trait3a = list()
for(i in seq(n_markers)){
  current_marker = paste0('G', i)
  f2_data[[paste0(current_marker, '_D')]] = ifelse(f2_data[[current_marker]], 0, 1)
  model_formula = paste0("Trait3a ~ Sex + LSB + LSW + ", current_marker, "+", current_marker, "_D")
  if((i - interval) >= 1)
    model_formula = paste0(model_formula, "+ G", paste0(i - interval), "+ G", paste0(i - interval), "_D")
  if((i + interval) <= n_markers)
    model_formula = paste0(model_formula, "+ G", paste0(i + interval), "+ G", paste0(i + interval), "_D")
  fl_marker_fits_Trait3a[[i]] = lm(as.formula(model_formula), data = f2_data)
}

n_markers = 31
interval = 6 # Play around with this value. 
fl_marker_fits_Trait2a = list()
for(i in seq(n_markers)){
  current_marker = paste0('G', i)
  f2_data[[paste0(current_marker, '_D')]] = ifelse(f2_data[[current_marker]], 0, 1)
  model_formula = paste0("Trait2a ~ Sex + LSB + LSW + ", current_marker, "+", current_marker, "_D")
  if((i - interval) >= 1)
    model_formula = paste0(model_formula, "+ G", paste0(i - interval), "+ G", paste0(i - interval), "_D")
  if((i + interval) <= n_markers)
    model_formula = paste0(model_formula, "+ G", paste0(i + interval), "+ G", paste0(i + interval), "_D")
  fl_marker_fits_Trait2a[[i]] = lm(as.formula(model_formula), data = f2_data)
}

LPR_1_int = ldply(fl_marker_fits_Trait2a, function(x) summary(x)$coefficients[c(5, 6), 'Pr(>|t|)'])%>% 
  ggplot(aes(1:n_markers, -log10(G1))) + geom_line() + labs(x = "marker", y = "LPR") +
  scale_x_continuous(breaks = seq(1, n_markers, by = 2))+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

LPR_2_int = ldply(fl_marker_fits_Trait3a, function(x) summary(x)$coefficients[c(5, 6), 'Pr(>|t|)'])%>% 
  ggplot(aes(1:n_markers, -log10(G1))) + geom_line() + labs(x = "marker", y = "LPR") +
  scale_x_continuous(breaks = seq(1, n_markers, by = 2))+ 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
multipleQTLint = plot_grid(LPR_1_int + ggtitle("Um QTL"), LPR_2_int + ggtitle("Dois QTLs"))
save_plot("apresentacao/images/multiQTL_int.png", multipleQTLint, base_height = 4, base_aspect_ratio = 1, ncol = 2, bg = "transparent")


ggplot(data.frame(x = c(-0, 1)), aes(x)) +
  stat_function(fun = function(x) dbeta(x, 0.5, 0.5), geom = "line")

par(bg=NA) 
pdf("apresentacao/images/beta_dist.pdf")
curve(dbeta(x, 0.5, 0.5), bty="l", xlab = expression('w'[j]), ylab = "Density")
dev.off()

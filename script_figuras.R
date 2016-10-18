library(ggplot2)
library(cowplot)
library(plyr)

set.seed(42)
genotypes = c("LL", "LS", "SS")
phenotypes = c(1, 0.5, -1)
p = 0.5; q = 1 - p
N = 300
pop_index = c(rep(1, floor(N*p^2)), rep(2, floor(N*2*p*q)), rep(3, floor(N*q^2)))
pop = phenotypes[pop_index] + 3 + rnorm(N, 0, 0.2)
medias = ddply(pop_df, .(genótipos), numcolwise(mean))
pop_df = data.frame(genótipos = factor(genotypes[pop_index], levels = c("SS", "LS", "LL")), 
                    Fenótipo = pop)

set.seed(42)
anova_plot = ggplot(pop_df, aes(genótipos, Fenótipo)) + geom_jitter(size = 2, width = 0.2) + 
  geom_point(data = medias, size = 5, color = "red")
print(anova_plot)

set.seed(42)
centered_anova_plot = ggplot(pop_df, aes(genótipos, Fenótipo)) + geom_jitter(size = 2, width = 0.2) + geom_hline(yintercept = 3) +
  scale_y_continuous(breaks = c(1, 0.5, 0, -1)+3, labels = c("a", "d", "0", "-a"))+ 
  geom_point(data = medias, size = 5, color = "red")
print(centered_anova_plot)

set.seed(42)
additive_regression = ggplot(pop_df, aes(genótipos, Fenótipo)) + geom_jitter(size = 2, width = 0.2) +
  scale_y_continuous(breaks = c(1, 0.5, 0, -1)+3, labels = c("a", "d", "0", "-a"))+ 
  scale_x_discrete(labels = c(-1, 0, 1)) + labs (x = "Genótipos", y = "Fenótipo") +
  geom_point(data = medias, size = 5, color = "red") + ggtitle("Regressão dos efeitos aditivos") + 
  geom_abline(slope = 1, intercept = 1)
print(additive_regression)

pop_df$dominance = ifelse(pop_df$genótipos == "LS",1,0)
dominance_regression = ggplot(pop_df, aes(ifelse(genótipos == "LS",1,0) , Fenótipo)) + 
  geom_jitter(size = 2, width = 0.1) +
  scale_y_continuous(breaks = c(1, 0.5, 0, -1)+3, labels = c("a", "d", "0", "-a"))+ 
  scale_x_continuous(breaks = c(0, 1))+ 
  labs (x = "Genótipos", y = "Fenótipo") +
  geom_point(data = medias, size = 5, color = "red") + 
  ggtitle("Regressão dos efeitos de dominância") + 
  geom_abline(slope = 0.5, intercept = 3)
print(dominance_regression)
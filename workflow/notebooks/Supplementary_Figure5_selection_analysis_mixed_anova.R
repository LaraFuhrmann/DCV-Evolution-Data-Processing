# A mixed ANOVA was performed where the viral passage is treated as
# within-subject factor and protein and pi-type (piN or piS) are treated as between-subject
# factors. To ensure approximate normal distribution of the residuals, the data the arcsin(sqrt(x))-transformation was applied. 
# Equal variance was ensured using Levene's test. 

base_path = "/Users/lfuhrmann/Documents/Projects/DCV-Evolution/DCV-Evolution-Data-Processing/workflow/notebooks"
df <- read.csv(paste(base_path,"/df_pis_protein_long.csv", sep = ""))

df$pi_cuberoot <- (df$pi)^(1/3)
df$id <- paste(df$genotype, df$replicate, df$protein, df$type) # subject id

library(ggpubr)
library(rstatix)

# check normality
# nice resource on qq-plots: https://onlinestatbook.com/2/advanced_graphs/q-q_plots.html 
ggqqplot(df[(df$protein=="1A") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="1A") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="2A") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="2A") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="2B") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="2B") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="2C") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="2C") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="3A") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="3A") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="3C") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="3C") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="RdRp") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="RdRp") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="VP2") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="VP2") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="VP4") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="VP4") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="VP3") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="VP3") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

ggqqplot(df[(df$protein=="VP1") & (df$type=="S"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")
ggqqplot(df[(df$protein=="VP1") & (df$type=="N"),], "pi_cuberoot", ggtheme = theme_bw()) + facet_grid(passage ~genotype, labeller = "label_both")

# check equal variance 
df %>%
  group_by(passage) %>%
  levene_test(pi_cuberoot ~ factor(protein)*factor(type)*factor(genotype))

## mixed ANOVA
res.aov <- anova_test(
  data = df, 
  dv = pi_cuberoot, 
  wid = id,
  within = passage, 
  between = c(genotype, protein, type)
)
get_anova_table(res.aov,  correction = "GG")
## result
##                 Effect   DFn   DFd     F        p p<.05   ges
##2                        protein 10.00 264.00 25.405 1.84e-33     * 0.272000
##3                           type  1.00 264.00 14.007 2.24e-04     * 0.020000
##4                        passage  1.44 378.85 10.192 3.35e-04     * 0.023000
##7                   protein:type 10.00 264.00 18.247 3.01e-25     * 0.212000
##10                  protein:passage 14.35 378.85  1.748 4.30e-02     * 0.039000
##14          type:passage  1.44 378.85  9.034 7.69e-04     * 0.021000
##protein:type:passage 14.35 378.85  2.554 1.00e-03     * 0.056000


pwc <- df %>%
  group_by(passage) %>%
  tukey_hsd(pi_cuberoot ~ protein:type, p.adjust.method = "bonferroni")


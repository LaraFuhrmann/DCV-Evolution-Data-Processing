# Figure 1C -- viral titers
# use mixed ANOVA
# within-subject factor: passage
# indepentent between factor: genotype 
# dependent variable: dcv titers

base_path = "/Users/lfuhrmann/Documents/Projects/DCV-Evolution/DCV-Evolution-Data-Processing/resources"
df <- read.csv(paste(base_path,"/measurements/20220926-RNAi_evolution_Titers_overview_formatted.csv", sep = ""))

df$DCV_titer_log <- 1/(df$DCV_titer)
df$id <- paste(df$genotype, df$passage, df$replicate)

df %>% 
  group_by(genotype, passage) %>% 
  get_summary_stats(DCV_titer)

library(rstatix)
library(ggpubr)

df %>% 
  anova_test(dv = DCV_titer_log, wid = id, between = genotype, within = passage, type = 3)


#test anova assumptions
# check for outliers
df %>%
  group_by(passage, genotype) %>%
  identify_outliers(DCV_titer_log)
# --> since we have outliers do robust anova
library(WRS2)
bwtrim(formula = DCV_titer ~ genotype * passage, id = id, data = df)
## between groups effect only (MOM estimator)
sppbi(DCV_titer_log ~ genotype*passage, id, data = df)



ggqqplot(df, "DCV_titer_log", ggtheme = theme_bw()) +
  facet_grid(passage ~ genotype)

df %>%
  group_by(passage) %>%
  levene_test(DCV_titer_log ~ genotype)

box_m(df[, "DCV_titer_log", drop = FALSE], df$genotype)

df$id <- paste(df$genotype, df$passage)

# Two-way mixed ANOVA test
res.aov <- anova_test(
  data = df, 
  dv = DCV_titer_log,
  wid = id,
  between = genotype, 
  within = passage
)
get_anova_table(res.aov)

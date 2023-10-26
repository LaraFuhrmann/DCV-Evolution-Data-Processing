base_path = "/Users/lfuhrmann/Documents/Projects/DCV-Evolution/DCV-Evolution-Data-Processing/workflow/notebooks"
df <- read.csv(paste(base_path,"/df_pis_protein_long.csv", sep = ""))

df$pi_log <- log10(df$pi+1)

#fit ANOVA
model_formula <- pi ~ factor(type) * (factor(passage) * factor(genotype) * factor(protein) * factor(type) )
fitted_anova <- aov(model_formula, data=df)

### check ANOVA assumptions

# check if residuals follow approximalty normal distribution 
hist(residuals(fitted_anova))
qqnorm(rstandard(fitted_anova))
qqline(rstandard(fitted_anova))

#conduct Levene's Test for equality of variances
library(car)
leveneTest(model_formula, data=df)


#view summary of nested ANOVA
summary(fitted_anova)

library(rstatix)
tuk<-tukey_hsd(fitted_anova, which="factor(type):factor(protein)") 


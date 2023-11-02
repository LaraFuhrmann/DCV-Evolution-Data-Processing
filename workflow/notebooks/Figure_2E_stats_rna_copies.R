# Using a nested ANOVA to account for the technical replicates, we are comparing 
# the RNA copies of the viral populations parental stock to the passage 10 evolved 
# virues in the three fly genotypes. 

# Load data
base_path = "/Users/lfuhrmann/Documents/Projects/DCV-Evolution/DCV-Evolution-Data-Processing/workflow/notebooks"
df <- read.csv(paste(base_path,"/RNA_copies_for_R_script.csv", sep = ""))

#exlcude passage 1
df <- df[!df$passage==1,]

df$RNA_copies_log <- log(df$RNA_copies)

df$replicate <- paste(df$genotype, df$biological_replicate)
#df$id <- paste(df$genotype,  df$biological_replicate)

#fit nested ANOVA
#nested_repeated_formula <- RNA_copies_log ~ factor(genotype) * factor(passage)*factor(replicate)  / factor(Technical.replicate) + Error(factor(id)/factor(passage))
nested_formula <- RNA_copies_log ~ factor(genotype) * factor(replicate)  / factor(Technical.replicate)
nest <- aov(data=df,nested_formula)

#view summary of nested ANOVA
summary(nest)

library(rstatix)
tuk<-tukey_hsd(nest, which="factor(genotype)")
tuk

### check ANOVA assumptions
# check if residuals follow approximalty normal distribution 
hist(residuals(nest))
qqnorm(rstandard(nest))
qqline(rstandard(nest))

#conduct Levene's Test for equality of variances
library(car)
leveneTest(nested_formula, data=df)

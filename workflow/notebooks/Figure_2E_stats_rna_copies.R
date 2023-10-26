
# Load data
base_path = "/Users/lfuhrmann/Documents/Projects/DCV-Evolution/DCV-Evolution-Data-Processing/workflow/notebooks"
df <- read.csv(paste(base_path,"/RNA_copies_for_R_script.csv", sep = ""))

# exlcude passage 1
#df <- df[!df$passage==1,]

df$RNA_copies_log <- log(df$RNA_copies)

#fit nested ANOVA
model_formula <- RNA_copies_log ~factor(genotype) * factor(passage) * factor(biological_replicate)  / factor(Technical.replicate)
nest <- aov(data=df,model_formula)

#view summary of nested ANOVA
summary(nest)

tuk<-tukey_hsd(nest, which="factor(genotype):factor(passage)")

### check ANOVA assumptions
# check if residuals follow approximalty normal distribution 
hist(residuals(nest))
qqnorm(rstandard(nest))
qqline(rstandard(nest))

#conduct Levene's Test for equality of variances
library(car)
leveneTest(model_formula, data=df)

####### R codes for Scientists data ####### 
# Last updated by Minhui Li: Mar 18, 2024

# set working directory
setwd("") 

# load data
df1<-read.csv("scientist_data_20240318.csv",header=TRUE,sep=',')

# calculate index
Social_Influence_Index <- rowMeans(df1[,21:24])
df2 <- cbind(df1, Social_Influence_Index)
str(df2)


########## Correlation analyses ##########
# Kruskal-Wallis tests for categorical variables
dis_variables <- df2[,c(2,18,19,20)]
kruskal.test(Social_Influence_Index ~ pri_con_area, data = dis_variables)
kruskal.test(Social_Influence_Index ~ pri_collaborator, data = dis_variables)
kruskal.test(Social_Influence_Index ~ pri_com_channel, data = dis_variables)
t.test(Social_Influence_Index ~ sex, data = dis_variables)


# test of profile variables
shapiro.test(x=df2$age) #age
t.test(Social_Influence_Index ~ sex, data = dis_variables) # sex
shapiro.test(x=df2$years) # # of years of experience
shapiro.test(x=df2$staff) # # of staff
kruskal.test(Social_Influence_Index ~ staff, data = df2) # staff
shapiro.test(x=df2$X.female) # female
kruskal.test(Social_Influence_Index ~ prop, data = df2) # dedication
kruskal.test(Social_Influence_Index ~ pri_con_area, data = dis_variables)  # Primary concerned area


# Correlation analyses for continuous/ordered variables
library(ggcorrplot)
library(caret)
con_variables <- df2[,-c(2,18,19,20,21,22,23,24,25)]
corr2 <- round(cor(con_variables, method = "spearman"), 2)
findCorrelation(corr2)
p.mat2 <- cor_pmat(con_variables)
ggcorrplot(corr2, p.mat = p.mat2, hc.order = F, type = "lower", lab = TRUE,insig = "blank") # plot initial correlation heat map

corr_variables <- con_variables[, c(3,4,7,8,10,11,13,15,16,17)] # select variables that were significantly correlated with the SII 
colnames(corr_variables) <- c("# of years of experience", "# of staff", "Annual fund", 
                              "# of projects", "Frequency of collaboration", "Conversion rate", 
                              "Expert consultation", "Conversion approaches", "Communication capability",
                              "Potential Social Influence Index")
corr_v <- round(cor(corr_variables, method = "spearman"), 2)
p.mat_v <- cor_pmat(corr_variables)
ggcorrplot(corr_v, p.mat = p.mat_v, hc.order = F, type = "lower", lab = TRUE,insig = "blank", 
           colors = c("#6D9EC1", "white", "#E46726"),
           tl.cex = 18,
           legend.title = "Correlation coefficient")  # plot final correlation heat map

# Save the data
write.csv(corr_v,"Correlation_matrix_sci.csv", row.names=T, quote=F)
write.csv(p.mat_v,"Corr_pvalue_matrix_sci.csv", row.names=T, quote=F)
write.csv(corr_variables,"Corr_variables_data_sci.csv", row.names=T, quote=F)


########## Validation ##########
library(ggplot2)
library(ggpmisc)
library(ggpubr)

p = ggplot(aes(sel_eva, Social_Influence_Index), data=df2)+
  geom_point(size=4,alpha=0.3,color="#6baed6")+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+ 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )+
  xlab("Self-evaluation scores of social influence")+
  ylab("PSII of scientists")
p + stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=5)

ggsave("sci_validation.pdf",width = 16,height = 12,units = "cm")


########## Generalized linear models ##########
# Choose all related variable to build glm model
colnames(con_variables)
glm1<-glm(Social_Influence_Index ~ years
          + staff
          + fund
          + X.projects
          + co_frq
          + X.conv
          + exp_con
          + X.con_app
          + com_cap, 
          data = con_variables, family = gaussian) 
summary(glm1)

# Derive optimal GLM based on backward stepwise selection
glm2<-step(glm1, direction='backward', scope=formula(glm1), trace=0)
glm2
summary(glm2)
glm2$anova # see the variables that were filtered out of glm2

# Scale variables using z-score
library("broom.mixed")
library(jtools)
z_sci_estimate <- summ(glm2, scale=T, robust="HC1") # scale variables using z-score
z_sci_estimate

# Plot the optimal GLM coefficients
plot_summs(z_sci_estimate, plot.distributions = TRUE, inner_ci_level = .9, scale=T)
z_sci_estimate <- tidy(z_sci_estimate)

# Save the data
write.csv(z_sci_estimate, "z_sci_estimate.csv", row.names=T, quote=F)


########## Pathway analyses ##########
# Conduct path analyses to examine the potential casual structure
my_packages <- c("ape", "caper", "nlme", "lavaan","tidySEM","processx","piecewiseSEM") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]  # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                
# Load required libraries
library(ape) #Version 3.3
library(caper) # Vresion 0.5.2
library(nlme) # Version 3.1.122
library(lavaan) # Version 0.5.19
library(tidySEM)
library(piecewiseSEM) # Version 1.0.0


# Read in data
sci_data<-read.csv("scientist_data_20240318.csv",header=TRUE,sep=',') # load data
Social_Influence_Index <- rowMeans(sci_data[,21:24])
scientist <- cbind(sci_data, Social_Influence_Index)
str(scientist)
con_variables <- scientist[,-c(18,19,20,21,22,23,24)]
corr_variables <- con_variables[, c(4,5,8,9,11,12,14,16,17,18)]
scientist = na.omit(corr_variables)
colnames(corr_variables)

# Now fit piecewise model with random effect
psem1<-psem(lm(fund ~ years + staff + co_frq, data = scientist),
            lm(exp_con ~ co_frq , data = scientist),
            lm(com_cap ~ co_frq + years + fund, data = scientist), 
            lm(X.projects ~ years, data = scientist),
            lm(Social_Influence_Index ~ years + exp_con + com_cap + fund + X.conv + X.projects, data = scientist))

# Evaluate model fitting
summary(psem1)
sum_psem1 <- summary(psem1)

# Save the data
capture.output(sum_psem1, file = "sci_psem_summary.txt")

# Plot the model
plot(psem1,
     node_attrs = data.frame(shape = "rectangle", color = "blue", fillcolor = "white"),
     edge_attrs = data.frame(style = "solid", color = "red"),
     ns_dashed = T,
     alpha = 0.05,
     show = "std",
     digits = 3,
)



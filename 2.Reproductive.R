#estimation of length at maturity of fish
#data required (lenth in mm, sex, maturity stages in binary form (0=immature & 1= mature))


# Specify libraries required for analysis
library(MASS)
#install.packages("psyphy")
library(psyphy)
library(boot)
#install.packages("dplyr")
library(dplyr)
library(ggplot2)

#data <- read.csv('length_maturity.csv')
#head(data)
maturity_data<-kawakawa%>%select(FL_cm, Sex..m.f.,Binary_maturity)
colnames(maturity_data)<-c("FL","Sex","Maturity.Stage")
head(maturity_data)

########Calculate sex ratio (M:F)############
sex_count <- table(maturity_data$Sex)

F_no <- 899
M_no <- 534

sex_ratio <- M_no / F_no

cat("Males =", M_no, "\n")
cat("Females =", F_no, "\n")
cat("Sex ratio (M:F) =", round(sex_ratio, 2), ":1\n")

#Chi-square test for sex ratio significant
observed_sex <- c(M_no, F_no)

chisq.test(
  x = observed_sex,
  p = c(0.5, 0.5)
)

######### Estimate Lm50 ###########
# Specify sex: "m", "f" or "both"

sex="both"

# original data (data1)
maturity_data1<-maturity_data
maturity_data<-if(sex=="both"){maturity_data}else{maturity_data[maturity_data$Sex==sex,]}
maturity_data<-na.omit(maturity_data)
attach(maturity_data)

# Fit a generalised linear model with a logit link function. STL is 
# stretched total length in mm. Maturity stage is the binomial response 
# variable; 1 for mature and 0 for immature
Maturity.Model<-glm(Maturity.Stage~FL,family=binomial(logit))

# A quick way for obtaining the length at which a specified proportion of
# the population is mature is using the dose.p function, where p = proportion
# mature
L50<- dose.p(Maturity.Model,p=c(0.5))
L95<- dose.p(Maturity.Model,p=c(0.95))

# Alternatively, the fitted parameters of the model can be re-arranged
# so that the model can be parameterised in terms of the key parameters
# of interest (the length at which 50% and 95% of population is mature)

L50<-(Maturity.Model$coef[1]*Maturity.Model$coef[2]^-1)
L95<-(1/Maturity.Model$coef[2])*log(1/0.05-1)-
  Maturity.Model$coef[1]/Maturity.Model$coef[2]


# Tablulate the raw maturity data into some length inverval (Lint), say 2 cm
Lint=2
min.x<-min(maturity_data$FL[maturity_data$Sex==Sex],na.rm=T)
max.x<-max(maturity_data$FL[maturity_data$Sex==Sex],na.rm=T)
max.x.plot<-ceiling(max.x/Lint)*Lint

Maturity.Table<-table(Maturity.Stage,cut(FL,seq(0,max.x.plot,Lint)))
# Save tabulated data
#write.csv(t(Maturity.Table),"Maturity.Table_both.csv")

Observed.Proportions<-Maturity.Table[2,]/(Maturity.Table[1,]+Maturity.Table[2,])

# Everything below here is just bootstrap-resampling and plotting

# Set up a vector of x-values over which to predict length at maturity
#if want to keep min proportion of maturity in y axis replace '0.0' to 'min.x'
new <- data.frame(FL = seq(0.0,max.x,0.01))

plot(new$FL,predict(Maturity.Model, data.frame(FL=new$FL),type="response"),
     type="l",las=1,ylab="Proportion mature",xlab="Fork Length (Cm)",bty="l")

# Add observed proportions calculated above if desired
points(seq(Lint/2,max.x.plot,Lint),Observed.Proportions,pch=19)

# Function for bootstrapping the data. This simply repeats the above analysis. 
# The parameters are also expressed in terms of L50 and L95 
boot.kawakawa <- function(data, i){
  
  data <- data[i, ]
  
  mod <- glm(
    Maturity.Stage ~ FL,
    family = binomial(link = "logit"),
    data = data
  )
  
  L50 <- dose.p(mod, p = 0.5)
  L95 <- dose.p(mod, p = 0.95)
  
  c(coef(mod), L50, L95)
}


# Run bootstrap

R=1000

Maturity.Bootstraps <- boot(maturity_data,
  statistic = boot.kawakawa,
  R = R
)
# Function for plotting confidence intervals around curves using
# bootstrap outputs
CI.Plotting <- Confints(
  maturity_data,
  Maturity.Bootstraps,
  Maturity.Model,
  R
)

######## Lm Plot #############
par(mfrow = c(1,1))

#Lm values from model
Lm50 <- 35.1
Lm95 <- 61.7
# recreating CI and Obs data frame
CI.df <- data.frame(
  FL = CI.Plotting$Length,
  Fit = CI.Plotting$fitted.maturity,
  Lower = CI.Plotting$boundsmatrix[1,],
  Upper = CI.Plotting$boundsmatrix[2,]
)

Obs.df <- data.frame(
  FL = seq(Lint/2, max.x.plot, Lint),
  Prop = Observed.Proportions
)

#check both data raw length
length(seq(Lint/2, max.x.plot, Lint))
length(Observed.Proportions)

head(Obs.df)

#####plot
p3<-ggplot() +
  
  # 95% bootstrap confidence band
  geom_ribbon(
    data = CI.df,
    aes(x = FL, ymin = Lower, ymax = Upper),
    fill = "red",
    alpha = 0.20
  ) +
  
  # maturity ogive
  geom_line(
    data = CI.df,
    aes(x = FL, y = Fit),
    colour = "red",
    linewidth = 1.3
  ) +
  
  # observed proportions
  geom_point(
    data = Obs.df,
    aes(x = FL, y = Prop),
    colour = "darkgreen",
    size = 2.5
  ) +
  
  # Lm50 and Lm95
  geom_vline(
    xintercept = Lm50,
    colour = "blue",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = Lm95,
    colour = "brown",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # 50% maturity line
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +
  
  # Lm labels top-left
  annotate(
    "text",
    x = 0.5,
    y = 0.98,
    hjust = 0,
    label = paste0(
      "Lm50 = ", round(Lm50,1), " cm\n",
      "Lm95 = ", round(Lm95,1), " cm"
    ),
    size = 5
  ) +
  
  labs(
    x = "Fork Length (cm)",
    y = "Proportion mature"
  ) +
  
  coord_cartesian(ylim = c(0,1)) +
  
  theme_classic(base_size = 14) +
  
  theme(
    #axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )
##save plot
ggsave(
  "Kawakawa_Maturity_Ogive.jpg",
  plot = p3,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)



#####plot withot observation point
p4<-ggplot() +
  
  # 95% bootstrap confidence band
  geom_ribbon(
    data = CI.df,
    aes(x = FL, ymin = Lower, ymax = Upper),
    fill = "red",
    alpha = 0.20
  ) +
  
  # maturity ogive
  geom_line(
    data = CI.df,
    aes(x = FL, y = Fit),
    colour = "red",
    linewidth = 1.3
  ) +
  
  # observed proportions
  #geom_point(
    #data = Obs.df,
    #aes(x = FL, y = Prop),
    #colour = "darkgreen",
    #size = 2.5
  #) +
  
  # Lm50 and Lm95
  geom_vline(
    xintercept = Lm50,
    colour = "blue",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = Lm95,
    colour = "brown",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # 50% maturity line
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +
  
  # Lm labels top-left
  annotate(
    "text",
    x = 0.5,
    y = 0.98,
    hjust = 0,
    label = paste0(
      "Lm50 = ", round(Lm50,1), " cm\n",
      "Lm95 = ", round(Lm95,1), " cm"
    ),
    size = 5
  ) +
  
  labs(
    x = "Fork Length (cm)",
    y = "Proportion mature"
  ) +
  
  coord_cartesian(ylim = c(0,1)) +
  
  theme_classic(base_size = 14) +
  
  theme(
    #axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )
##save plot
ggsave(
  "Kawakawa_Maturity_Ogive1.jpg",
  plot = p4,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#############fecundity#######

#import data#
fecundity_data<- read.csv("fecun_data.csv")
head(fecundity_data)
colnames(fecundity_data)<- c("TL_cm","FL_cm","Fecundity")

# Identify rows with missing FL
missing_fl1 <- is.na(fecundity_data$FL_cm)

# Predict FL from TL with round zero decimal
fecundity_data$FL_cm[missing_fl1] <- round(
  predict(
    model1,
    newdata = fecundity_data[missing_fl1, ]
  ),
  0
)

# View updated data
View(fecundity_data)

#remove NA raw
fecundity_data<-na.omit(fecundity_data)

#Min, Max and avg fecundity
fecundity.min<- min(fecundity_data$Fecundity)
fecundity.max<- max(fecundity_data$Fecundity)
fecundity.avg<- round(mean(fecundity_data$Fecundity),0)

####Est. relationship b/w FL and fecundity
head(fecundity_data)

fecun_model <- lm(Fecundity ~ FL_cm, data = fecundity_data)
summary(fecun_model)

# Extract coefficients
a1 <- coef(fecun_model)[1]
a1<-exp(a1)
b1 <- round(coef(fecun_model)[2],2)

#########plot preparation#######
# Extract model statistics

r2 <- round(summary(fecun_model)$r.squared, 4)
n <- nrow(fecun_model)

eq <- paste0(
  "Fecundity = ", a1, " + ", b1, " FL_cm",
  "\nR² = ", r2,
  "\nn = ", n
)
  

#ggplot of log-transformed model
#library(ggplot2)

p8 <- ggplot(fecun_model, aes(x = FL_cm, y = Fecundity)) +
  
  geom_point(
    colour = "black",
    size = 2,
    alpha = 0.7
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "blue",
    fill = "pink",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = min(fecundity_data$FL_cm),
    y = max(fecundity_data$Fecundity),
    label = eq,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  
  labs(
    x = "Fork Length cm",
    y = "Fecundity",
    title = ""
  ) +
  
  theme_bw(base_size = 14)

p8

#save plot
ggsave(
  filename = "Kawakawa_L-fecundity.jpg",
  plot = p8,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

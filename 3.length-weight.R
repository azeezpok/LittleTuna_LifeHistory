# Source data check
head(kawakawa)

# prepare dataset
LWdata <-kawakawa%>%select(FL_cm, Total.weight_g)
colnames(LWdata)<-c("length","weight")
min(LWdata1$weight); max(LWdata1$weight)
#weight column contains blank values ("") rather than NA, 
#first identify them and then remove those rows.

#Check the structure
str(LWdata)

#If weight is a character column, use:
  
# Remove rows where weight is blank
LWdata <- LWdata[LWdata$weight != "", ]

# Convert weight to numeric
LWdata$weight <- as.numeric(LWdata$weight)


LWdata1<-na.omit(LWdata)
head(LWdata1)
#View(LWdata1)

# Log-transform
LWdata1$logL <- log(LWdata1$length)
LWdata1$logW <- log(LWdata1$weight)

# Fit log-log regression: logW = a + b*logL
LWmodel <- lm(logW ~ logL, data = LWdata1)
summary(LWmodel)

# Extract coefficients
a <- coef(LWmodel)[1]
a<-round(exp(a),4)
b <- round(coef(LWmodel)[2],2)

#########plot preparation#######
# Extract model statistics

r2 <- round(summary(LWmodel)$r.squared, 4)
n <- nrow(LWdata1)

eq <- paste0(
  "log(W) = ", a, " + ", b, " log(L)",
  "\nR² = ", r2,
  "\nn = ", n
)


#ggplot of log-transformed model
library(ggplot2)

p5 <- ggplot(LWdata1, aes(x = logL, y = logW)) +
  
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
    x = min(LWdata1$logL),
    y = max(LWdata1$logW),
    label = eq,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  
  labs(
    x = "log Length",
    y = "log Weight",
    title = ""
  ) +
  
  theme_bw(base_size = 14)

p5

#save plot
ggsave(
  filename = "Kawakawa_L-W.jpg",
  plot = p5,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

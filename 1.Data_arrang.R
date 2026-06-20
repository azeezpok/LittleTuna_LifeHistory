
setwd("/Users/abdulazeez/Documents/Conference and symposium/SS3_training/Biology data/E.affinis_bio")
library(ggpubr)
library(ggplot2)

# Read data
kawakawa <- read.csv("Kawakawa_data.csv")

# Check missing values
summary(kawakawa)

#########Build Regression Model#############
# Data with both TL and FL available
length_data <- na.omit(kawakawa[, c("TL_cm", "FL_cm")])

# Linear regression model
model1 <- lm(FL_cm ~ TL_cm, data = length_data)

# Model summary
summary(model1)

#Diagnostic Plots
par(mfrow = c(2,2))
plot(model1)

##Scatter Plot with Regression Line

ggplot(length_data, aes(x = TL_cm, y = FL_cm)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Total Length (cm)",
    y = "Fork Length (cm)",
    title = ""
  ) +
  theme_bw(base_size = 14)

#Display Regression Equation and R² on Plot
library(ggpubr)

p1<-ggplot(length_data, aes(TL_cm, FL_cm)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  
  stat_regline_equation(
    label.x.npc = "left",
    label.y.npc = 0.95
  ) +
  
  stat_cor(
    label.x.npc = "left",
    label.y.npc = 0.85
  ) +
  
  theme_bw()
#or
# Extract coefficients
intercept <- round(coef(model1)[1], 2)
slope <- round(coef(model1)[2], 3)

# Extract R²
r2 <- round(summary(model1)$r.squared, 3)

# Create equation label
eq <- paste0(
  "FL = ", intercept,
  " + ", slope, " × TL",
  "\nR² = ", r2
)

p2 <- ggplot(length_data, aes(TL_cm, FL_cm)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate(
    "text",
    x = min(length_data$TL_cm),
    y = max(length_data$FL_cm),
    label = eq,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Total Length (cm)",
    y = "Fork Length (cm)",
    title = ""
  )


#save plot
ggsave(
  filename = "Kawakawa_TL_FL_Regression1.jpg",
  plot = p2,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

########Predict Missing Fork Length Values####
# Identify rows with missing FL
missing_fl <- is.na(kawakawa$FL_cm)

# Predict FL from TL with round one decimal
kawakawa$FL_cm[missing_fl] <- round(
  predict(
    model1,
    newdata = kawakawa[missing_fl, ]
  ),
  1
)

# View updated data
head(kawakawa)

write.csv(kawakawa,"kawakawa_FL_imputed.csv",
  row.names = FALSE
)

#####create a new Date column by combining the Year, Month, and a fixed day (15)
# Create Date column
kawakawa$Date <- as.Date(
  paste(kawakawa$Year, kawakawa$Month, "15"),
  format = "%Y %B %d"
)

# Convert to d/m/y format
kawakawa$Date <- format(kawakawa$Date, "%d/%m/%Y")

# View result
head(kawakawa)

#save file
write.csv(kawakawa,"kawakawa_date_imputed.csv",
          row.names = FALSE)

###Preparing LF data for SS3##########
library(dplyr)

kawakawa_lfreq <- kawakawa %>%
  mutate(
    Unique_Trip = paste(State, Year, sep = "_")
  )
head(kawakawa_lfreq)

kawakawa_lfreq<-kawakawa_lfreq%>%dplyr::select(Year, FL_cm, Unique_Trip)
min(kawakawa_lfreq$Length)
max(kawakawa_lfreq$Length)

#To remove all rows where length is 10 cm or less
kawakawa_lfreq <- kawakawa_lfreq %>%
  dplyr::filter(FL_cm > 10)

kawakawa_lfreq1<-colnames(kawakawa_lfreq)<-c("Year","Length", "UniqueTrip")

head(kawakawa_lfreq1)

#save file
write.csv(kawakawa_lfreq1,"kawakawa_lfreq.csv",
          row.names = FALSE)

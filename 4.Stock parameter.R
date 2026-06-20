#_____________________________________________________________________________________________
###### Packages ##########
#_____________________________________________________________________________________________
library(dplyr)
library(TropFishR)
library(fishboot)
library(ks)
library(ggplot2)
# master branch
remotes::install_github(repo = "rschwamborn/fishboot",force = TRUE)

# dev branch (unstable)
remotes::install_github(repo = "rschwamborn/fishboot", ref = "dev")
library(fishboot)
#_____________________________________________________________________________________________
###### L-F arrangement ##########
#_____________________________________________________________________________________________
#source the data from 1. Data_arrang.R file
head(kawakawa)
kawa1<- kawakawa%>%dplyr::select(Date,FL_cm)
kawa1<-na.omit(kawa1)
colnames(kawa1)<-c("date","length")
kawa1$date<- as.Date(kawa1$date, format = "%d/%m/%Y")

#filter data
kawa1 <- kawa1 %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  filter(date >= as.Date("2021-01-15"))

kawa1$date <- format(kawa1$date, "%d/%m/%Y")
head(kawa1)

#View(mydata1)
ggplot(kawa1, aes(length)) +
  geom_histogram()

#To remove all rows where length is 10 cm or less
kawa1 <- kawa1 %>%
  filter(length > 10)

L_mean= mean(kawa1$length)
min(kawa1$length)
max(kawa1$length)

#jpeg("Histogram_length.jpg", res = 600,height = 6,width = 9,units = "in")
hist(kawa1$length, xlim = c(10,90),ylim = c(0,600),xlab= "Length (cm)",main = "")
abline(v=40.67,col="blue")
text(x = mean(kawa1$length), y = 600, label = paste("Mean L:", round(mean(kawa1$length), 2)))
dev.off()


meanL <- mean(kawa1$length, na.rm = TRUE)
Lp95 <- quantile(kawa1$length, probs = 0.95, na.rm = TRUE)

Lp95

p6 <- ggplot(mydata1, aes(x = length)) +
  geom_histogram(
    binwidth = 2,
    colour = "black",
    fill = "lightgrey"
  ) +
  
  geom_vline(
    xintercept = meanL,
    colour = "blue",
    linewidth = 1
  ) +
  
  annotate(
    "text",
    x = meanL + 5,
    y = 280,
    label = paste0("Mean L = ", round(meanL, 2), " cm"),
    colour = "blue",
    hjust = 0,
    size = 5
  ) +
  
  scale_x_continuous(
    limits = c(10, 90)
  ) +
  
  scale_y_continuous(
    limits = c(0, 300)
  ) +
  
  labs(
    x = "Length (cm)",
    y = "Frequency"
  ) +
  
  theme_bw(base_size = 14)

p6 

#save plot
ggsave(
  "Histogram_length.jpg", plot = p6,width = 8,
  height = 6,units = "in", dpi = 600)

########create LF
kawa1$date<- as.Date(kawa1$date, format = "%d/%m/%Y")

kawa1_2<- lfqCreate(data = kawa1, Lname = "length", Dname = "date", bin_size = 2,
                    species = "Euthynnus affinis", stock = "India")
plot(kawa1_2, Fname = "catch") #best one

kawa1_5<- lfqModify(kawa1_2, bin_size = 5) 
plot(affinis5, Fname = "catch")

kawa1_10<- lfqModify(kawa1_2, bin_size = 10) 
plot(affinis10, Fname = "catch")

#LF restructuring with MA
kawa1_5_ma3<-lfqRestructure(kawa1_5, MA=3)
kawa1_5_ma5<-lfqRestructure(kawa1_5, MA=5) 
kawa1_5_ma7<-lfqRestructure(kawa1_5, MA=7) #best one
kawa1_5_ma9<-lfqRestructure(kawa1_5, MA=9)
plot(kawa1_5_ma9,hist.sc = 0.75)
plot(kawa1_5_ma9, Fname = "rcounts")


#_____________________________________________________________________________________________
###### Linf and K estimation ##########
#_____________________________________________________________________________________________

####### Powell Wetherall plot
par(mfcol = c(1,1))
PW_kawa <- powell_wetherall(param = kawa1_5_ma7,
                            catch_columns = 1:ncol(kawa1_5_ma7$catch),
                            reg_int = c(5,10))
PW_kawa$Linf_est
PW_kawa$confidenceInt_Linf

#########K scan
Lmax<- max(kawa1$length)
Linf_cal<-Lmax/0.95

KScan_kawa <- ELEFAN(kawa1_5_ma7, method = "cross", Linf_fix = Linf_cal, 
                     K_range = seq(0.1,1.0,by=0.05),
                     cross.date = kawa1_5_ma7$dates[1],cross.midLength = kawa1_5_ma7$midLengths[1])
#MA=7, addl.sqrt = TRUE, hide.progressbar = TRUE)
#K scan with optimise method
KScan1_kawa <- ELEFAN(kawa1_5_ma7, method = "optimise", Linf_range = seq(65,85,by=0.5), 
                      K_range = seq(0.01,0.6,by=0.05))
# show results
KScan1_kawa$par
KScan1_kawa$Rn_max

###ELEFAN
#test the best model based on Rn value
set.seed(123)

res_SA <- ELEFAN_SA(kawa1_5_ma7)
res_GA <- ELEFAN_GA(kawa1_5_ma7)

res_SA$par
res_GA$par

res_SA$Rn_max
res_GA$Rn_max

## check how many fish >Linf
sum(kawa1$length > 65, na.rm = TRUE)
sum(kawa1$length > 70, na.rm = TRUE)
sum(kawa1$length > 75, na.rm = TRUE)

#SA performed better

res_SA_kawa <- ELEFAN_SA(kawa1_5, MA = 7, seasonalised = F, 
                         init_par = list(Linf = 70, K = 0.25, t_anchor = 0.5, C=0.5, ts = 0.5),
                         low_par = list(Linf = 67, K = 0.05, t_anchor = 0, C = 0, ts = 0),
                         up_par = list(Linf = 80, K = 0.5, t_anchor = 1, C = 1, ts = 1))

# show results
res_SA_kawa$par; res_SA_kawa$Rn_max


#ELEFAN_SA bootstrapped
lfq<-kawa1_5
MA <- 7
init_par <- NULL
low_par <- list(Linf = 65, K = 0.05, t_anchor = 0)
up_par <- list(Linf = 85, K = 0.7, t_anchor = 1)
SA_time <- 15
SA_temp <- 1e5
nresamp <- 500

# parallel version 
library(parallel)
t1 <- Sys.time()
res_kawa <- ELEFAN_SA_boot(lfq=lfq, MA = MA, seasonalised = F,
                           init_par = init_par, up_par = up_par, low_par = low_par,
                           SA_time = SA_time, SA_temp = SA_temp,
                           nresamp = nresamp, parallel = TRUE, 
                           seed = 1)
t2 <- Sys.time()
t2 - t1
res_kawa

# plot resulting distributions
univariate_density(res_kawa, use_hist = F)

# plot scatterhist of Linf and K
LinfK_scatterhist(res_kawa)

#jpeg("density_LinfK.jpg", res = 600,height = 6,width = 9,units = "in")
univariate_density(res_kawa, use_hist = F)
dev.off()

#plot the LF with res result
#par such as Linf, K & t_anchor are adapted from bootstrap; otherwise par=res_SA$par
#jpeg("LFcurve.jpg", res = 600,height = 6,width = 9,units = "in")
plot(res_SA_kawa,draw =F)
tmp <-lfqFitCurves(res_SA_kawa,par =list(Linf = 66.8, K = 0.29, t_anchor = 0.75),
                   col=4,lty=2,draw =TRUE)
legend("top",ncol=1,legend =c("estimated"),col=4,lty=1)
dev.off()

#######calculate t0#########
Linf <- 66.8
K <- 0.29

log_t0 <- -0.3922 -
  0.2752 * log10(Linf) -
  1.038  * log10(K)

t0 <- -10^(log_t0)

t0

############# assign estimates to the data list##############
kawa1_5_1 <- c(kawa1_5, res_SA_kawa$par)
#t0<-0 #not required in the TropfishR package as it gives any meaning
#coilia_wc1.0_1$t0=t0
###parameter updated from boostrap density
kawa1_5_1$Linf<-66.8
kawa1_5_1$K<-0.29
kawa1_5_1$t_anchor<-0.75
kawa1_5_1$phiL<-3.21
class(kawa1_5_1) <- "lfq"

#_____________________________________________________________________________________________
###### Estimate mortality & E ##########
#_____________________________________________________________________________________________

# estimation of M
tmax=t0+3/kawa1_5_1$K
tmax1<- res_SA_kawa$agemax

Ms <- M_empirical(Linf = kawa1_5_1$Linf, K_l = kawa1_5_1$K, method = "Then_growth")
Ms.pauly <- M_empirical(Linf = kawa1_5_1$Linf, K_l = kawa1_5_1$K,temp = 28, 
                        method = "Pauly_Linf")
M.list<-M_empirical(Linf = res_SA_kawa$par$Linf, K_l = res_SA_kawa$par$K,temp = 28, tmax = tmax,
                    method = c("Then_growth","Pauly_Linf","AlversonCarney","Hoenig",
                               "Then_tmax"))
Ms
Ms.pauly
M.list
kawa1_5_1$M <- as.numeric(Ms)
# show results M
paste("M =", as.numeric(Ms))

# run catch curve
# summarise catch matrix into vector and add plus group which is smaller than Linf
kawa1_5_2 <- lfqModify(kawa1_5_1, vectorise_catch = TRUE, plus_group = 67.5)
kawa1_5_2$catch<-as.matrix(rowMeans(kawa1_5_2$catch))

res_cc <- catchCurve(kawa1_5_2,  catch_columns = 1:ncol(kawa1_5_2$catch), 
                     reg_int = c(5,10), calc_ogive = T)
res_cc

#jpeg("Z_LCCC.jpg", res = 600,height = 6,width = 9,units = "in")
catchCurve(kawa1_5_2,  catch_columns = 1:ncol(kawa1_5_2$catch), 
           reg_int = c(5,10), calc_ogive = F)
dev.off()

#conver age to length est. pop of catch at 50% Length

library(ggplot2)

Lc50 <- round(res_cc$L50, 2)
Lc95<- round(res_cc$L95,2)

p7 <- ggplot(
  data.frame(Length = res_cc$midLengths,Selectivity = res_cc$Sest),
  aes(x = Length, y = Selectivity)
) +
  geom_line(colour = "red",linewidth = 1.2) +
  geom_vline(xintercept = res_cc$L50,colour = "blue",linewidth = 1,
    linetype = "dashed") +
  # Lm50 and Lm95
  geom_vline(
    xintercept = Lc50,
    colour = "blue",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = Lc95,
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
      "Lc50 = ", round(Lc50,2), " cm\n",
      "Lc95 = ", round(Lc95,2), " cm"
    ),
    size = 5
  ) +
  
  labs(x = "Fork Length (cm)", y = "Proportion") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14)
p7

#save plot
ggsave("Lc50.jpg",plot = p7,width = 9,height = 6,units = "in",
  dpi = 600)

#calculate F and E from Z and M
f=res_cc$Z-Ms
f
res_cc$Z
expl.rate=f/res_cc$Z
expl.rate

# assign the estimated parameters to the data
kawa1_5_2$M <- as.numeric(Ms)
kawa1_5_2$Z <- res_cc$Z
kawa1_5_2$currF <- as.numeric(kawa1_5_2$Z - kawa1_5_2$M)
kawa1_5_2$E <- kawa1_5_2$currF/kawa1_5_2$Z
kawa1_5_2$a<-0.0166
kawa1_5_2$b<-2.99
str(kawa1_5_2)


#_____________________________________________________________________________________________
###### Estimate VPA ##########
#_____________________________________________________________________________________________
set.seed(1000)
vpa <- VPA(param = kawa1_5_2,catch_columns = 1:ncol(kawa1_5_2$catch), 
           terminalF = kawa1_5_2$currF, terminalE = 0.5,analysis_type = "VPA", plot=T)
vpa
#plot
#jpeg("VPA.jpg", res = 600,height = 6,width = 9,units = "in")
VPA(param = kawa1_5_2,catch_columns = 1:ncol(kawa1_5_2$catch), 
    terminalF = kawa1_5_2$currF, terminalE = 0.5,analysis_type = "VPA", plot=T)
dev.off()


#_____________________________________________________________________________________________
###### TB model ##########
#_____________________________________________________________________________________________
vpa$FM_calc
kawa1_5_2$FM<-vpa$FM_calc

#add market value per kg for each length group
#rlfdata1$meanValue<-c(30,30,50,80,80,130,130,160,160,160,180,180,180,180,200,200)
#code for adding repeating values as alternative to above one
kawa1_5_2$meanValue<-c(rep(80,2),rep(120,2),rep(150,3),rep(180,3))
#coilia_wc1.0_2$meanValue<-NULL

# Thompson and Bell model with changes in F
TB1 <- predict_mod(kawa1_5_2, type = "ThompBell",
                   FM_change = seq(0,10,0.05), stock_size_1 = 1,FM_relative = T, 
                   curr.E = kawa1_5_2$E, plot = T, hide.progressbar = TRUE)

# Thompson and Bell model with changes in F and Lc
TB2 <- predict_mod(kawa1_5_2, type = "ThompBell",
                   FM_change = seq(0,5,0.1), Lc_change = seq(15,67.5,2),
                   stock_size_1 = 1,
                   curr.E = kawa1_5_2$E, curr.Lc = res_cc$L50,
                   s_list = list(selecType = "trawl_ogive",
                                 L50 = res_cc$L50, L75 = res_cc$L75),
                   plot = T, hide.progressbar = TRUE)

TB1_frm<-as.data.frame(cbind(TB1$meanB,TB1$totY,TB1$FM_change))
names(TB1_frm)<-c("Biomass","Yield","Fishing mortality")
TB1_frm
# plot results

#jpeg("TB_model1.jpg", res = 600,height = 7,width = 9,units = "in")
par(mfrow = c(2,1), mar = c(4,5,2,8)+0.1, oma = c(1,0,0,0))
plot(TB1, mark = TRUE)
mtext("(a)", side = 3, at = -1, line = 0.6)
plot(TB2, type = "Isopleth", xaxis1 = "FM", mark = TRUE, contour = 6)
mtext("(b)", side = 3, at = -0.1, line = 0.6)
dev.off()
# Biological reference levels
TB1$df_Es

# Current yield and biomass levels
TB1$currents
#_____________________________________________________________________________________________
###### THE END ##########
#_____________________________________________________________________________________________
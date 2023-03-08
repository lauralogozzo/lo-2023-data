# Fe-DOM Coupling L&O MS
# Created: 12/21/2022
# Modified: 12/21/2022
# Author: Laura Logozzo

# ----------- # Run End-Member Mixing Analysis (EMMA) on PARAFAC components ------------

# Following: https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/End-member-modelling-analysis/The-EMMA-algorithm/The-EMMAgeo-package/index.html

## Steps:
# Check the data
# Bracket weight transformation limits
# Bracket numbers of end-members
# Calculate and extract robust end-members
# Run EMMA with optimal parameters
# Evaluate model quality
# Estimate model uncertainty

# Load packages

library("EMMAgeo")
library("data.table")
library("ggplot2")
library("ggpubr")
library("xlsx")
library("dplyr")
library("tidyr")
library("reshape2")
library("finalfit")
source("https://lauralogozzo.github.io/assets/myCols.R.txt") # Custom functions for picking colors
source("https://lauralogozzo.github.io/assets/stat_cor_r2.R.txt") # Custom function for adding R2 to ggplot

# Set GGplot Aesthetics

theme_set(theme_light() + 
            theme(
              panel.grid = element_blank(),
              panel.border = element_rect(color = "black"),
              axis.title = element_text(size = 14),
              axis.text = element_text(color = "black", size = 14),
              axis.ticks = element_line(color = "black"),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12)
            )
)

# SET TO DATA WORKING DIRECTORY

setwd("") # Navigate to data folder lo-2023-data

# Load data

d <- read.csv("01-EMMA/Fe-DOM-data.csv") # Data from 01-EMMA

# Create dataframe with just PARAFAC components

emma_data <- d[,c("DT","Q_cms","C1_0.2_per","C2_0.2_per","C3_0.2_per","C4_0.2_per","C5_0.2_per")]

# Create metadata file linking sample id to sites and dates

emma_data$sample_id <- 1:nrow(emma_data)
emma_data_meta <- emma_data[,c("sample_id","DT","Q_cms")]
emma_data <- subset(emma_data, select = -c(DT,Q_cms))

# Remove missing values

emma_data <- emma_data[!is.na(emma_data$C1_0.2),]

# Set rownames to be sample id

rownames(emma_data) <- emma_data$sample_id

emma_data <- emma_data %>%
  dplyr::select(-sample_id)

# Data needs to be as a matrix:

emma_data <- data.matrix(emma_data)
colnames(emma_data) <- 1:5

# Divide by 100 to have values equal to 1 (not strictly necessary)

emma_data <- emma_data/100

# Data needs to be numeric. Is data numeric?

is.numeric(emma_data)

# Are there missing values?

sum(is.na(emma_data)) # 0 means NO missing values!
sum(apply(emma_data, MARGIN = 2, FUN = sum) == 0) # 0 means NO missing values!
apply(emma_data, MARGIN = 1, FUN = sum) # All rows should sum to the same value (in this case 1)

## Plot the Data

par(mfrow = c(1,1))
plot(rep(1:5,nrow(emma_data)),emma_data[,1:5], type = 'n')
for(i in 1:nrow(emma_data)){
  lines(1:5,emma_data[i,])
}
legend("topleft",lty = 1, "Samples")

# number of samples in analysis?

dim(emma_data)[1]

## Bracket transformation weight limits
# A column-wise weight transformation is applied by scaling the columns based on percentiles with lower (l) and upper (100−l) boundaries as weights.

l.seq <- seq(from = 0, to = 0.5, by = 0.01)
l <- test.l(emma_data, l = l.seq); l
l.max <-l$l.max # Max value that the transformation remains stable
l.max

## Bracket numbers of end members
# i.e., HOW MANY ENDMEMBERS ARE THERE? This is generally unknown!

seq.q <- seq(from = 2, to = 5, by = 1) # set sequence of end-members
seq.l <- seq(from = 0, to = l.max, by= 0.01) # sequence of weight transformations

params <- test.parameters(emma_data, 
                          q = seq.q, 
                          l = seq.l, 
                          rotation = "Varimax",
                          plot = "mRt", 
                          # "mRt" (mean relative total error), 
                          # "mEm" (mean absolute row-wise error), 
                          # "mEn" (mean absolute column-wise error), 
                          # "mRm" (mean relative row-wise error), 
                          # "mRn" (mean relative column-wise error), 
                          # "ol" (number of overlapping end-members). 
                          legend = 'bottomright', 
                          cex = 1)

# Request optimal number of end-members

num.end.members <- params$q.max
num.end.members # 2
q <- max(num.end.members, na.rm = T)

# Check data matrix for consistency

check.data(X = emma_data, num.end.members, seq.l, c = 1)

## Calculate and extract robust end-members
# Tests all possible combinations of end members and weight tansformation limits iteratively

xlabs <- c(expression(paste("Component Number",
                            sep = "")),
           expression(paste("Component Number",
                            sep = "")))

seq.l <- seq(from = 0, to = l.max, by = 0.02)
robust.params <- test.robustness(X = emma_data,
                                 q = unique(na.exclude(num.end.members)), # end-members
                                 l = seq.l, # weight transformation limits
                                 plot = TRUE,
                                 colour = c('blue', 'green'), 
                                 xlab = xlabs)

# Distinguish between modes

hist(robust.params$modes,
     breaks = 10,
     main = "Mode positions",
     xlab = "Class", 
     xaxt = 'n')
axis(side = 1, at = seq(0.5,5.5,0.5))

# Get limits

# Extract end-member loadings whose modes fall into specified limits

limits = cbind(c(1, 4.5), # lower limit: first bar start, second bar start (on x axis)
               c(1.5, 5.0)) # upper limit: first bar end, second bar end (on x axis)

rEM <- robust.EM(robust.params, Vqsn = robust.params$Vqsn,  # samples 
                 limits = limits,
                 Vqn = robust.params$Vqn, # normalised factor loadings
                 plot = TRUE,
                 legend = "topright",
                 cex = 0.85,
                 colour = c("red","blue"),
                 median = TRUE)

Vqn.mean <- rEM$Vqn$mean # robust end-members mean
Vqn.sd <- rEM$Vqn$sd # robust end-members standard deviation
Vqn.res <- residual.EM(Vqn.mean) # residual end-member

## Plot the variability explained by the end members

# Create dataset with results

Vqn <- rbind(Vqn.mean, Vqn.res)

# Rename rows

rownames(Vqn) <- c("EM 1", "EM 2","Residuals")
rownames(Vqn.sd) <- c("EM 1","EM 2")

# Reshape

Vqn <- reshape2::melt(Vqn)
Vqn.sd <- reshape2::melt(Vqn.sd)

# Rename columns

colnames(Vqn) <- c("EM","Component","Mean")
colnames(Vqn.sd) <- c("EM","Component","SD")

# Make PARAFAC component #s factors

Vqn$Component <- factor(Vqn$Component)
Vqn.sd$Component <- factor(Vqn.sd$Component)

# Merge datasets

Vqn <- merge(Vqn, Vqn.sd, by = c("EM","Component"),all.x = T)

# Get color palette

myCols <- getCols("dune")

# Plot with residuals

ggplot(data = Vqn, aes(x=Component,y=Mean, fill = factor(EM))) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), position = position_dodge(width=0.9), width = 0.25) +
  scale_fill_manual(values = c(myCols[c(5,6)],"grey60"), name = "") +
  xlab("Component") +
  ylab("Loadings")

# Find optimal weight transformation limit

l_opt <- get.l.opt(X = emma_data, 
                   l = seq.l, 
                   quality = "mRt", 
                   Vqn = Vqn.mean, 
                   plot = TRUE)[1]
# Not much difference across the weight transformations. So will use l = 0

## Run EMMA with optimal parameters

EM <- EMMA(X = emma_data, #samples
           q = nrow(Vqn.mean),  # number of end-members
           l = 0, 
           Vqn = Vqn.mean, # robust end-member loadings
           plot = TRUE,
           EM.ID = c("EM 1","EM 2"),
           xlab = c(expression(paste("Component Number")), "Sample ID"),
           cex = 0.7,
           col = myCols[c(5,7)])

# Percent sample-wise variance explained is poor, but we can evaluate against the flux from GPP below to evaluate the model

# Save model results

save(EM, file = "Data/EMMA/EM.Rdata")

## Evaluate model quality

# Examine variance

EM$Mqs.var / 100

# RMSE
EM$RMSEm

# Convert results to long format for plotting

emma.wide <- data.frame(EM$loadings)
colnames(emma.wide) <- c("C1","C2","C3","C4","C5")
emma.wide$EM <- c("Allochthonous-like","Autochthonous-like")

emma.long <- emma.wide %>%
  tidyr::gather(Components,Value,-c(EM))

emma.long$sd <- matrix(rEM$loadings$sd)

# Plot Loadings

myCols <- getCols("dune")
ggplot(data=emma.long, aes(y=Value, x = Components, fill = EM)) + 
  geom_bar(stat = "identity", position = "dodge2") +
  scale_fill_manual(values = c("Allochthonous-like" = myCols[5],"Autochthonous-like" = myCols[6]), name = "End-Members") +
  geom_errorbar(aes(ymin = Value - sd, ymax = Value + sd), position = position_dodge(width=0.9), width = 0.25) +
  ylab("Normalized Loading") +
  theme(legend.position = c(0.5,0.8))

# Summarize scores

scores <- data.frame(EM[["scores"]],row.names = rownames(emma_data)); colnames(scores) <- c("Allo_per","Auto_per")
complete_data <- merge(emma_data_meta,scores, by.x = "sample_id", by.y = "row.names", all = T)

# Merge with original dataframe

d <- merge(d,complete_data[,c("DT","Allo_per","Auto_per")], by = "DT")

# Calculate Allo and Auto Concentrations from Percentages

d$Allo <- d$Allo_per * d$DOC0.2
d$Auto <- d$Auto_per * d$DOC0.2

# Compare to Auto GPP to evaluate model

ggplot(data=d, aes(x=DOC.GPP, y=Auto)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor_r2()
# Fairly good fit

# Create new dataset for input into LOADEST

loadest_input <- d[,c("DT","Q_cfs","Fe0.02mgL","FeCmgL","DOC0.2","Allo")] # Include only columns to run in loadest
loadest_input$CDATE <- strftime(d$DT,"%Y%m%d") # Reformat date
loadest_input$CTIME <- strftime(d$DT,"%H%M") # Reformat time
loadest_input$DOC0.2 <- round_tidy(loadest_input$DOC0.2,1)
loadest_input$Allo <- round_tidy(loadest_input$Allo,1)
loadest_input$Fe0.02mgL <- round_tidy(loadest_input$Fe0.02mgL,3)
loadest_input$FeCmgL <- round_tidy(loadest_input$FeCmgL,3)
loadest_input <- subset(loadest_input, select = -c(DT)) # Remove original Date and Time column
loadest_input[loadest_input == "NA"] <- round_tidy(-9999,0)
loadest_input <- loadest_input[,c("CDATE","CTIME","Q_cfs","Fe0.02mgL","FeCmgL","DOC0.2","Allo")]
colnames(loadest_input) <- c("CDATE","CTIME","CFLOW","sFe","cFe","doc","allo")

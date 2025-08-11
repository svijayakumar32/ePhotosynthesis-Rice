# Set the working directory by choosing a folder interactively

# Install the rChoiceDialogs package
install.packages("rChoiceDialogs")

# Load rChoiceDialogs
library(rChoiceDialogs)

# Set the working directory
setwd(rchoose.dir(default = getwd(),
            caption = "Select Directory"))

# Install A/Ci fitting package if necessary

# msuRACiFit is a package of R scripts used to fit CO2 response curves developed by Thomas Sharkeys lab at Michigan State University.

# Source Code:
# https://github.com/poales/msuRACiFit

# Publication:
# <https://doi.org/10.1111/pce.14153>

# Install msuRACiFit package using devtools or remotes:

library(devtools)
devtools::install_github("poales/msuRACiFit") # OR
# remotes::install_github("poales/msuRACiFit")

# Load the here, msuRACiFit and readxl libraries:

library(here)
library(readxl)
library(msuRACiFit)

# Using 'here' roots the file path to the R environment.

# Specify an output directory:

output_path <- here("Outputs", "Parameters")

# Import csv files containing gas exchange data:

A330502WT1 <- read.csv(here("Data","IR64-A009-07-33-05-02_Wildtype1.csv"))
A330502WT2 <- read.csv(here("Data","IR64-A009-07-33-05-02_Wildtype2.csv"))
A330504WT1 <- read.csv(here("Data","IR64-A009-07-33-05-04_Wildtype1.csv"))
A330504WT2 <- read.csv(here("Data","IR64-A009-07-33-05-04_Wildtype2.csv"))
A330506WT1 <- read.csv(here("Data","IR64-A009-07-33-05-06_Wildtype1.csv"))
A330506WT2 <- read.csv(here("Data","IR64-A009-07-33-05-06_Wildtype2.csv"))
A330509WT1 <- read.csv(here("Data","IR64-A009-07-33-05-09_Wildtype1.csv"))
A330509WT2 <- read.csv(here("Data","IR64-A009-07-33-05-09_Wildtype2.csv"))

# Store the files in a list:

datasets <- list(c(A330502WT1),
                 c(A330502WT2),
                 c(A330504WT1),
                 c(A330504WT2),
                 c(A330506WT1),
                 c(A330506WT2),
                 c(A330509WT1),
                 c(A330509WT2))

# Gamma Star can be calculated for different temperatures using the calcGammaStar function.

# msuRACiFit uses preset in-vitro values from Hermida-Carrera et al (2016), where everything is presented in ppm, except for N. tabacum and A. thaliana.

# Since the A/Ci responses are analysed in Pa, the resulting Gamma\*, Kc and Ko must be adjusted from ppm to Pa units before running calculations.

# $c$ and $\Delta H_a$ for Gamma\* presented in the app/ preset table are for ppm, and Gamma\* is corrected as follows:

# Gamma* (Pa) = Gamma*(ppm) * 1000 / 1000000 * Patm,
# Kc (Pa) = Kc(umol mol -1) * 1000 / 1000000 * Patm
# Kcair (Pa) = Kcair(umol mol -1) * 1000 / 1000000 * Patm,
# Ko (kPa) = O(Pa) / ((Kcair(Pa) / Kc(Pa)) - 1)/1000

# Create output vectors of mean and SEM for temperature and pressure:

mean_temps <- numeric(length = 8) 
mean_press <- numeric(length = 8)
sem_mean_temps <- numeric(length = 8) 
sem_mean_press <- numeric(length = 8)

# Calculate mean and SEM for temperature and pressure:

for (i in seq_along(datasets)) { 
  mean_temps[i] <- mean(datasets[[i]]$Tleaf) 
  mean_press[i] <- mean(datasets[[i]]$Press)
  sem_mean_temps[i] <- sd(mean_temps) / sqrt(length(mean_temps))
  sem_mean_press[i] <- sd(mean_press) / sqrt(length(mean_press))
}

# Call the calcgammastar function with different temperature inputs and get outputs saved in a single vector for 8 files:

Gstars_rice<-msuRACiFit::calcGammaStar(13.7,24.6,mean_temps,21) # Rice c and dHa at leaf temperature
sem_Gstars_rice<- sd(Gstars_rice)/ sqrt(length(Gstars_rice))
Gstars_rice_25<-msuRACiFit::calcGammaStar(13.7,24.6,25,21) # Rice c and dHa at 25C

# Next, account for the variation in Press for each dataset and converting GammaStar back to Pa:

Gstars_Pa_rice<-Gstars_rice*(1000 / 1000000)*mean_press 
sem_Gstars_Pa_rice<- sd(Gstars_Pa_rice)/ sqrt(length(Gstars_Pa_rice))

# Save these gamma star values in a .csv:

Gstars_rice_file <- file.path(output_path, "Gstars_rice.csv")
write.csv(Gstars_rice,Gstars_rice_file) #ppm
Gstars_Pa_rice_file <- file.path(output_path, "Gstars_Pa_rice.csv")
write.csv(Gstars_Pa_rice,Gstars_Pa_rice_file) #Pa

# Alternatively, to fit data using the Shiny app:
# Run msuRACiFit::genApp() 
## Select *O.sativa* as preset plant
## Fill in the average Tleaf and Pa for each curve
## Generate guesses and fits
## Save outputs as .png and csv tables

# Calculate Kc and Kcair for our Tleaf values in umol mol-1, then convert them into Pa before working out Ko.

# Make sure all constants are named differently within the function to what they are called outside the function.

# Import data for Kc_c and Kc_dHa (from Hermida-Carrera et al, 2016):

Kc_c <- 38.9   # unitless
Kc_dHa <- 83.1 # kJ mol -1

# Import data for Kcair_c and Kcair_dHa from Hermida-Carrera et al, 2016:

Kcair_c <- 30.5   # unitless
Kcair_dHa <- 60.5 # kJ mol -1

# Write function to use the mean temperature in celsius from each curve for calculating Kc and Kcair:

kinetic_params <- function(tempinC) 
{ ideal_gas <- 0.008314 # Ideal gas constant 
kelvin <- 273.15 # Kelvin constant 
Kc_scaling <- Kc_c # Kc scaling constant (unitless)
Kc_activation <- Kc_dHa # Kc energy of activation (kJ mol -1) 
Kcair_scaling <- Kcair_c # Kcair scaling constant (unitless) 
Kcair_activation <- Kcair_dHa # Kcair energy of activation (kJ mol -1) 
O2 <- 21000 # Oxygen concentration (Pa) 
# Calculate Kc in umol at a specific temperature 
Kc <- exp(Kc_scaling - (Kc_activation / (ideal_gas * (kelvin + tempinC)))) 
# Calculate Kcair in umol at a specific temperature 
Kcair <- exp(Kcair_scaling - (Kcair_activation / (ideal_gas * (kelvin + tempinC))))      
return(c(Kc,Kcair)) }

# Create an empty dataframe to store all Kc and Kcair outputs:

all_param <- data.frame(Kc = numeric(length(mean_temps)), 
                        Kcair = numeric(length(mean_temps)))

# Loop through kinetic_params function, calculate Kc and Ko values at each temperature in mean_temps and convert values from ppm to Pa according to each mean_press:

for (j in seq_along(mean_temps)) { 
  # Call the kinetic_params function to calculate Kc and Kcair for the current temperature and convert them to Pa based on mean pressure 
  param <- kinetic_params(mean_temps[j])*(1000 / 1000000)*mean_press[j] 
  # Store the results in the output dataframe 
  all_param[j, ] <- param }

# Currently Kc and Kcair are in Pa

# Call the kinetic_params function to calculate $K_{cair}$ in $umol \ mol^{-1}$ :

Kcair_umol <- vector("list", length(mean_temps)) 
for (j in seq_along(mean_temps)) { 
  Kcair_umol[[j]] <- kinetic_params(mean_temps[j])[2]}

# Save Kcair and temperature values in .csv files:

Kcair_file <- file.path(output_path, "Kcair_umol_rice.csv")
temp_file <- file.path(output_path, "temps_rice.csv")

write.csv(Kcair_umol, Kcair_file)
write.csv(mean_temps,temp_file)

# Append a column for calculating Ko from Kc and Kcair and converting these values to kPa:

all_param$Ko<-(21000 / ((all_param$Kcair / all_param$Kc) - 1))/1000 
all_param$Temperature <- mean_temps 
all_param$Press <- mean_press
mean_Kc<-mean(all_param$Kc)
mean_Ko<-mean(all_param$Ko)
mean_Kcair<-mean(all_param$Kcair)
sem_Kc <- sd(all_param$Kc)/ sqrt(length(all_param$Kc))
sem_Ko <- sd(all_param$Ko)/ sqrt(length(all_param$Ko))
sem_Kcair <- sd(all_param$Kcair)/ sqrt(length(all_param$Kcair))

# Save parameters in a file:

all_param_file <- file.path(output_path, "kinetic_parameters_rice.csv")
write.csv(all_param,all_param_file)

# Input Kc and Ko values from all_param into fitComplete

# Call fitComplete with different inputs and use a list of lists to store outputs

# Initialise an empty list to store output of fitComplete:

complete_Fits_rice_list <- vector("list", length(datasets)) 

# This loop indexes through vectors (mean_temps,mean_press,Gstars_Pa_rice) and lists (assimilation_rates, intercellular_CO2) to fit 8 curves using fitComplete and save their outputs

# Vector of forceValues/bound_l/bound_h in order: VcMax, J, TPU, gm, Rd, aG, aS

for (k in seq_along(datasets)) {  
  data_k <- datasets[[k]]
  if (hasName(data_k, "data")) {
    data_k <- data_k$data
  }
  data_k <- as.data.frame(data_k)
  
  complete_Fits_rice <- fitComplete(
    data = data_k,  
    name_assimilation = "Photo",
    name_ci = "Ci",
    gammastar = Gstars_Pa_rice[k],
    O2 = 21,
    pressure = mean_press[k],
    tleaf = mean_temps[k],
    initialGuess = NA,
    forceValues = c(NA, NA, NA, NA, NA, NA, NA),
    bound_l = c(1, 1, 1, .001, .001, 0, 0),
    bound_h = c(1000, 1000, 1000, 30, 10, 1, .75),
    ignoreTPU = F,
    maxiter = 500,
    Kc = all_param[k, "Kc"],
    Ko = all_param[k, "Ko"]
  )
  
  complete_Fits_rice_list[[k]] <- complete_Fits_rice
}

# Save parameter outputs:

for (o in seq_along(complete_Fits_rice_list)) {
  filename <- sprintf("completeFits_rice_params_%d.csv", o)
  output_filename <- file.path(output_path,filename)
  write.csv(complete_Fits_rice_list[[o]][[1]], 
            file = output_filename, 
            row.names = FALSE)
}

# Calculate $g_{m}$ in $umol \ m^{-2} \ s^{-1}$:

gm_umol <- vector("list", length(mean_press)) 
for (j in seq_along(mean_press)) { 
  gm_umol[[j]] <- complete_Fits_rice_list[[j]][[1]]$gm*(mean_press[j]/1000)}

# Save gm values in a file:

gm_file <- file.path(output_path, "gm_umol_rice.csv")
write.csv(gm_umol,gm_file)

# Initialise an empty list to store table outputs and sums of squared residuals:

complete_Tables_rice_list <- vector("list", length(datasets))
SSR <- numeric(length(complete_Tables_rice_list))

# Produce data tables using parameters fitted for each curve and the SSR associated with each fit:

for (k in seq_along(datasets)) {  
  data_k <- datasets[[k]]
  if (hasName(data_k, "data")) {
    data_k <- data_k$data
  }
  data_k <- as.data.frame(data_k)
  
  complete_Tables_rice <- reconstituteTable(
    data = data_k,  
    fitParams = complete_Fits_rice_list[[k]][[1]],
    name_assimilation = "Photo",
    name_ci = "Ci",
    gammastar = Gstars_Pa_rice[k],
    O2 = 21,
    pressure=mean_press[(k)], 
    tleaf = mean_temps[(k)], 
    ignoreTPU=F, 
    Kc=all_param[k, "Kc"], 
    Ko=all_param[k, "Ko"])
  
  complete_Tables_rice_list[[k]] <- complete_Tables_rice
  SSR[k] <- sum(complete_Tables_rice_list[[k]]$`res^2`)
}

# Save data tables in separate files:

for (p in seq_along(complete_Tables_rice_list)) {
  filename <- sprintf("completeFits_rice_dataTable_%d.csv", p)
  output_filename <- file.path(output_path,filename)
  write.csv(complete_Tables_rice_list[[p]], 
            file = output_filename, 
            row.names = FALSE)
}

# Extract the parameter values from each list and combine them in a single dataframe, adding an identifier column 1-9 (based on k) as the first column:

combined_results <- do.call(rbind, 
                            lapply(seq_along(complete_Fits_rice_list), function(k) {
                              if(!is.null(complete_Fits_rice_list[[k]][[1]])) {
                                df <- complete_Fits_rice_list[[k]][[1]]  
                                df$Dataset_ID <- k  
                                df <- df[, c("Dataset_ID", 
                                             setdiff(names(df), 
                                                     "Dataset_ID"))]  
                                return(df)
                              } else {
                                return(NULL)
                              }
                            }))

# Compute mean and SEM of parameter values:

VcMax <- mean(combined_results$VcMax)
J <- mean(combined_results$J)
TPU <- mean(combined_results$TPU)
Rd <-  mean(combined_results$rL)
gm <-  mean(combined_results$gm)
sem_VcMax <- sd(combined_results$VcMax)/ sqrt(length(combined_results$VcMax))
sem_J <- sd(combined_results$J)/ sqrt(length(combined_results$J))  
sem_TPU <- sd(combined_results$TPU)/ sqrt(length(combined_results$TPU)) 
sem_Rd <- sd(combined_results$rL)/ sqrt(length(combined_results$rL))
sem_gm <- sd(combined_results$gm)/ sqrt(length(combined_results$gm))   

# Save combined parameter values:

parameters_file <- file.path(output_path, "complete_Fits_rice_params_all.csv")
write.csv(combined_results,parameters_file)

# Save SSR values:

SSR_file <- file.path(output_path, "SSR_values_rice.csv")
write.csv(SSR,SSR_file)

# Initialise an empty list to store graph outputs:

complete_Graphs_rice_list <- vector("list", length(datasets))

# Produce graphs of limitations using parameters fitted for each curve:

for (k in seq_along(datasets)) {  
  data_k <- datasets[[k]]
  if (hasName(data_k, "data")) {
    data_k <- data_k$data
  }
  data_k <- as.data.frame(data_k)
  
  complete_Graphs_rice <- reconstituteGraph(
    data = data_k,  
    fitParams = complete_Fits_rice_list[[k]][[1]],
    name_assimilation = "Photo",
    name_ci = "Ci",
    gammastar = Gstars_Pa_rice[k],
    O2 = 21,
    pressure=mean_press[(k)], 
    tleaf = mean_temps[(k)], 
    ignoreTPU=F, 
    Kc=all_param[k, "Kc"], 
    Ko=all_param[k, "Ko"])
  
  complete_Graphs_rice_list[[k]] <- complete_Graphs_rice
}

# Save graphs:

for (p in seq_along(complete_Graphs_rice_list)) {
  filename <- sprintf("completefits_rice_Graph_%d.png", p)
  output_filename <- file.path(output_path, filename)
  
  png(output_filename, width = 519, height = 457)
  plot(complete_Graphs_rice_list[[p]])
  dev.off()
}

# Create output for limitation data:

limitations_list <- vector("list", length(complete_Graphs_rice_list))

# Export plot coordinates for limitations:

for (q in seq_along(complete_Graphs_rice_list)) {
  graph <- complete_Graphs_rice_list[[q]]
  Ac <- graph$plot_env$cdat
  Aj <- graph$plot_env$jdat
  Ap <- graph$plot_env$pdat
  
  limitations_list[[q]] <- list(Ac, Aj, Ap)
}

# Write limitations to output file:

limitations_file <- file.path(output_path, "completeFit_rice_limitations.csv")
write.csv(limitations_list,limitations_file)

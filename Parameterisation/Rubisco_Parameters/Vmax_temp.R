# Set the working directory by choosing a folder interactively

# Install the rChoiceDialogs package
install.packages("rChoiceDialogs")

# Load rChoiceDialogs
library(rChoiceDialogs)

# Set the working directory
setwd(rchoose.dir(default = getwd(),
            caption = "Select Directory"))

# Load libraries:

library(ggplot2)
library(minpack.lm)

# Bernacchi temperature functions

# Currently the e-Photosynthesis model sets PrV111 = PsV1 * 0.24, according to the known value for tobacco at 25°C (Whitney et al., 1999).
# However, this ratio varies across different crop species and temperatures.
# Seeing as there is no existing temperature functions for rice, we will use the equations by Bernacchi (2001) and values from Makino (1988) to compare the temperature responses of the parameters.

# Specify constants:

T25 <- 25 + 273.15                 # Reference temperature in K 
Tp <- 28.9310407291759 + 273.15    # Average measurement temperature in K
R <- 0.008314                     # Ideal gas constant in kJ K-1 mol-1
c_Vc <- 26.35                       # Scaling constant, Vcmax             (N. tabacum L. cv. W38, Bernacchi et al 2001)
dHa_Vc <- 65.33                    # Activation energy, Vcmax, kJ mol -1 (N. tabacum L. cv. W38, Bernacchi et al 2001)
c_Vo <- 22.98                       # Scaling constant, Vomax             (N. tabacum L. cv. W38, Bernacchi et al 2001)
dHa_Vo <- 60.11                     # Activation energy, Vomax, kJ mol -1 (N. tabacum L. cv. W38, Bernacchi et al 2001)

# Use the Arrhenius equation to calculate the Vo/Vc ratio at 25°C for N.tabacum (c and dHa from Bernacchi et al., 2001):
# Parameter = exp(c- dHa/ R * T_k)

Bernacchi_Vc_25 <- exp(c_Vc - dHa_Vc / (R * T25))
Bernacchi_Vo_25 <- exp(c_Vo - dHa_Vo / (R * T25))
Bernacchi_PrPs_ratio_25 <- Bernacchi_Vo_25 / Bernacchi_Vc_25
print(Bernacchi_PrPs_ratio_25)

# Calculate Vo/Vc ratios across a range of temperatures (10-40°C) :

Temp_Range <- data.frame(Temp_C = 10:40)                       # Temperatures in °C
Temp_Range$Temp_K <- Temp_Range$Temp_C + 273.15                # Temperatures in Kelvin
Temp_Range$Vc <- exp(c_Vc - dHa_Vc / (R * Temp_Range$Temp_K))  # Vc values
Temp_Range$Vo <- exp(c_Vo - dHa_Vo / (R * Temp_Range$Temp_K))  # Vo values
Temp_Range$Vo_Vc_ratio <- Temp_Range$Vo / Temp_Range$Vc        # Vo/Vc ratios

# Plot Vo/Vc data:

ggplot(Temp_Range, aes(x = Temp_C, y = Vo_Vc_ratio)) +
                   geom_line(color = "blue") +
                   geom_point(color = "blue", size = 2) +
                   labs(x = expression("Temperature ("*degree*C*")"), y = "Vo/Vc") +
                   theme_minimal()
ggsave("Outputs/Temp_vs_VcVo.png")

# Plot Vo over temperature:

ggplot(Temp_Range, aes(x = Temp_C, y = Vo)) +
                   geom_line(color = "green") +
                   geom_point(color = "green", size = 2) +
                   labs(x = expression("Temperature ("*degree*C*")"), y = "Vo") +
                   theme_minimal()
ggsave("Outputs/Temp_vs_Vo.png")

# Plot Vc over temperature:

ggplot(Temp_Range, aes(x = Temp_C, y = Vc)) +
                   geom_line(color = "red") +
                   geom_point(color = "red", size = 2) +
                   labs(x = expression("Temperature ("*degree*C*")"), y = "Vc") +
                   theme_minimal()
ggsave("Outputs/Temp_vs_Vc.png")

# Calculate Bernacchi values at measurement temperature:

Bernacchi_Vc_Tp <- exp(c_Vc - dHa_Vc / (R * Tp))
Bernacchi_Vo_Tp <- exp(c_Vo - dHa_Vo / (R * Tp))
Bernacchi_PrPs_ratio_Tp <- Bernacchi_Vo_Tp / Bernacchi_Vc_Tp
print(Bernacchi_PrPs_ratio_Tp)

# Calculate Makino values at 25°C:

Makino_Vc_25 <- 1.77                                # Carboxylase activity (µmol (mg enzyme)⁻¹ min⁻¹)
Makino_Vo_25 <- 0.58                                # Oxygenase activity (µmol (mg enzyme)⁻¹ min⁻¹)
Makino_PrPs_ratio_25 <- Makino_Vo_25 / Makino_Vc_25 # Vo/Vc ratio
print(Makino_PrPs_ratio_25)

# Adjusting the value Vcmax at 25°C to measurement temperature Tp

# We can use one of two approaches to calculate the Makino value of Vo/Vc at Tp:

# 1)  Use the ratios of Makino and Bernacchi at 25°C and multiply by Bernacchi parameter estimate at a given temperature
# 2)  Fit an Arrhenius equation to a non-linear model to find the Bernacchi temperature response of (plotted above)

# Use linear scaling of Makino to Bernacchi ratios:

Makino_PrPs_ratio_Tp <- Makino_PrPs_ratio_25*Bernacchi_PrPs_ratio_Tp/Bernacchi_PrPs_ratio_25
print(Makino_PrPs_ratio_Tp)

# Non-linear fits for Vo/Vc

# Define predictor and response variables:

T <- Temp_Range$Temp_K      # Temperature in K
y <- Temp_Range$Vo_Vc_ratio # Vo/Vc ratio

# So far, we have defined the Arrhenius equation as:

# Parameter = exp (c- dHa/R * T_k)

# c is meant to scale the amplitude of the function whereas dHa represents the temperature dependence relationship.

# Using the term c−dHa combined these two effects into one, making them harder to interpret individually.
# Without c, the model assumes that dHa depends solely on temperature and ignores other potential variations.

# This means the equation can also be rewritten as:

# y = c * exp(-dHa/ R * T_k)

# to ensure that the exponent has a negative sign.

# Define an exponential model with non-linear parameters using this form of the Arrhenius equation.

# Create a fit type for a custom nonlinear model using nlsLM()
 
# nlsLM() uses the Levenberg-Marquardt algorithm for nonlinear least-squares estimation

# Compare the values by re-scaling the Vo/Vc ratio from Bernacchi so that it has a value of 1 at 25°C and then fitting the Arrhenius function using these values.
# Divide all Vo/Vc values calculated at different temperatures by the value at 25°C:

normalised_ratios <- Temp_Range$Vo_Vc_ratio/Bernacchi_PrPs_ratio_25;

# Obtain Arrhenius constants using the re-scaled values:
nls_norm_model <- nlsLM(normalised_ratios ~ c_VoVc * exp(-dHa_VoVc / (R * T)),
                                            data = Temp_Range,
                                            start = list(c_VoVc = 0, 
                                                         dHa_VoVc = 2))
summary(nls_norm_model)

# Extract the coefficients from this second model:

norm_coefficients <- coef(nls_norm_model)
c_norm <- as.numeric(norm_coefficients["c_VoVc"])
dHa_norm <- as.numeric(norm_coefficients["dHa_VoVc"])

# Predict the Bernacchi ratio at Tp using the parameters obtained from the non-linear fit of normalised values:

Bernacchi_PrPs_ratio_Tp <- (c_norm * exp(-dHa_norm / (R * Tp)))
print(Bernacchi_PrPs_ratio_Tp)

# Multiply the Bernacchi ratio at Tp by the Makino ratio at 25°C:

Makino_PrPs_ratio_Tp2 <- Makino_PrPs_ratio_25*Bernacchi_PrPs_ratio_Tp
print(Makino_PrPs_ratio_Tp2)

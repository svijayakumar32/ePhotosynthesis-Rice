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

# Data Extraction
# To modify the regulation of Rubisco in the e-Photosynthesis model, we need to fit a model to be able to predict the response of activation state in rice to changes in temperature.

# There is published data in Fig 4A of [Makino and Sage (2007)](https://academic.oup.com/pcp/article/48/10/1472/1938331),
# which plots the change in activation state across temperatures for rice: Data/Makino_Sage_Fig_4-01.png.

# To be able to extract this data, the points can be estimated manually or using measuring tools or software such as:
# [GraphReader](https://www.graphreader.com/)
# [PlotDigitizer](https://plotdigitizer.com/app) 
# [WebPlotDigitizer](https://automeris.io/).

# The points can then be saved in a .csv file to be loaded later (`Fig_4A-plot-data.csv`).

# Define the data points:

# Use manual estimates
# x <-c(16,20,25,32,37,41) # Temperature (°C)
# y <-c(85,85,80,73,69,48) # Activation state (%)

# Alternatively, you can load data from the .csv file compiled using one of the automated tools: 

# These points are extracted using Web Plot Digitizer:

plot_data <- read.csv("Data/Fig_4A_plot_data.csv", 
                       header = FALSE)

# Convert activation state from percentage to decimal fraction
plot_data$V2<- plot_data$V2/100

# Specify variables:

x <- plot_data$V1 # Temperature (°C)
y <- plot_data$V2 # Activation state (proportion)

print(plot_data)

# Plot the points, setting y-axis breaks at intervals of 0.2:

ggplot(plot_data, 
       aes(x, y)) +
       geom_line(color = "black", 
                 linewidth = 1) +
       geom_point(color = "black", 
                  fill = "black", 
                  shape = 21, 
                  size = 2) +
       scale_x_continuous(limits = c(10, 45)) +
       scale_y_continuous(limits = c(0, 1.2), 
                          breaks = seq(0.2, 1.2, by = 0.2)) +
       labs(x = "Leaf temperature (°C)", 
            y = "Rubisco activation") +
       theme_minimal()


# Polynomial curve fitting
 
# Try fitting the data to quadratic, cubic and quartic models:

poly_model2 <- lm(y ~ poly(x, 
                           2, 
                           raw = TRUE),
                      data = plot_data)
poly_model3 <- lm(y ~ poly(x, 
                           3, 
                           raw = TRUE),
                      data = plot_data)
poly_model4 <- lm(y ~ poly(x, 
                           4, 
                           raw = TRUE),
                      data = plot_data)

# Evaluation of fitted models

# Check the summary statistics for the three fits:

summary_quad <- summary(poly_model2)
summary_quad

summary_cubic <- summary(poly_model3)
summary_cubic

summary_quartic <- summary(poly_model4)
summary_quartic

# Extract the fit coefficients for the quadratic and cubic models:

coef_quad <- coef(poly_model2)
coef_quad
coef_cubic <- coef(poly_model3)
coef_cubic

# Extract the residuals for the quadratic and cubic models:

residuals_quad <- residuals(poly_model2)
residuals_cubic <- residuals(poly_model3)

# Save the polynomial functions including the model coefficients for the quadratic and cubic models:

quad_func <- function(x) coef_quad[1] + 
                          coef_quad[2] * x + 
                          coef_quad[3] * x^2
cubic_func <- function(x) coef_cubic[1] + 
                          coef_cubic[2] * x + 
                          coef_cubic[3] * x^2 + 
                          coef_cubic[4] * x^3

# Obtain formulae for the polynomial functions:

quad_formula <- bquote(Quadratic: y == 
                       .(round(coef_quad[3], 4)) * x^2 + 
                       .(round(coef_quad[2], 3)) * x + 
                       .(round(coef_quad[1], 2)))

cubic_formula <- bquote(Cubic: y == 
                        .(round(coef_cubic[4], 5)) * x^3 + 
                        .(round(coef_cubic[3], 4)) * x^2 + 
                        .(round(coef_cubic[2], 3)) * x + 
                        .(round(coef_cubic[1], 2)))
quad_formula
cubic_formula

# The `bquote` function ensures creation of the formulae with LaTeX formatting.
#' `.()` is used to evaluate the expressions and insert the coefficients dynamically into the formulae.

# Plot the polynomial models on top of the data points using stat_function:

ggplot(plot_data, 
       aes(x, y)) +
  geom_line(color = "black", 
            linewidth = 1) +
  geom_point(color = "black", 
             fill = "black", 
             shape = 21, 
             size = 2) +
  stat_function(fun = quad_func, 
                aes(color = "Quadratic"), 
                linewidth = 1, 
                linetype = "dashed",
                xlim = c(0, 45)) +  
  stat_function(fun = cubic_func, 
                aes(color = "Cubic"), 
                linewidth = 1, 
                linetype = "dashed",
                xlim = c(0, 45)) +  
  scale_x_continuous(limits = c(10, 45)) +
  scale_y_continuous(limits = c(0, 1.2), 
                     breaks = seq(0.2, 1.2, by = 0.2)) +
  labs(x = "Leaf temperature (°C)", 
       y = "Rubisco activation",
       color = "Polynomial Fit") +
  theme_minimal() +
  scale_color_manual(
    values = c("Quadratic" = "red", 
               "Cubic" = "blue"),
    labels = c(quad_formula, 
               cubic_formula))

# Obtain predictions using the polynomial models:

y_pred_quad <- predict(poly_model2)
y_pred_cubic <- predict(poly_model3)

# Plot the residuals against temperature for the quadratic model:

ggplot(plot_data, 
       aes(x = x)) +
       geom_point(aes
                  (y = residuals_quad), 
                  color = "red") +
       geom_hline(yintercept = 0,
                  linetype = "dashed") +
       ggtitle("Residuals of Quadratic Model vs. Temperature") +
       xlab("Leaf temperature (°C)") +
       ylab("Residuals") +
       theme_minimal()

# Plot the residuals against temperature for the cubic model:

ggplot(plot_data, 
       aes(x = x)) +
       geom_point(aes
                  (y = residuals_cubic), 
                  color = "blue") +
       geom_hline(yintercept = 0,
                  linetype = "dashed") +
       ggtitle("Residuals of Cubic Model vs. Temperature") +
       xlab("Leaf temperature (°C)") +
       ylab("Residuals") +
       theme_minimal()

# Plot the residuals against fitted values for the quadratic model:

ggplot(plot_data, 
       aes(x = y_pred_quad)) +
       geom_point(aes
                  (y = residuals_quad), 
                  color = "red") +
       geom_hline(yintercept = 0,
                  linetype = "dashed") +
       ggtitle("Residuals of Quadratic Model vs. Fitted Values") +
       xlab("Fitted Values") +
       ylab("Residuals") +
       theme_minimal()

# Plot the residuals against fitted values for the cubic model:

ggplot(plot_data, 
       aes(x = y_pred_cubic)) +
       geom_point(aes
                  (y = residuals_cubic), 
                  color = "blue") +
       geom_hline(yintercept = 0,
                  linetype = "dashed") +
       ggtitle("Residuals of Cubic Model vs. Fitted Values") +
       xlab("Fitted Values") +
       ylab("Residuals") +
       theme_minimal()

# Compute MSE, RMSE and AIC for the models:

mse_quad <- mean((y - y_pred_quad)^2)
mse_quad
mse_cubic <- mean((y - y_pred_cubic)^2)
mse_cubic
rmse_quad <- sqrt(mean((y - y_pred_quad)^2))
rmse_quad
rmse_cubic <- sqrt(mean((y - y_pred_cubic)^2))
rmse_cubic
aic_quad <- AIC(poly_model2)
aic_quad
aic_cubic <- AIC(poly_model3)
aic_cubic

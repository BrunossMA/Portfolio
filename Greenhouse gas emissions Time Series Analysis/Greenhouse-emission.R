# Load necessary libraries
library(fpp2)
library(imputeTS)
library(TSA)
library(forecast)
library(tidyverse)
library(lubridate)
library(tseries)
library(qpcR)
library(ICglm)

# Source additional functions
source("My_functions.R")
set.seed(224)

# Read data from the CSV file
dat <- read.csv("C:/Users/mauro/OneDrive/Desktop/Project Time Series/Wolrd/per-capita-ghg-emissions.csv", header = TRUE)

# Filter data for France and reset row names
data <- dat[dat[, 1] == "France", ]
row.names(data) <- NULL

# Select relevant columns and rename them
data <- data[, c(3, 4)]
colnames(data) <- c("date", "x")

# Convert the emissions data to numeric
data$x <- as.numeric(data$x)

# Function to create a standardized mean-variance plot
stdmean_plot <- function(dat, lambda = "No", blocks = 10, title = "Original Data") {
  if (lambda != "No") {
    dat <- BoxCox(dat, lambda = lambda)
  }
  
  means <- rep(0, blocks)
  stds <- rep(0, blocks)
  
  idx <- seq(0, 205, length(dat) %/% blocks)
  idx[1] <- 1
  
  for (i in 2:length(idx)) {
    means[i - 1] <- mean(dat[idx[i - 1]:idx[i] - 1])
    stds[i - 1] <- sd(dat[idx[i - 1]:idx[i] - 1])
  }
  
  plot(means, stds, ylim = c(0.00001, max(dat) / 2), main = title)
}

# Filter data to include only records from 1950 onward
data <- data[data$date > 1950, ]

# Define the time series object
ts <- ts(data$x, frequency = 1, start = c(1950))
autoplot(ts, col = "red") +
  xlab('Date') +
  ggtitle('Per Capita Greenhouse Gas Emissions') +
  ylab("Tonnes")

# Plot moving averages
autoplot(ts, series = "Data") +
  autolayer(ma(ts, 5), series = "5-MA") +
  autolayer(ma(ts, 10), series = "10-MA") +
  autolayer(ma(ts, 15), series = "15-MA") +
  xlab("Year") + 
  ylab("Tonnes") +
  ggtitle("Per Capita Greenhouse Gas Emissions") +
  scale_colour_manual(values = c("Data" = "grey50", "5-MA" = "red", "10-MA" = "blue", "15-MA" = "green"),
                      breaks = c("Data", "5-MA", "10-MA", "15-MA"))

# Residual calculation and analysis
A <- ma(ts, order = 5)  # Example of a moving average for comparison
RES <- ts - A
checkresiduals(RES)

# Lowess smoothing
autoplot(ts, series = "Data") +
  autolayer(A, series = "Lowess") +
  xlab("Year") + 
  ylab("Tonnes") +
  ggtitle("Per Capita Greenhouse Gas Emissions") +
  scale_colour_manual(values = c("Data" = "grey50", "Lowess" = "red"),
                      breaks = c("Data", "Lowess"))

# Standard deviation-mean plots
par(mfrow = c(1, 2))
stdmean_plot(ts, "No", 10)
stdmean_plot(ts, BoxCox.lambda(ts), 6)
stdmean_plot(ts, 0, 6)
stdmean_plot(ts, BoxCox.lambda(ts, method = "guerrero"), 10, "Transformed Data")

# Box-Cox transformation
ts_box <- BoxCox(ts, BoxCox.lambda(ts, method = "guerrero"))
autoplot(ts_box, col = "blue") +
  xlab('Date') + 
  ylab('Tonnes') + 
  ggtitle("Per Capita Greenhouse Gas Emissions (Box-Cox Transformed)") +
  theme_get()

# Autocorrelation plots
ggAcf(ts) + ggtitle('Autocorrelation Function')
ggPacf(ts) + ggtitle('Partial Autocorrelation Function')
ggAcf(ts_box) + ggtitle('Autocorrelation Function (Transformed)')
ggPacf(ts_box) + ggtitle('Partial Autocorrelation Function (Transformed)')

# Differencing the series
dxt <- diff(ts_box, 1)
dxt2 <- diff(dxt, 1)

# Plotting differenced series
gridExtra::grid.arrange(
  autoplot(ts_box) + ggtitle(expression("Time Series w"[t])) + ylab(""),
  autoplot(dxt) + ggtitle(expression(Delta * w[t])) + ylab(""), 
  autoplot(dxt2) + ggtitle(expression(Delta^2 * w[t])) + ylab(""), ncol = 2
)

# ACF and PACF of differenced series
gridExtra::grid.arrange(
  ggAcf(dxt) + ggtitle(expression(Delta * w[t])),
  ggAcf(dxt2) + ggtitle(expression(Delta^2 * w[t])), ncol = 2
)

gridExtra::grid.arrange(
  ggPacf(dxt) + ggtitle(expression(Delta * w[t])),
  ggPacf(dxt2) + ggtitle(expression(Delta^2 * w[t])), ncol = 2
)

# Variance analysis
var(ts_box)
var(dxt) / var(ts_box)
var(dxt2) / var(ts_box)

# Augmented Dickey-Fuller test for stationarity
adf.test(dxt)
adf.test(dxt2)

# ACF and PACF for second differenced series
gridExtra::grid.arrange(
  ggAcf(dxt2, lag.max = 24) + ggtitle(expression("ACF " * Delta^2 * "w"[t])),
  ggPacf(dxt2, lag.max = 24) + ggtitle(expression("PACF "*Delta^2 * "w"[t])), ncol = 1
)

# Forecasting
train <- window(ts_box, start = c(1950), end = c(2010))
test <- window(ts_box, start = c(2011))

# Model selection using BIC and AIC
all_fit <- expand.grid(p = 0:3, d = 2, q = 0:3, SBIC = 0, HQIC = 0, LOG = 0, V = 0)
n <- 0
for (i in 1:nrow(all_fit)) {
  skip_to_next <- FALSE
  tryCatch(
    {
      fit <- Arima(train, order = c(all_fit[i, 1], 2, all_fit[i, 3]))
      all_fit[i, 4] <- BICadj(fit)
      all_fit[i, 5] <- HQIC(fit)
      all_fit[i, 6] <- fit$loglik
      all_fit[i, 7] <- fit$sigma2
      n <- n + 1
      print(n)
    }, error = function(e) { skip_to_next <<- TRUE }
  )
  
  if (skip_to_next) {
    next
  }
}

# Model selection based on BIC and AIC
selection <- all_fit
selection[order(selection$SBIC, decreasing = FALSE), ]
selection[order(selection$HQIC, decreasing = FALSE), ]

# Fit selected models
fit_t1 <- Arima(train, order = c(0, 2, 1))
fit_t2 <- Arima(train, order = c(1, 2, 1))
fit_t3 <- Arima(train, order = c(3, 2, 2))

# Residual analysis for selected models
gridExtra::grid.arrange(
  ggAcf(fit_t1$residuals, lag.max = 12) + ggtitle("ACF ARIMA(0,2,1)"),
  ggPacf(fit_t1$residuals, lag.max = 12) + ggtitle("PACF ARIMA(0,2,1)"),
  ggAcf(fit_t2$residuals, lag.max = 12) + ggtitle("ACF ARIMA(1,2,1)"),
  ggPacf(fit_t2$residuals, lag.max = 12) + ggtitle("PACF ARIMA(1,2,1)"),
  ggAcf(fit_t3$residuals, lag.max = 12) + ggtitle("ACF ARIMA(3,2,2)"),
  ggPacf(fit_t3$residuals, lag.max = 12) + ggtitle("PACF ARIMA(3,2,2)"), ncol = 2
)

# Residual analysis plots
gridExtra::grid.arrange(
  autoplot(fit_t1$residuals) + ylab("Residuals") + ggtitle("ARIMA(0,2,1)"),
  autoplot((fit_t1$residuals)^2) + ylab("Squared Residuals") + ggtitle("ARIMA(0,2,1)"), ncol = 2
)

# QQ plot for normality check
qqnorm(fit_t1$residuals)
grid()  # Add grid to the plot
qqline(fit_t1$residuals, lwd = 2, col = "red")  # Add a QQ line

# Check residuals and perform Jarque-Bera test
checkresiduals(fit_t3)
jarque.bera.test(fit_t2$residuals)

# Coefficient test for the fitted model
coeftest(fit_t3)

# Make predictions
pred <- forecast(fit_t1, h = 10)
pred

# Calculate accuracy metrics
RMSE <- sqrt(mean((test - pred$mean)^2))
MAE <- mean(abs(pred$mean - test))
MPE <- 100 * mean((test - pred$mean) / test)
MAPE <- 100 * mean(abs(pred$mean - test) / abs(test))

# Inverse Box-Cox transformation for predictions
pred_x <- InvBoxCox(pred$mean, BoxCox.lambda(ts, method = "guerrero"))
pred_upper <- InvBoxCox(pred$upper, BoxCox.lambda(ts, method = "guerrero"))
pred_lower <- InvBoxCox(pred$lower, BoxCox.lambda(ts, method = "guerrero"))

# Plot predictions
autoplot(InvBoxCox(pred$mean, BoxCox.lambda(ts, method = "guerrero")), series = "Prediction") +
  autolayer(ts, colour = TRUE, series = "Time Series") +
  geom_ribbon(aes(ymin = InvBoxCox(pred$lower[, 2], BoxCox.lambda(ts, method = "guerrero")), 
                  ymax = InvBoxCox(pred$upper[, 2], BoxCox.lambda(ts, method = "guerrero"))), 
              alpha = 0.2, fill = "green") +
  ggtitle('Prediction with an ARIMA(1,1,1)(0,1,1)[12]') + 
  ylab("Tonnes")

# Subjective model selection based on ACF-PACF plots
fit_g <- Arima(train, order = c(0, 1, 0), include.constant = TRUE)
autoplot(fit_g$residuals)

# QQ plot for model residuals
qqnorm(fit_g$residuals)
grid()  # Add grid to the plot
qqline(fit_g$residuals, lwd = 2, col = "red")  # Add a QQ line

# Check residuals and perform Jarque-Bera test
checkresiduals(fit_g$residuals)
jarque.bera.test(fit_g$residuals)

# Coefficient test for the fitted model
coeftest(fit_g)

# Make predictions with the subjective model
pred_g <- forecast(fit_g, h = 10)
pred_g

# Calculate accuracy metrics for the subjective model
RMSE <- sqrt(mean((test - pred_g$mean)^2))
MAE <- mean(abs(pred_g$mean - test))
MPE <- 100 * mean((test - pred_g$mean) / test)
MAPE <- 100 * mean(abs(pred_g$mean - test) / test)

# Inverse Box-Cox transformation for predictions from the subjective model
pred_x <- InvBoxCox(pred_g$x, BoxCox.lambda(ts, method = "guerrero"))
pred_upper <- InvBoxCox(pred_g$upper, BoxCox.lambda(ts, method = "guerrero"))
pred_lower <- InvBoxCox(pred_g$lower, BoxCox.lambda(ts, method = "guerrero"))

# Plot predictions for the subjective model
autoplot(InvBoxCox(pred_g$mean, BoxCox.lambda(ts, method = "guerrero")), series = "Prediction") +
  autolayer(ts, colour = TRUE, series = "Time Series") +
  geom_ribbon(aes(ymin = InvBoxCox(pred_g$lower[, 2], BoxCox.lambda(ts, method = "guerrero")), 
                  ymax = InvBoxCox(pred_g$upper[, 2], BoxCox.lambda(ts, method = "guerrero"))), 
              alpha = 0.2, fill = "green") +
  ggtitle('Prediction with an ARIMA(1,1,1)(0,1,1)[12]') + 
  ylab("Tonnes")

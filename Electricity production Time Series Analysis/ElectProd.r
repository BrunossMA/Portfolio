# Load necessary libraries
library(fpp2)            # For forecasting and time series analysis
library(forecast)        # For ARIMA modeling and forecasting
library(tidyverse)       # For data manipulation and visualization
library(lubridate)       # For date manipulation
library(gridExtra)       # For arranging multiple grid layouts
library(kableExtra)      # For creating formatted tables

# Set seed for reproducibility
set.seed(224)

# Read the CSV file
dat <- read.csv("/Electricity production Time Series Analysis/MES062023.csv", header = TRUE)

# Function to generate a standardized mean plot
stdmean_plot <- function(dat, lambda = "No", blocks = 10, title) {
  if (lambda != "No") {
    dat <- BoxCox(dat, lambda = lambda)
  }
  
  means <- numeric(blocks)
  stds <- numeric(blocks)
  idx <- seq(0, 205, length.out = length(dat) %/% blocks)
  idx[1] <- 1
  
  for (i in 2:length(idx)) {
    means[i - 1] <- mean(dat[idx[i - 1]:idx[i] - 1])
    stds[i - 1] <- sd(dat[idx[i - 1]:idx[i] - 1])
  }
  
  plot(means, stds, ylim = c(0, max(dat) / 2), main = title)
}

# Function to extract data based on country and type
extract <- function(country, type) {
  data <- dat[dat[, 1] == country & str_detect(dat[, 4], type), -c(1, 3, 4, 6)]
  row.names(data) <- NULL
  colnames(data) <- c("date", type)
  data$date <- dmy(paste("01", data$date))  # Convert date format
  data[, 2] <- rev(as.numeric(data[, 2]))   # Reverse the order of values
  return(data)
}

# Extract data for France
data <- extract("France", "Total Renewables")
data$`Total Renewables` <- as.numeric(data$`Total Renewables`)

# Define time series
ts <- ts(data[, 2], frequency = 12, start = c(2010, 1), end = c(2023, 6))
autoplot(ts, col = "blue") +
  xlab('Date') +
  ggtitle('French Electricity Production from Renewable Sources') +
  theme_get()

# Seasonal plots
ggseasonplot(ts, polar = FALSE) + 
  ggtitle('Seasonal Plot of Electricity Production') + 
  ylab('GWh')

ggsubseriesplot(ts) + 
  theme_bw() + 
  ylab("GWh") + 
  ggtitle("Sub-series Plot of Electricity Production")

# Box-Cox transformation analysis
par(mfrow = c(1, 2))

stdmean_plot(ts, "No", 10, "Original Data")
stdmean_plot(ts, BoxCox.lambda(ts), 10)

stdmean_plot(ts, BoxCox.lambda(ts, method = "guerrero"), 10, "Transformed Data")
stdmean_plot(ts, BoxCox.lambda(ts, method = "loglik"), 10)

ts_box <- BoxCox(ts, BoxCox.lambda(ts, method = "guerrero"))

autoplot(ts_box, col = "red") +
  xlab('Date') + 
  ylab('GWh') + 
  ggtitle('French Electricity Production from Renewable Sources (Box-Cox Transformed)') +
  theme_get()

# Autocorrelation functions
ggAcf(ts) + ggtitle('Autocorrelation Function')
ggPacf(ts) + ggtitle('Partial Autocorrelation Function')
ggAcf(ts_box) + ggtitle('Autocorrelation Function (Box-Cox)')
ggPacf(ts_box) + ggtitle('Partial Autocorrelation Function (Box-Cox)')

# Additive decomposition
additive <- decompose(ts, type = "additive")
autoplot(additive)

# Check residuals
checkresiduals(remainder(additive))

# Normality test
shapiro.test(remainder(additive))

# Visualization of decomposition components
autoplot(ts, series = "Data") +
  autolayer(seasadj(additive), series = "Seasonally Adjusted Series", size = 1) +
  autolayer(trendcycle(additive), series = "Trend Cycle", size = 1) +
  theme(legend.position = "bottom")

# Multiplicative decomposition
multiplicative <- decompose(ts, type = "multiplicative")
autoplot(multiplicative)

checkresiduals(remainder(multiplicative))
shapiro.test(remainder(multiplicative))

# STL decomposition
stl_box <- mstl(ts, lambda = NULL)
autoplot(stl_box)
checkresiduals(remainder(stl_box))
shapiro.test(remainder(stl_box))

# Differencing to remove trend and seasonality
dxt <- diff(ts_box, 1)
dd12xt <- diff(dxt, 12)
d12xt <- diff(ts_box, 12)

# Plot differences
gridExtra::grid.arrange(
  autoplot(ts_box) + ylab("") + ggtitle(expression("Time Series w"[t])),
  autoplot(dxt) + ylab("") + ggtitle(expression(Delta * "w"[t])),
  autoplot(d12xt) + ylab("") + ggtitle(expression(Delta[12] * "w"[t])),
  autoplot(dd12xt) + ylab("") + ggtitle(expression(Delta * Delta[12] * "w"[t])), ncol = 2
)

# ACF and PACF of differences
gridExtra::grid.arrange(
  ggAcf(ts_box) + ylab("") + ggtitle(expression("Time Series w"[t])),
  ggAcf(dxt) + ylab("") + ggtitle(expression(Delta * "w"[t])),
  ggAcf(d12xt) + ylab("") + ggtitle(expression(Delta[12] * "w"[t])),
  ggAcf(dd12xt) + ylab("") + ggtitle(expression(Delta * Delta[12] * "w"[t])), ncol = 2
)

gridExtra::grid.arrange(
  ggPacf(ts_box) + ggtitle("Original Series"),
  ggPacf(dxt) + ggtitle("dXt"),
  ggPacf(d12xt) + ggtitle("D12Xt"),
  ggPacf(dd12xt) + ggtitle("dD12Xt"), ncol = 2
)

# Variance analysis for differencing
var_ts_box <- var(ts_box)
var_dxt <- var(dxt) / var_ts_box
var_d12xt <- var(d12xt) / var_ts_box
var_dd12xt <- var(dd12xt) / var_ts_box

dat_variance <- data.frame(
  Variance = c(var(ts_box), var(dxt), var(d12xt), var(dd12xt)),
  Fraction_Variance = c("-", var_dxt, var_d12xt, var_dd12xt)
)
rownames(dat_variance) <- c(
  expression("Time Series w"[t]),
  expression(Delta * "w"[t]),
  expression(Delta[12] * "w"[t]),
  expression(Delta * Delta[12] * "w"[t])
)

# Display variance metrics
kableExtra::kable(dat_variance)

# ACF and PACF plots for differencing
gridExtra::grid.arrange(
  ggAcf(dxt, lag.max = 24) + ggtitle(expression("ACF " * Delta * "w"[t])),
  ggPacf(dxt, lag.max = 24) + ggtitle(expression("PACF " * Delta * "w"[t])), ncol = 1
)

# Augmented Dickey-Fuller test for stationarity
adf.test(dxt, k = 24)

# Training and testing data split
train <- window(ts_box, start = c(2010, 1), end = c(2020, 12))
test <- window(ts_box, start = c(2021, 1))

# Model selection for ARIMA
all_fit <- expand.grid(
  ar = 0:3, d = 1, ma = 0:3, AR = 0:3, D = 0, MA = 0:3,
  BICadj = 0, HQIC = 0, LOG = 0, V = 0
)
n <- 0

for (i in 1:nrow(all_fit)) {
  skip_to_next <- FALSE
  tryCatch(
    {
      fit <- Arima(
        train,
        order = c(all_fit[i, 1], 1, all_fit[i, 3]),
        seasonal = list(order = c(all_fit[i, 4], 0, all_fit[i, 6]), period = 12)
      )
      all_fit[i, 7] <- BICadj(fit)
      all_fit[i, 8] <- HQIC(fit)
      all_fit[i, 9] <- logLik(fit)
      all_fit[i, 10] <- var(residuals(fit))
      n <- n + 1
    },
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  if (skip_to_next) {
    next
  }
}

# Results
all_fit <- all_fit[all_fit$BICadj != 0, ]
all_fit %>%
  arrange(BICadj) %>%
  slice(1:5)

# Selection and evaluation of the best model
best_fit <- Arima(train, order = c(3, 1, 2), seasonal = list(order = c(2, 0, 3), period = 12))
summary(best_fit)
checkresiduals(best_fit)

# Residual analysis
autoplot(best_fit$residuals) +
  ggtitle("Residuals of the Best Model") +
  ylab("Residuals")

qqnorm(best_fit$residuals)
grid()  # Add grid
qqline(best_fit$residuals, lwd = 2, col = "red")  # Reference line

checkresiduals(best_fit)
jarque.bera.test(best_fit$residuals)

coeftest(best_fit)

# Forecasting
pred <- forecast(best_fit, h = 30)
autoplot(pred) +
  autolayer(test, series = "Test Set") +
  ggtitle("Forecast of Electricity Production in France") +
  xlab("Date") +
  ylab("GWh")

# Evaluate forecasting accuracy
RMSE <- sqrt(mean((test - pred$mean)^2))
MAE <- mean(abs(pred$mean - test))
MPE <- 100 * mean((test - pred$mean) / test)
MAPE <- 100 * mean(abs(pred$mean - test) / test)

# Display error metrics
error_metrics <- data.frame(RMSE, MAE, MPE, MAPE)
kableExtra::kable(error_metrics)

# Inverse Box-Cox transformation for forecasts
pred_x <- InvBoxCox(pred$x, BoxCox.lambda(ts, method = "guerrero"))
pred_upper <- InvBoxCox(pred$upper, BoxCox.lambda(ts, method = "guerrero"))
pred_lower <- InvBoxCox(pred$lower, BoxCox.lambda(ts, method = "guerrero"))

autoplot(InvBoxCox(pred$mean, BoxCox.lambda(ts, method = "guerrero")), series = "Prediction") +
  autolayer(ts, colour = TRUE, series = "Time Series") +
  geom_ribbon(aes(
    ymin = InvBoxCox(pred$lower[, 2], BoxCox.lambda(ts, method = "guerrero")),
    ymax = InvBoxCox(pred$upper[, 2], BoxCox.lambda(ts, method = "guerrero"))
  ), alpha = 0.2, fill = "green") +
  ggtitle('Forecast with ARIMA(3,1,2)(2,0,3)[12]') +
  ylab("GWh")

# Repeat analysis with a second ARIMA model
fit_g <- Arima(train, order = c(3, 1, 1), seasonal = list(order = c(1, 1, 1), period = 12))
autoplot(fit_g$residuals) +
  ggtitle("Residuals of the Second Model") +
  ylab("Residuals")

qqnorm(fit_g$residuals)
grid()  # Add grid
qqline(fit_g$residuals, lwd = 2, col = "red")

checkresiduals(fit_g$residuals)
jarque.bera.test(fit_g$residuals)
BICadj(fit_g)
HQIC(fit_g)

coeftest(fit_g)

# Forecasting with the second model
pred_g <- forecast(fit_g, h = 30)
autoplot(pred_g) +
  autolayer(test, series = "Test Set") +
  ggtitle("Forecast of Electricity Production in France with Second Model") +
  xlab("Date") +
  ylab("GWh")

# Calculate error metrics for the second model
RMSE_g <- sqrt(mean((test - pred_g$mean)^2))
MAE_g <- mean(abs(pred_g$mean - test))
MPE_g <- 100 * mean((test - pred_g$mean) / test)
MAPE_g <- 100 * mean(abs(pred_g$mean - test) / test)

# Display error metrics for the second model
error_metrics_g <- data.frame(RMSE_g, MAE_g, MPE_g, MAPE_g)
kableExtra::kable(error_metrics_g)

# Inverse Box-Cox transformation for the second model forecasts
pred_g_x <- InvBoxCox(pred_g$x, BoxCox.lambda(ts, method = "guerrero"))
pred_g_upper <- InvBoxCox(pred_g$upper, BoxCox.lambda(ts, method = "guerrero"))
pred_g_lower <- InvBoxCox(pred_g$lower, BoxCox.lambda(ts, method = "guerrero"))

autoplot(InvBoxCox(pred_g$mean, BoxCox.lambda(ts, method = "guerrero")), series = "Prediction") +
  autolayer(ts, colour = TRUE, series = "Time Series") +
  geom_ribbon(aes(
    ymin = InvBoxCox(pred_g$lower[, 2], BoxCox.lambda(ts, method = "guerrero")),
    ymax = InvBoxCox(pred_g$upper[, 2], BoxCox.lambda(ts, method = "guerrero"))
  ), alpha = 0.2, fill = "green") +
  ggtitle('Forecast with ARIMA(3,1,1)(1,1,1)[12]') +
  ylab("GWh")

# Clean up unnecessary variables
rm(dat, data, dxt, dd12xt, ts_box)

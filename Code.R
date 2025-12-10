# MATH 271.2 Advanced Time Series and Forecasting
# Group 5 (Cheng, De Leon, Montemayor)
# Project 1 (Linear Time Series Models)

# A. Importing Libraries and Datasets

# Time Series Libraries
# library(TSA)
library(tseries)
# library(itsmr)
library(forecast)
# library(astsa)

# Data Visualization Libraries
library(ggplot2)
library(patchwork)

df <- read.csv(file.choose()) # find disease_pidsr_totals
df

# B. Data Preprocessing 

# 1. Data Filtering
disease_location_filter <- function(df_name, disease_name, location_name) {
  df_filtered <- df_name[df_name$disease_common_name == disease_name & 
                           df_name$adm3_pcode == location_name, names(df)]
  df_filtered <- aggregate(case_total ~ date, data = df_filtered, sum)
  df_filtered <- df_filtered[order(df_filtered$date), ]
  ts <- ts(df_filtered$case_total, frequency = 52, start = c(2008, 1))
  return(ts)
}

# Location codes: Refer to location.csv
# PH015518000: Dagupan City, Pangasinan (DAG) 
# PH034919000: Palayan City, Nueva Ecija (PAL)
# PH050506000: Legazpi City, Albay (LEG)
# PH063022000: Iloilo City, Iloilo (ILO)
# PH072230000: Mandaue City, Cebu (MDE)
# PH083747000: Tacloban City, Leyte (TAC) 
# PH097332000: Zamboanga City, Zamboanga del Sur (ZAM)
# PH104305000: Cagayan de Oro City, Misamis Oriental (CAG) 
# PH112402000: Davao City, Davao del Sur (DAV)
# PH137401000: City of Mandaluyong, NCR, Second District (MLG)
# PH137503000: City of Navotas, NCR, Third District (NAV)
# PH137603000: City of Muntinlupa, NCR, Fourth District (MUN)

dengue_dag <- disease_location_filter(df, "DENGUE FEVER", "PH015518000")
dengue_pal <- disease_location_filter(df, "DENGUE FEVER", "PH034919000")
dengue_leg <- disease_location_filter(df, "DENGUE FEVER", "PH050506000")
dengue_ilo <- disease_location_filter(df, "DENGUE FEVER", "PH063022000")
dengue_mde <- disease_location_filter(df, "DENGUE FEVER", "PH072230000")
dengue_tac <- disease_location_filter(df, "DENGUE FEVER", "PH083747000")
dengue_zam <- disease_location_filter(df, "DENGUE FEVER", "PH097332000")
dengue_cag <- disease_location_filter(df, "DENGUE FEVER", "PH104305000")
dengue_dav <- disease_location_filter(df, "DENGUE FEVER", "PH112402000")
dengue_mlg <- disease_location_filter(df, "DENGUE FEVER", "PH137401000")
dengue_nav <- disease_location_filter(df, "DENGUE FEVER", "PH137503000")
dengue_mun <- disease_location_filter(df, "DENGUE FEVER", "PH137603000")

dengue_all <- cbind(time(dengue_dag), dengue_dag, dengue_pal, dengue_leg, dengue_ilo, dengue_mde, dengue_tac,
                    dengue_zam, dengue_cag, dengue_dav, dengue_mlg, dengue_nav, dengue_mun)
colnames(dengue_all) <- c("Time", "Dagupan", "Palayan", "Legazpi", "Iloilo", "Mandaue", "Tacloban",
                   "Zamboanga", "Cagayan", "Davao", "Mandaluyong", "Navotas", "Muntinlupa")
dengue_all <- as.data.frame(dengue_all)

# 2. Plotting dengue cases in all regions [Figure 1]
plots <- lapply(1:12, function(i) {
  ggplot(dengue_all, aes(x = .data[["Time"]], 
                         y = .data[[colnames(dengue_all)[i + 1]]])
         ) +
    geom_line() +
    labs(y = colnames(dengue_all)[[i+1]]) + # add y-axis label
    theme_minimal() +
    theme( # remove x-axis ticks and label
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
})

for (i in 9:12) {
  plots[[i]] <- plots[[i]] +
    theme( # add x-axis ticks and label for bottom plots
      axis.title.x = element_text(),
      axis.text.x = element_text()
    ) +
    labs(x = "Time")
 
}

# Arrange plots
dengue_all_plot <- (plots[[1]] / plots[[5]] / plots[[9]]) | 
  (plots[[2]] / plots[[6]] / plots[[10]]) | 
  (plots[[3]] / plots[[7]] / plots[[11]]) | 
  (plots[[4]] / plots[[8]] / plots[[12]])
dengue_all_plot <- dengue_all_plot + plot_annotation(
  title = "Dengue Cases in Filipino Cities",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
)
dengue_all_plot

# 3. Plotting dengue cases in Iloilo [Figure 2]
# Convert time series to data frame with a time index
df_ilo <- data.frame(
  Time = as.numeric(time(dengue_ilo)),
  Cases = as.numeric(dengue_ilo)
)

dengue_ilo_plot <- ggplot(df_ilo, aes(x = Time, y = Cases)) +
  geom_line() +
  scale_x_continuous(breaks = seq(2008, 2022, by = 1)) +  # yearly ticks
  labs(x = "Time", y = "Cases") +
  theme_minimal() +
  plot_annotation(
    title = "Dengue Cases in Iloilo",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
dengue_ilo_plot

# C. Hypothesis Testing 
# 1. No transformations
adf.test(dengue_ilo)                 # stationary
Box.test(dengue_ilo, type="Ljung") # autocorrelated
# However, the plot evidently has seasonality, which we can estimate to be annual (52 weeks).

# 2. log
dengue_ilo_log = log(dengue_ilo+1)  # +1 to avoid log(0)
plot(dengue_ilo_log)
adf.test(dengue_ilo_log)               # stationary
Box.test(dengue_ilo_log, type="Ljung") # autocorrelated

# 3. log, s=52, D=1
dengue_ilo_log_s52_D1 = diff(log(dengue_ilo+1), 52)
plot(dengue_ilo_log_s52_D1)
adf.test(dengue_ilo_log_s52_D1)               # not stationary
Box.test(dengue_ilo_log_s52_D1, type="Ljung") # autocorrelated

# 4. log, s=52, D=1, d=1
dengue_ilo_log_s52_D1_d1 = diff(diff(log(dengue_ilo+1), 52), 1)
plot(dengue_ilo_log_s52_D1_d1)
adf.test(dengue_ilo_log_s52_D1_d1)               # stationary
Box.test(dengue_ilo_log_s52_D1_d1, type="Ljung") # autocorrelated

# 5. log, s=156, D=1
dengue_ilo_s156 <- ts(dengue_ilo, frequency = 156, start = c(2008, 1))
dengue_ilo_log_s156 = log(dengue_ilo_s156+1)
dengue_ilo_log_s156_D1 = diff(dengue_ilo_log_s156, 156)
adf.test(dengue_ilo_log_s156_D1)                  # not stationary
Box.test(dengue_ilo_log_s156_D1, type = "Ljung")  # autocorrelated

# 6. log, s=156, D=1, d=1
dengue_ilo_log_s156_D1_d1 <- diff(dengue_ilo_log_s156_D1, 1)
adf.test(dengue_ilo_log_s156_D1_d1)                   # stationary
Box.test(dengue_ilo_log_s156_D1_d1, type = "Ljung")   # autocorrleated

# D. Order Determination
# 1. s=52, D=1, d=1 [Figure 3]
pacf(dengue_ilo_log_s52_D1_d1, lag = 520, main = "PACF of Time Series with One-Year Season")
acf(dengue_ilo_log_s52_D1_d1, lag = 520, main = "ACF of Time Series with One-Year Season")

# 2. s=156, D=1, d=1 [Figure 4]
pacf(dengue_ilo_log_s156_D1_d1, lag = 780, main = "PACF of Time Series with Three-Year Season") 
acf(dengue_ilo_log_s156_D1_d1, lag = 780, main = "ACF of Time Series with Three-Year Season") 

# E. Modeling
# AIC for each model is commented
# 1. s=52
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(0, 1, 0), period = 52)) # 1392.31
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 0), period = 52)) # 1263.52
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 0), period = 52)) # 1286.95
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 0), period = 52)) # 1265.1
Arima(dengue_ilo_log, order = c(2, 1, 0), seasonal = list(order = c(0, 1, 0), period = 52)) # 1263.68
Arima(dengue_ilo_log, order = c(2, 1, 1), seasonal = list(order = c(0, 1, 0), period = 52)) # 1264.64
Arima(dengue_ilo_log, order = c(3, 1, 0), seasonal = list(order = c(0, 1, 0), period = 52)) # 1264.24
Arima(dengue_ilo_log, order = c(3, 1, 1), seasonal = list(order = c(0, 1, 0), period = 52)) # 1262.51
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(0, 1, 1), period = 52)) # 1010.34
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52)) # 894.26
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 1), period = 52)) # 920.99
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52)) # 896.08
Arima(dengue_ilo_log, order = c(2, 1, 0), seasonal = list(order = c(0, 1, 1), period = 52)) # 897.09
Arima(dengue_ilo_log, order = c(2, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52)) # 896.72
Arima(dengue_ilo_log, order = c(3, 1, 0), seasonal = list(order = c(0, 1, 1), period = 52)) # 896.83
Arima(dengue_ilo_log, order = c(3, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52)) # 898.72

# 2. s=156
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(0, 1, 0), period = 156)) # 1123.30
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 0), period = 156)) # 1009.51
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 0), period = 156)) # 1041.22
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 0), period = 156)) # 1010.89
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(0, 1, 1), period = 156)) # 991.56
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 156)) # 882.96 [BEST]
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 1), period = 156)) # 913.17
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 1), period = 156)) # 883.47
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(1, 1, 0), period = 156)) # 1018.30
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(1, 1, 0), period = 156)) # 912.24
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(1, 1, 0), period = 156)) # Non-finite finite-difference value
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 0), period = 156)) # Non-finite initial vmmin value
Arima(dengue_ilo_log, order = c(0, 1, 0), seasonal = list(order = c(1, 1, 1), period = 156)) # 993.42
Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(1, 1, 1), period = 156)) # 884.92
Arima(dengue_ilo_log, order = c(1, 1, 0), seasonal = list(order = c(1, 1, 1), period = 156)) # Non-finite initial vmmin value
Arima(dengue_ilo_log, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 1), period = 156)) # Non-finite initial vmmin value

# 3. Functional form of the best model
model <- Arima(dengue_ilo_log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 156))

# F. Forecasting
# 1. Exploration
pred <- forecast(model, h=52*2)   # two years
time_vals <- time(pred$mean)      # time values of the prediction
plot(pred, ylab = "Cases", xlab = "Time", main = "Log-Transformed Dengue Cases Forecast")

# 2. Convert mean and standard error of log-transformed series to those of original series
mean_log <- pred$mean                             # mean of the log-transformed series
upper95 <- pred$upper[, "95%"]                    # upper bound of the 95% confidence interval
se_log <- (upper95 - mean_log) / qnorm(0.975)     # standard error of the log-transformed series (95% confidence interval)
mean_orig <- exp(mean_log + 0.5 * se_log^2)       # mean of the original series
se_orig <- mean_orig * sqrt(exp(se_log^2) - 1)    # standard error of the original series
lower_bound <- mean_orig - qnorm(0.995) * se_orig # lower bound of the 99% confidence interval of the original series 
upper_bound <- mean_orig + qnorm(0.995) * se_orig # upper bound of the 99% confidence interval of the original series 

# 3. QQ plot for normality [Figure 5]
qqline(model$residuals, col="red")

# 4. Plot of past and predicted data (mean) [Figure 6]
complete_ts <- round(ts(c(dengue_ilo, mean_orig), start = start(dengue_ilo), frequency = 52))
df_forecast_with_historical <- data.frame(
  Time = time(complete_ts),
  Cases = as.numeric(complete_ts),
  Source = c(rep("Past", length(dengue_ilo)), 
             rep("Forecast", length(complete_ts) - length(dengue_ilo)))
)

ggplot(df_forecast_with_historical, aes(x = Time, y = Cases, color = Source)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Past" = "black", "Forecast" = "#1874CD")) +
  scale_x_continuous(breaks = seq(2008, 2026, by = 1)) +  # yearly ticks
  labs(x = "Time", y = "Cases") +
  theme_minimal() +
  theme(legend.position = "none") +
  plot_annotation(
    title = "Historical and Forecasted Dengue Cases in Iloilo",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 5. Plot of predicted data only (mean and confidence interval) [Figure 7]
df_forecast_with_ci <- data.frame(
  Time = time_vals,
  Mean = mean_orig,
  Lower = lower_bound,
  Upper = upper_bound
)

ggplot(df_forecast_with_ci, aes(x = Time, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#1874CD", alpha = 0.2) +
  geom_line(color = "#1874CD", linewidth = 1) +
  scale_x_continuous(
    breaks = c(2023, 2023.25, 2023.5, 2023.75, 2024, 2024.25, 2024.5, 2024.75, 2025),
    labels = c("2023-01", "2023-04", "2023-07", "2023-10",
               "2024-01", "2024-04", "2024-07", "2024-10", "2025-01")
  ) +
  labs(x = "Time", y = "Cases") +
  theme_minimal() +
  theme(legend.position = "none") +
  plot_annotation(
    title = "Forecasted Dengue Cases in Iloilo with Confidence Intervals",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 6. Cumulative predicted dengue cases in 2024
df_confirmation <- data.frame(
  Week = time(complete_ts[833:884]),
  CurrentCases = df_forecast_with_historical[833:884,"Cases"],
  CumulativeCases = cumsum(df_forecast_with_historical[833:884,"Cases"])
)






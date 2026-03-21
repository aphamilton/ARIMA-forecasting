# Introduction ------------
#The goal of this program is to create a model for forecasting custodial populations.

setwd ("C:/Users/Arthur/Documents/R/Custodial Forecasting Project/Data Files/")
getwd()

# Loading libraries and data --------------------

#install.packages("VGAMextra", dependencies = TRUE)

library(ggplot2)
library(urca)
library(TSstudio)
library(aTSA)
library(tseries)
library(forecast)
library(ggfortify)
library(lmtest)
library(bayesforecast)
#library(VGAMextra)
#library(DescTools)
#library(ggthemes)
#library(fable)
#library(reshape2)

#Make all plot titles centered horizontally
theme_update(plot.title = element_text(hjust = 0.5))

load_dataset <- function(file_name)
{
  loaded_data <- read.csv(file_name, fileEncoding = "UTF-8-BOM",
                              header = TRUE)
  return(loaded_data)
}

#Data for 2002-2018 - I am no longer using it
#sc_init_data = load_dataset("StatsCan_initial_data.csv")

#Data for linear regression
sc_all_variables = load_dataset("StatsCan_all_variables_1979_2011.csv")
sc_fewer_variables = load_dataset("StatsCan_fewer_variables_1979_2018.csv")

#Data for ARIMA
actual_in_data = load_dataset("Actual_in_count_SK.csv")

# My own assumption testing functions ----------------------

acf_plot_xticks = seq(0,16,2)
acf_plot_yticks = seq(-1,1,0.2)

ggacf_and_ggpacf <- function(data)
{
  acf_plot <- ggacf(data, title=NULL) +
                    scale_x_continuous(breaks = acf_plot_xticks) +
                    scale_y_continuous(breaks = acf_plot_yticks)
  pacf_plot <- ggpacf(data, title=NULL) +
                      scale_x_continuous(breaks = acf_plot_xticks) +
                      scale_y_continuous(breaks = acf_plot_yticks)
  acf_and_pacf_plots <- list(acf_plot, pacf_plot)
  acf_and_pacf_plots
  return(acf_and_pacf_plots)
}

my_ggqqplot <- function(data)
{
  df <- data.frame(data)
  ggplot(df, aes(sample = data)) +
    stat_qq(colour="blue", size=2) + stat_qq_line(colour="purple", size=1) +
    labs(x="Theoretical Quantiles", y="Sample Quantiles")
}

test_stationarity <- function(time_series)
{
  acf_and_pacf_plots <- ggacf_and_ggpacf(time_series)
  pp_test_results <- pp.test(time_series)
  kpss_test_results <- kpss.test(time_series)
  stationarity_tests <- list(acf_and_pacf_plots,
                             pp_test_results,
                             kpss_test_results)
  stationarity_tests
}

# My own ARIMA functions ----------------------------

aicc <- function(arima_model)
{
  n = arima_model$nobs
  p = length(arima_model$coef)
  aicc = arima_model$aic + 2*p*(p+1)/(n-p-1)
  return(aicc)
}

print_ics <- function(arima_model, time_series)
{
  aic=arima_model$aic
  aicc=aicc(arima_model)
  bic=AIC(arima_model, k=log(length(time_series)))
  return(aicc)
}

my_arima_drift <- function(time_series, p, d, q)
{
  arima_model = Arima(time_series, c(p,d,q), method="ML", include.drift = TRUE)
  aicc = print_ics(arima_model, time_series)
  return(aicc)
}

my_arima_no_drift <- function(time_series, p, d, q)
{
  arima_model = Arima(time_series, c(p,d,q), method="ML", include.drift = FALSE)
  aicc = print_ics(arima_model, time_series)
  return(aicc)
}

find_best_arima_aicc <- function(time_series, d, max_p, max_q)
{
  best_aicc = 1000000
  best_ar = 1000000
  best_ma = 1000000
  best_has_drift = TRUE
  
  #Test models without drift
  for (ar_order in 0:max_p)
  {
    for (ma_order in 0:max_q)
    {
      current_aicc = my_arima_no_drift(time_series, ar_order, d, ma_order)
      if (current_aicc < best_aicc)
      {
        best_aicc = current_aicc
        best_ar = ar_order
        best_ma = ma_order
        best_has_drift = FALSE
      }
    }
  }
  
  #Test models with drift (if d < 2)
  if (d < 2)
  {
    for (ar_order in 0:max_p)
    {
      for (ma_order in 0:max_q)
      {
        current_aicc = my_arima_drift(time_series, ar_order, d, ma_order)
        if (current_aicc < best_aicc)
        {
          best_aicc = current_aicc
          best_ar = ar_order
          best_ma = ma_order
          best_has_drift = TRUE
        }
      }
    }
  }
  
  best_arima <- list("AICc" = best_aicc, "AR" = best_ar,
                     "MA" = best_ma, "Drift" = best_has_drift)
  return(best_arima) 
}

print_best_arima <- function(best_arima_output)
{
  print("AR of best ARIMA model:")
  print(best_arima_output$"AR")
  print("MA of best ARIMA model:")
  print(best_arima_output$"MA")
  print("AICc of best ARIMA model:")
  print(best_arima_output$"AICc")
  print("Best ARIMA model has drift?")
  print(best_arima_output$"Drift")
}

# Simple linear regression (1979-2018) ------------------

sc_simple_lm = lm(Incarceration_rate ~ Year, data = sc_fewer_variables)
summary(sc_simple_lm)
sc_simple_fitted = predict.lm(sc_simple_lm)

ggplot(data=sc_fewer_variables, aes(Year, Incarceration_rate)) +
  geom_point(color="blue", size=2) +
  geom_line(aes(y=sc_simple_fitted), color="purple", size=1) +
  labs(y="Incarceration rate per 100,000")

sc_simple_res = resid(sc_simple_lm)
ggacf_and_ggpacf(sc_simple_res)
dwtest(sc_simple_lm)
my_ggqqplot(sc_simple_res)
shapiro.test(sc_simple_res)

# Multiple linear regression on all variables (1979-2011) -----------
sc_all_lm1 = lm(Incarceration_rate ~ Year + GINI + CPI +
                 Low_income_percent + Police_rate + YM_unemployment_rate,
               data = sc_all_variables)
summary(sc_all_lm1)
sc_all_fitted1 = predict.lm(sc_all_lm1)

sc_all_lm2 = lm(Incarceration_rate ~ Year + CPI +
                  Low_income_percent + Police_rate + YM_unemployment_rate,
                data = sc_all_variables)
summary(sc_all_lm2)
sc_all_fitted2 = predict.lm(sc_all_lm2)

sc_all_lm3 = lm(Incarceration_rate ~ Year + CPI +
                  Low_income_percent + YM_unemployment_rate,
                data = sc_all_variables)
summary(sc_all_lm3)
sc_all_fitted3 = predict.lm(sc_all_lm3)

sc_all_lm4 = lm(Incarceration_rate ~ Year + CPI + YM_unemployment_rate,
                data = sc_all_variables)
summary(sc_all_lm4)
sc_all_fitted4 = predict.lm(sc_all_lm4)

sc_all_lm5 = lm(Incarceration_rate ~ Year + CPI, data = sc_all_variables)
summary(sc_all_lm5)
sc_all_fitted5 = predict.lm(sc_all_lm5)

cor(sc_all_variables)
cor.test(sc_all_variables$Year, sc_all_variables$CPI)

ggplot(data=sc_all_variables, aes(Year, Incarceration_rate)) +
  geom_point(color="blue", size=2) +
  geom_line(aes(y=sc_all_fitted5), color="purple", size=1) +
  labs(y="Incarceration rate per 100,000")

sc_all_res5 = resid(sc_all_lm5)
ggacf_and_ggpacf(sc_all_res5)
dwtest(sc_all_lm5)
my_ggqqplot(sc_all_res5)
shapiro.test(sc_all_res5)

# Multiple linear regression on fewer variables (1979-2018) -----------
sc_fewer_lm1 = lm(Incarceration_rate ~ Year + GINI + CPI +
                  Police_rate + YM_unemployment_rate,
                  data = sc_fewer_variables)
summary(sc_fewer_lm1)
sc_fewer_fitted1 = predict.lm(sc_fewer_lm1)

sc_fewer_lm2 = lm(Incarceration_rate ~ Year + CPI +
                  Police_rate + YM_unemployment_rate,
                  data = sc_fewer_variables)
summary(sc_fewer_lm2)
sc_fewer_fitted2 = predict.lm(sc_fewer_lm2)

sc_fewer_lm3 = lm(Incarceration_rate ~ Year + CPI + Police_rate,
                  data = sc_fewer_variables)
summary(sc_fewer_lm3)
sc_fewer_fitted3 = predict.lm(sc_fewer_lm3)

sc_fewer_lm4 = lm(Incarceration_rate ~ Year + CPI, data = sc_fewer_variables)
summary(sc_fewer_lm4)
sc_fewer_fitted4 = predict.lm(sc_fewer_lm4)

ggplot(data=sc_fewer_variables, aes(Year, Incarceration_rate)) +
  geom_point(color="blue", size=2) +
  geom_line(aes(y=sc_fewer_fitted4), color="purple", size=1) +
  labs(y="Incarceration rate per 100,000")

cor(sc_fewer_variables)
cor.test(sc_fewer_variables$Year, sc_fewer_variables$CPI)

sc_fewer_res4 = resid(sc_fewer_lm4)
ggacf_and_ggpacf(sc_fewer_res4)
dwtest(sc_fewer_lm4)
my_ggqqplot(sc_fewer_res4)
shapiro.test(sc_fewer_res4)

# ARIMA on total population -------------------

total_ts = ts(actual_in_data$Total, start = c(1978), frequency = 1)
test_stationarity(total_ts)
total_ts1 = diff(total_ts, 1, 1)
test_stationarity(total_ts1)

total_best_arima = find_best_arima_aicc(total_ts, 1, 3, 3)
print_best_arima(total_best_arima)
total_ts_arima = Arima(total_ts, c(total_best_arima$"AR", 1,
                        total_best_arima$"MA"), method="ML",
                        include.drift = total_best_arima$"Drift")
print(total_ts_arima)
coeftest(total_ts_arima)

total_ts_arima_res = resid(total_ts_arima)
ggacf_and_ggpacf(total_ts_arima_res)
Box.test(total_ts_arima_res, lag = 20, type = "Ljung-Box",
         fitdf = total_best_arima$"AR" + total_best_arima$"MA")

total_ts_arima_forecast = forecast(total_ts_arima, h=10, level=c(95))
autoplot(total_ts_arima_forecast, ts.size = 1, ts.colour = 'blue',
         predict.size = 0.75, predict.linetype = 'dashed', predict.colour = 'red',
         conf.int = TRUE, conf.int.fill = 'orange') +
         labs(x="Year", y="Actual-in count")

# ARIMA on remand population ------------------------

remand_ts = ts(actual_in_data$Remand, start = c(1978), frequency = 1)
test_stationarity(remand_ts)

remand_ts1 = diff(remand_ts, 1, 1)
test_stationarity(remand_ts1)

remand_best_arima = find_best_arima_aicc(remand_ts, 1, 3, 3)
print_best_arima(remand_best_arima)
remand_ts_arima = Arima(remand_ts, c(remand_best_arima$"AR", 1,
                        remand_best_arima$"MA"), method="ML",
                        include.drift = total_best_arima$"Drift")
print(remand_ts_arima)
coeftest(remand_ts_arima)

remand_ts_arima_res = resid(remand_ts_arima)
ggacf_and_ggpacf(remand_ts_arima_res)
Box.test(remand_ts_arima_res, lag = 20, type = "Ljung-Box",
         fitdf = remand_best_arima$"AR" + remand_best_arima$"MA")

remand_ts_arima_forecast = forecast(remand_ts_arima, h=10, level=c(95))
autoplot(remand_ts_arima_forecast, ts.size = 1, ts.colour = 'blue',
         predict.size = 0.75, predict.linetype = 'dashed', predict.colour = 'red',
         conf.int = TRUE, conf.int.fill = 'orange') +
         labs(x="Year", y="Actual-in count")

# ARIMA on sentenced population ------------------------

sentenced_ts = ts(actual_in_data$Sentenced, start = c(1978), frequency = 1)
test_stationarity(sentenced_ts)

sentenced_ts1 = diff(sentenced_ts, 1, 1)
test_stationarity(sentenced_ts1)

#Commented out because MA term in (0,1,1) was not stat significant
#sentenced_best_arima = find_best_arima_aicc(sentenced_ts, 1, 3, 3)
#print_best_arima(sentenced_best_arima)
#sentenced_ts_arima = Arima(sentenced_ts, c(sentenced_best_arima$"AR", 1,
#                           sentenced_best_arima$"MA"), method="ML",
#                           include.drift = sentenced_best_arima$"Drift")

sentenced_ts_arima = Arima(sentenced_ts, c(0,1,0), method="ML",
                           include.drift = FALSE)
print(sentenced_ts_arima)

#Commented out because there are no AR/MA terms and there is no drift
#coeftest(sentenced_ts_arima)

sentenced_ts_arima_res = resid(sentenced_ts_arima)
ggacf_and_ggpacf(sentenced_ts_arima_res)
Box.test(sentenced_ts_arima_res, lag = 20, type = "Ljung-Box",
         fitdf = 0)

sentenced_ts_arima_forecast = forecast(sentenced_ts_arima, h=10, level=c(95))
autoplot(sentenced_ts_arima_forecast, ts.size = 1, ts.colour = 'blue',
         predict.size = 0.75, predict.linetype = 'dashed', predict.colour = 'red',
         conf.int = TRUE, conf.int.fill = 'orange') +
         labs(x="Year", y="Actual-in count")

# EVERYTHING BELOW HERE DID NOT END UP BEING USED -----------------

# Old version of stationarity tests (not in use) ---------------

# ts_test_station <- function(time_series)
# {
#   ts_acf = acf(time_series, plot=FALSE)
#   ggplot(data=ts_acf, mapping=aes(x=lag, y=acf)) +
#          geom_bar(stat = "identity", position = "identity") +
#          labs(title=time_series, x="Lag", y="ACF")
#   ts_pacf = pacf(time_series, plot=FALSE)
#   ggplot(data=ts_pacf, mapping=aes(x=lag, y=pacf)) +
#          geom_bar(stat = "identity", position = "identity") +
#          labs(title=time_series, x="Lag", y="PACF")
#   pp.test(time_series)
#   kpss.test(time_series) 
# }

# Old version of finding best ARIMA model (didn't work) -----------

# best_aic = 1000000
# best_aic_ar = 1000000
# best_aic_ma = 1000000
# for (ar_order in 0:3)
# {
#   for (ma_order in 0:3)
#   {
#     arima_result <- arima(total_ts, c(ar_order, 1, ma_order))
#     print("Printing ARIMA result")
#     print(arima_result)
#     if (arima_result.aIc < best_aic)
#     {
#       best_aic = arima_result.aIc
#       best_aic_ar = ar_order
#       best_aic_ma = ma_order
#     }
#   }
# }
#
# print("The best AIC is the following, at the following orders of AR and MA:")
# print(best_aic)
# print(best_aic_ar)
# print(best_aic_ma)

# Regression on initial data (2002-2018) ------------

#Simple linear regression on initial data
# sc_init_lm = lm(Total_custody_per_100000 ~ Year, data = sc_init_data)
# summary(sc_init_lm)
# sc_init_fitted = predict.lm(sc_init_lm)
# print(sc_init_fitted)
# 
# plot(sc_init_data$Total_custody_per_100000, sc_init_fitted)
# ggplot(data=sc_init_data, aes(Year, Total_custody_per_100000)) +
#       geom_point(color="blue", size=2) +
#       geom_line(aes(y=sc_init_fitted), color="purple", size=1)

#Multiple linear regression on initial data
# sc_init_lm = lm(Total_custody_per_100000 ~ Year + Real_GDP +
#                   Unemployment_rate + CPI, data = sc_init_data)
# summary(sc_init_lm)
# sc_init_fitted = predict.lm(sc_init_lm)
# print(sc_init_fitted)
# 
# #print(length(stats_can_fitted))
# #print(length(sc_init_data$Total_custody_per_100000))
# 
# plot(sc_init_data$Total_custody_per_100000, sc_init_fitted)
# ggplot(data=sc_init_data, aes(Year, Total_custody_per_100000)) +
#   geom_point(color="blue", size=2) +
#   geom_line(aes(y=sc_init_fitted), color="purple", size=1)

# Attempts (unsuccessful) to generalize regression and plotting -------

# slr_and_plot <- function(dataset, dv, iv)
# {
#   print("Dataset: ")
#   print(dataset)
#   #print("dv: ")
#   #print(dv)
#   print("dataset$dv: ")
#   print(dataset$dv)
#   the_lm = lm(dv ~ iv, data = dataset)
#   print("Got to here")
#   summary(the_lm)
#   fitted_values = predict.lm(the_lm)
#   print("Fitted values: ")
#   print(fitted_values)
#   
#   
#   plot(dataset$dv, fitted_values)
#   ggplot(data = dataset$dv, aes(iv, dv)) +
#                     geom_point(color = "blue", size = 2) +
#                     geom_line(aes(y = fitted_values), color = "purple", size = 1)
#   
#   return(the_lm)
# }
# 
# mlr_all_ivs_and_plot <- function(dataset, dv, iv_to_plot)
# {
#   the_lm = lm(dv ~ . - dv, data = dataset)
#   summary(the_lm)
#   fitted_values = predict.lm(the_lm)
#   print(fitted_values)
#   
#   #plot(dataset$dv, fitted_values)
#   ggplot(data = dataset$dv, aes(iv_to_plot, dv)) +
#         geom_point(color = "blue", size = 2) +
#         geom_line(aes(y = fitted_values), color = "purple", size = 1)
# }
# 
# sc_init_lm = slr_and_plot(sc_init_data,
#                           Total_custody_per_100000,
#                           Year)
#summary(sc_init_lm)

#mlr_and_plot(sc_init_data,
#             sc_init_data$Total_custody_per_100000,
#             sc_init_data$Year)

#mlr_all_ivs_and_plot(sc_all_variables,
#                     sc_all_variables$Incarceration_rate,
#                     sc_all_variables$Year)

# Plot of the total actual-in count time series with ggplot -------------------

# only_every_fourth <- function(values)
# {
#   values[seq(from = 2, to = length(values), by = 2)] <- ""
#   values[seq(from = 3, to = length(values), by = 4)] <- ""
#   #values[seq(from = 5, to = length(values), by = 8)] <- ""
#   return(values)
# }
#
# year_labels = only_every_fourth(actual_in_data$Year)
#
# more_year_labels = seq(1, 41, 4)
# 
# ggplot(data=actual_in_data, aes(x=Year, y=Total)) +
#        geom_point(color="blue", size=2) +
#        scale_x_discrete(breaks = actual_in_data$Year[more_year_labels],
#                         labels = actual_in_data$Year[more_year_labels]) +
#        scale_y_continuous(name = "Total Custodial Population")# +
#        #coord_cartesian(ylim = c(500,2500))

# Unused alternate ACF/PACF plots ---------------------

#total_ts_acf = acf(total_ts, plot=FALSE)
#total_ts_acf_df <- with(total_ts_acf, data.frame(lag, acf))
#ggplot(data=total_ts_acf_df, mapping=aes(x=lag, y=acf)) +
#  geom_bar(stat = "identity", position = "identity", fill = "royalblue2") +
#  labs(title="ACF Plot of the Total Time Series at I(0)", x="Lag", y="ACF")
#total_ts_pacf = pacf(total_ts, plot=FALSE)
#autoplot(total_ts_pacf)
#total_ts_pacf_df <- with(total_ts_pacf, data.frame(lag, pacf))
#ggplot(data=total_ts_pacf, mapping=aes(x=lag, y=pacf)) +
#  geom_bar(stat = "identity", position = "identity") +
#  labs(title="PACF plot of total_ts", x="Lag", y="PACF")
#
#total_ts1_acf = acf(total_ts1, plot=FALSE)
#total_ts1_acf_df <- with(total_ts1_acf, data.frame(lag, acf))
#ggplot(data=total_ts1_acf_df, mapping=aes(x=lag, y=acf)) +
#  geom_bar(stat = "identity", position = "identity", fill = "royalblue2") +
#  labs(title="ACF Plot of the Total Time Series at I(1)", x="Lag", y="ACF")
#total_ts1_pacf = pacf(total_ts1, plot=FALSE)
#autoplot(total_ts1_pacf)
#total_ts1_pacf_df <- with(total_ts1_pacf, data.frame(lag, pacf))
#ggplot(data=total_ts1_pacf, mapping=aes(x=lag, y=pacf)) +
#  geom_bar(stat = "identity", position = "identity") +
#  labs(title="PACF plot of total_ts1", x="Lag", y="PACF")

# Attempting alternate ARIMA models -----------------------

# #Attempting ARIMA(1,1,0) model
# sentenced_ts_arima110 = Arima(sentenced_ts, c(1,1,0), method="ML",
#                               include.drift = FALSE)
# sentenced_ts_arima110_res = resid(sentenced_ts_arima110)
# acf(sentenced_ts_arima110_res)
# pacf(sentenced_ts_arima110_res)
# Box.test(sentenced_ts_arima110_res, lag=20, type="Ljung-Box", fitdf=1)
# sentenced_ts_arima110_forecast = forecast(sentenced_ts_arima110,
#                                           h=10, level=c(95))
# plot(sentenced_ts_arima110_forecast)
# 
# #Attempting ARIMA(0,1,1) model
# sentenced_ts_arima011 = Arima(sentenced_ts, c(0,1,1), method="ML",
#                               include.drift = FALSE)
# sentenced_ts_arima011_res = resid(sentenced_ts_arima011)
# acf(sentenced_ts_arima011_res)
# pacf(sentenced_ts_arima011_res)
# Box.test(sentenced_ts_arima011_res, lag=20, type="Ljung-Box", fitdf=1)
# sentenced_ts_arima011_forecast = forecast(sentenced_ts_arima011,
#                                           h=10, level=c(95))
# plot(sentenced_ts_arima011_forecast)

# #Attempting ARIMA(1,1,0) model
# remand_ts_arima110 = Arima(remand_ts, c(1,1,0), method="ML",
#                            include.drift = TRUE)
# remand_ts_arima110_res = resid(remand_ts_arima110)
# acf(remand_ts_arima110_res)
# pacf(remand_ts_arima110_res)
# Box.test(remand_ts_arima110_res, lag=20, type="Ljung-Box", fitdf=1)
# remand_ts_arima110_forecast = forecast(remand_ts_arima110,
#                                        h=10, level=c(95))
# plot(remand_ts_arima110_forecast)
# 
# #Attempting ARIMA(0,1,1) model
# remand_ts_arima011 = Arima(remand_ts, c(0,1,1), method="ML",
#                            include.drift = TRUE)
# remand_ts_arima011_res = resid(remand_ts_arima011)
# acf(remand_ts_arima011_res)
# pacf(remand_ts_arima011_res)
# Box.test(remand_ts_arima011_res, lag=20, type="Ljung-Box", fitdf=1)
# remand_ts_arima011_forecast = forecast(remand_ts_arima011,
#                                        h=10, level=c(95))
# plot(remand_ts_arima011_forecast)

# Plotting all 3 time series together w/ forecasts (unsuccessfully) ---------------

# autoplot(remand_ts_arima_forecast, ts.size = 1, ts.colour = 'blue',
#          predict.size = 0.75, predict.linetype = 'dashed',
#          conf.int = TRUE, conf.int.fill = 'blue') +
#   autolayer(sentenced_ts_arima_forecast, ts.size = 1, ts.colour = 'red',
#             predict.size = 0.75, predict.linetype = 'dashed',
#             conf.int = TRUE, conf.int.fill = 'red', alpha = 0.5) +  
#   autolayer(total_ts_arima_forecast, ts.size = 1, ts.colour = 'green',
#             predict.size = 0.75, predict.linetype = 'dashed',
#             conf.int = TRUE, conf.int.fill = 'green', alpha = 0.5) +
#   labs(y="Actual-in count") +
#   labs(x="Year") +
#   ggtitle("Ten-Year Forecast for SK Custodial Population")

# the_year_labels = seq(1, 41, 4)
# ggplot(data=actual_in_data) +
#         geom_line(aes(x=Year, y=Total), color="blue", size=1) +
#         scale_x_discrete(breaks = actual_in_data$Year[the_year_labels],
#                           labels = actual_in_data$Year[the_year_labels]) +
#         scale_y_continuous(name = "Total Custodial Population") +
#         coord_cartesian(ylim = c(0,2500))

#joint_forecasts <- rbind(actual_ts_arima_forecast,
#                         sentenced_ts_arima_forecast,
#                         remand_ts_arima_forecast)

#plot(joint_forecasts)
#ggplot(joint_forecasts, aes(x, y, group=model)) + geom_point()

#plot(actual_ts_arima_forecast)
#plot(sentenced_ts_arima_forecast)
#plot(remand_ts_arima_forecast)
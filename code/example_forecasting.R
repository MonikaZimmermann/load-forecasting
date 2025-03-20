# Load Packages ---------------------------------------------------------------------------------
library(data.table)
library(here) # Finds the repository root
library(lubridate) # Date-time manipulation
library(ggplot2) # Visualization
library(colorBlindness)
library(scales)
library(mgcv) # GAM
library(mgcViz) # GAM visualization
library(smooth) # for ETS models
library(forecast) # for ETS models
library(zoo)
library(dplyr)
library(tidyr)


# Load Data -------------------------------------------------------------------------------------
# define zones
zones <- c("FR", "DE") # for the example we chose Germany and France
# load csv file data for these zones from data folder
data_path <- here("data")
csv_files <- setNames(file.path(data_path, paste0(zones, "_data.csv")), zones)
zone_data_list <- lapply(csv_files, fread)


# Set Hyperparameter -------------------------------------------------------------------------------------
S <- 24 # for hourly data
H <- 168 * 52 # forecasting horizon
D <- H * 4 # in-sample length

# for temperature model
talpha <- c(1, 14) * S # smoothing parameters for temperatures
order.max <- 8 * 168 # maximum order for AR

# for states-ETS model
ltimestamps <- list(c(
    seq(8, 19), # peak hours Monday
    seq(8, 19) + 24, # Tuesday
    seq(8, 19) + 24 * 2, # Wednesday
    seq(8, 19) + 24 * 3, # Thursday
    seq(8, 19) + 24 * 4 # Friday
)) # peak hours in the week for temporal aggregation of first stage GAM residuals
freq <- 168 # frequency reduction, e.g. for hourly data with daily slicing freq=24, with weekly slicing freq=168
hADAM <- 2 # in-sample TMAE-loss horizon
lagsADAM <- 52 # 1 for ANN model, period of seasonal component, otherwise, e.g. lagsADAM=52, does not allow multiple frequencies

# for GAM model
gamma <- log(D) / 2 # for smoother fits  in a BIC-like range

rho <- 0
# number of knots
khldfixed <- 8
khldweekday <- 8
khldp <- 8
ktemp <- 6
kHoD <- 24
kHoW <- 168
kHoY <- 12
kHoYHoD <- c(12, 8)
kHoYHoW <- c(12, 7)
kstates <- 8
# ensure integer knot placement
knots <- list(HoD = c(0, 24), HoY = c(0, 8760), HoW = c(0, 168))
knotsTemp <- list(HoY_utc = c(0, 8760), HoD_utc = c(0, 24))



# Forecasting -------------------------------------------------------------------------------------
# loop through zones and save the forecasted load, GAM and AR model
models <- list()
models_AR <- list()
load_forecast <- list()
for (zn in zones) {
    index <- 1:D
    data <- zone_data_list[[zn]][index, ]
    data_ext <- zone_data_list[[zn]][1:(D + H), ]
    # get impact
    IMPACT <- get.impact(data$time_local_tz[c(index, 1:H + tail(index, 1))], data$load_local_tz[c(index, 1:H + tail(index, 1))], idy = 1:D)
    # temp  model
    temp_lt <- forecast_temp(S = 24, talpha = talpha, temp = data$temp_utc, HoD_utc = data_ext$HoD_utc, HoY_utc = data_ext$HoY_utc, time_utc = data_ext$time_utc, H = H, HoD = data_ext$HoD_local_tz, knots = knotsTemp)$temp_lt
    # holiday columns
    HLDfixed_cols <- grep("HLDfixed_", colnames(data_ext), value = TRUE)
    HLDweekday_cols <- grep("HLDweekday_", colnames(data_ext), value = TRUE)
    HLDP_cols <- grep("HLDP_", colnames(data_ext), value = TRUE)
    HLDfixed <- data_ext[, ..HLDfixed_cols]
    HLDweekday <- data_ext[, ..HLDweekday_cols]
    HLDP <- data_ext[, ..HLDP_cols]
    # load model
    load <- forecast_AR_GAM(
        Y = data$load_local_tz, HoD = data_ext$HoD_local_tz, HoW = data_ext$HoW_local_tz, HoY = data_ext$HoY_local_tz, HLDfixed = HLDfixed, HLDweekday = HLDweekday, HLDP = HLDP, H = H, S = S, impact = IMPACT,
        order.max = order.max, temp_in = TRUE, temp_lt = temp_lt, ltimestamps = ltimestamps, ADAM = TRUE, freq = freq, lagsADAM = lagsADAM, hADAM = hADAM, rho = rho, gamma = gamma, khldfixed = khldfixed, khldweekday = khldweekday, khldp = khldp, ktemp = ktemp, kHoD = kHoD, kHoW = kHoW, kHoY = kHoY, kHoYHoW = kHoYHoW, kHoYHoD = kHoYHoD, kstates = kstates, knots = knots
    )
    load_forecast[[zn]] <- load$forecasts
    models[[zn]] <- load$mod
    models_AR[[zn]] <- load$modres
}


# Plotting -------------------------------------------------------------------------------------
# plot the forecasted and true load for the one year horizon

# define the path to the 'plots'
plot_path <- here("plots")

for (zn in zones) {
    plot_data <- data.frame(
        Time = zone_data_list[[zn]][(D + 1):(D + H), ]$time_local_tz,
        True = zone_data_list[[zn]][(D + 1):(D + H), ]$load_local_tz,
        Forecasted = load_forecast[[zn]]
    )

    p <- ggplot(plot_data, aes(x = Time)) +
        geom_line(aes(y = True, color = "True"), linewidth = 0.3, alpha = 0.3) +
        geom_line(aes(y = Forecasted, color = "Forecasted"), linewidth = 0.3, alpha = 0.3, linetype = "dashed") +
        scale_x_datetime(
            breaks = seq(min(plot_data$Time), max(plot_data$Time), by = "4 weeks"),
            labels = date_format("%b %Y"),
            expand = c(0, 0)
        ) +
        scale_y_continuous(name = "Load [GW]") +
        scale_color_manual(name = "Load", values = c("True" = "black", "Forecasted" = "blue")) +
        labs(title = "Forecasted and True Load - One-Year Horizon", x = "Time") +
        theme_minimal() +
        theme_big +
        theme(
            legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(size = 0.2, linetype = "dashed", colour = "lightgrey"),
            panel.grid.minor = element_blank()
        )

    ggsave(filename = file.path(plot_path, paste0("Forecast_vs_True_", zn, ".pdf")), plot = p, width = 15, height = 5)
}




# Plotting -------------------------------------------------------------------------------------
# plot the fitted load disaggregated by the components (without ar)

# define the path to the 'plots'
plot_path <- here("plots")

# define winter holiday periods
# 2015-2016
start_date_2015 <- as.POSIXct("2015-12-15", tz = "UTC")
end_date_2015 <- as.POSIXct("2016-01-15", tz = "UTC")
# 2016-2017
start_date_2016 <- as.POSIXct("2016-12-15", tz = "UTC")
end_date_2016 <- as.POSIXct("2017-01-15", tz = "UTC")
# 2017-2018
start_date_2017 <- as.POSIXct("2017-12-15", tz = "UTC")
end_date_2017 <- as.POSIXct("2018-01-15", tz = "UTC")

# define specific holiday dates
holidays <- c("2015-12-24", "2015-12-25", "2015-12-26", "2015-12-31", "2016-01-01", "2016-01-06")

# create holiday shading rectangles
holiday_rects <- data.frame(
    xmin = as.POSIXct(holidays, tz = "UTC"),
    xmax = as.POSIXct(holidays, tz = "UTC") + days(1),
    ymin = -Inf,
    ymax = Inf
)

# define plot colors (color-blind friendly)
col_hol <- Green2Magenta16Steps[c(9, 10, 12, 15)]
col_seas <- c("#abf6f7", "#bef1f2")
col_temp <- SteppedSequential5Steps[c(1, 4)]
col_level <- SteppedSequential5Steps[c(11, 14)]
colors <- c(col_hol[4], col_hol[3], col_hol[2], col_seas[1], col_temp[1], col_level[1])

# define big_theme for plots on page width

theme_big <- theme(
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)
)

for (zn in zones) {
    # extract smooth term labels
    term_labels <- unlist(lapply(models[[zn]]$smooth, function(x) x$label))

    # get smooth effects
    plot_data <- as.data.frame(predict.gam(models[[zn]], type = "terms", se.fit = TRUE)$fit)

    # aggregate smooth effects
    # socio-economic effects
    plot_data <- plot_data %>%
        mutate(
            Level = rowSums(select(., starts_with("s(HoD):")))
        ) %>%
        select(., -starts_with("s(HoD):"))

    # holiday effects
    plot_data <- plot_data %>%
        mutate(
            HLDF = rowSums(select(., ends_with("impact")))
        ) %>%
        select(., -ends_with("impact"))

    if (zn == "DE") {
        plot_data <- plot_data %>%
            mutate(
                HLDW = rowSums(select(., ends_with("ay)"))) + as.numeric(`s(HLDweekday_local_tz.CorpusChristi)`) + as.numeric(`s(HLDweekday_local_tz.Pentecost)`)
            ) %>%
            select(., -ends_with("ay)"), -`s(HLDweekday_local_tz.CorpusChristi)`, -`s(HLDweekday_local_tz.Pentecost)`)
    } else {
        plot_data <- plot_data %>%
            mutate(
                HLDW = rowSums(select(., ends_with("ay)")))
            ) %>%
            select(., -ends_with("ay)"))
    }
    # rename the "s(HLDP_Local_tz)" column to HLDP
    plot_data <- plot_data %>%
        rename(HLDP = `s(HLDP_Local_tz)`)

    # temp effects
    plot_data <- plot_data %>%
        mutate(
            Temp = rowSums(select(., starts_with("s(temp")))
        ) %>%
        select(., -starts_with("s(temp"))

    # seasonal effects
    plot_data <- plot_data %>%
        mutate(
            Seas = rowSums(select(., starts_with("s(Ho")))
        ) %>%
        select(., -starts_with("s(Ho"))

    plot_data <- plot_data %>%
        mutate(
            SeasInt = rowSums(select(., starts_with("ti")))
        ) %>%
        select(., -starts_with("ti"))

    plot_data <- plot_data %>%
        mutate(
            Periodics = rowSums(select(., starts_with("Seas")))
        ) %>%
        select(., -starts_with("Seas"))

    # align time
    Time <- zone_data_list[[zn]][index, ]$time_local_tz

    plot_data <- plot_data %>%
        mutate(Time = Time, Intercept = models[[zn]]$coefficients[1])

    # compute actual and fitted load
    plot_data <- plot_data %>%
        mutate(
            Load_act = models[[zn]]$fitted.values + models[[zn]]$residuals - Intercept,
            Load_fit = models[[zn]]$fitted.values - Intercept
        )

    # subset data for winter holiday periods
    subset_data <- rbind(
        plot_data %>%
            filter(Time >= start_date_2015 & Time <= end_date_2015),
        plot_data %>%
            filter(Time >= start_date_2016 & Time <= end_date_2016),
        plot_data %>%
            filter(Time >= start_date_2017 & Time <= end_date_2017)
    )
    subset_data <- subset_data %>% mutate(Time = rep(seq(from = start_date_2015, to = end_date_2015, by = "1 hour"), 3))

    subset_data <- subset_data %>% mutate(Year = c(rep("2015-2016", nrow(subset_data) / 3), rep("2016-2017", nrow(subset_data) / 3), rep("2017-2018", nrow(subset_data) / 3)))

    # reshape the data into long format
    long_data <- subset_data %>%
        pivot_longer(cols = -c(Time, Intercept, Year, Load_fit, Load_act), names_to = "Effects", values_to = "Value1")

    long_data <- long_data %>%
        pivot_longer(cols = -c(Time, Intercept, Year, Effects, Value1), names_to = "Load", values_to = "Value2")

    # reorder more sensibly
    long_data$Effects <- factor(long_data$Effects, levels = c("HLDW", "HLDF", "HLDP", "Periodics", "Temp", "Level"))

    long_data$Load <- factor(long_data$Load, levels = c("Load_act", "Load_fit"))

    # Create stacked plot
    p_stack <- ggplot(long_data, aes(x = Time, y = Value1, fill = Effects)) +
        geom_rect(data = holiday_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.7, inherit.aes = FALSE) +
        geom_area() +
        geom_line(aes(x = Time, y = Value2, color = Load), linewidth = 0.2, alpha = 0.9, linetype = "solid") +
        scale_fill_manual(values = colors, labels = c("HldW", "HldF", "HldP", "Seasonal (incl. interactions)", "Temp", "States")) +
        scale_color_manual(values = c("#7d7a7afc", "black"), labels = c("Actual", "Fit")) +
        facet_wrap(~Year, nrow = 3, scales = "free_y", strip.position = "right") +
        scale_x_datetime(labels = scales::date_format("%b %d"), breaks = seq(min(long_data$Time), max(long_data$Time), by = "2 days"), expand = c(0, 0)) +
        scale_y_continuous(name = "Load [GW]", breaks = seq(-30, 20, by = 5), expand = c(0, 0)) +
        theme_minimal() +
        theme_big +
        theme(
            legend.position = "top",
            axis.title.x = element_blank(),
            legend.key.width = unit(0.5, "lines"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.25, linetype = "dashed", colour = "lightgrey"),
            panel.grid.minor = element_line(size = 0.25, linetype = "dashed", colour = "lightgrey")
        ) +
        guides(fill = guide_legend(nrow = 1), color = guide_legend(override.aes = list(linewidth = 0.45)))

    # save plot
    filename <- paste0(plot_path, "/Stack_Winter_plot_", zn, ".pdf")
    ggsave(filename, plot = p_stack, width = 6, height = 5)
}

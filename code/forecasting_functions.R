# DST.trafo: function for clock adjustment on regression vectors ------------------------------

DST.trafo <- function(X, Xtime, Xtz = "CET", freq = NULL) {
    ## Xtime in UTC - ONLY for EUROPEAN DATA with change at UTC+1
    Xinit <- as.matrix(X)
    atime.init <- as.numeric(Xtime)
    if (is.null(freq)) freq <- as.numeric(names(which.max(table(diff(atime.init)))))
    S <- 24 * 60 * 60 / freq
    atime <- seq(atime.init[1], atime.init[length(atime.init)], freq)
    idmatch <- match(atime.init, atime)
    X <- array(, dim = c(length(atime), dim(Xinit)[2]))
    X[idmatch, ] <- Xinit


    xx <- as.POSIXct(atime, origin = "1970-01-01", tz = Xtz)
    days <- unique(as.Date(xx, tz = Xtz))
    if (Xtz %in% c("WET", "CET", "EET")) { # EUROPE
        DST.SPRING <- grepl("030", format(days, "%m%w")) & as.numeric(format(days, "%d")) > 24
        DST.FALL <- grepl("100", format(days, "%m%w")) & as.numeric(format(days, "%d")) > 24
    }
    if (Xtz %in% c("US/Hawaii", "US/Alaska", "US/Pacific", "UC/Mountain", "US/Central", "US/Eastern")) { # N-America
        DST.SPRING <- grepl("030", format(days, "%m%w")) & as.numeric(format(days, "%d")) > 7 & as.numeric(format(days, "%d")) < 15
        DST.FALL <- grepl("110", format(days, "%m%w")) & as.numeric(format(days, "%d")) < 8
    }
    DST <- !(DST.SPRING | DST.FALL)

    xxf.start <- format(xx[1:S], tz = Xtz, usetz = TRUE)
    xxf.end <- format(xx[seq.int(length(xx) - S + 1, length(xx))], tz = Xtz, usetz = TRUE)


    DLf <- format(days)
    Dlen <- length(DLf)

    Shift <- 2 # for CET
    if (Xtz == "WET") Shift <- 1
    if (Xtz == "EET") Shift <- 3

    Xout <- array(, dim = c(Dlen, S, dim(X)[2]))
    k <- 0
    ## first entry:
    i.d <- 1
    idx <- grep(DLf[i.d], xxf.start)
    if (DST[i.d]) {
        Xout[i.d, S + 1 - rev(1:length(idx)), ] <- X[k + 1:length(idx), ]
    } else {
        # 	print(i.d)
        tmp <- (S - rev(idx) + 1)
        if (DST.SPRING[i.d]) {
            ## MARCH
            for (i.S in seq_along(idx)) {
                if (tmp[i.S] <= Shift * S / 24) Xout[i.d, S - S / 24 - length(idx) + i.S, ] <- X[k + i.S, ]
                if (tmp[i.S] == Shift * S / 24) Xout[i.d, S - S / 24 - length(idx) + i.S + 1:(S / 24), ] <- t(X[k + i.S, ] + t(tcrossprod((1:(S / 24)) / (length(1:(S / 24)) + 1), (X[k + i.S + 1, ] - X[k + i.S, ]))))
                if (tmp[i.S] > Shift * S / 24) Xout[i.d, S - S / 24 - length(idx) + i.S + S / 24, ] <- X[k + i.S, ]
            } # i.S
        } else {
            ## October
            for (i.S in seq_along(idx)) {
                if (tmp[i.S] <= Shift * S / 24) Xout[i.d, S + S / 24 - length(idx) + i.S, ] <- X[k + i.S, ]
                if (tmp[i.S] %in% (Shift * S / 24 + 1:(S / 24))) Xout[i.d, S + S / 24 - length(idx) + i.S, ] <- 0.5 * (X[k + i.S, ] + X[k + i.S + S / 24, ])
                if (tmp[i.S] > (Shift + 2) * S / 24) Xout[i.d, S + S / 24 - length(idx) + i.S - S / 24, ] <- X[k + i.S, ]
            } # i.S
        }
    } #
    k <- k + length(idx)
    for (i.d in seq_along(DLf)[c(-1, -length(DLf))]) { ## first and last extra
        idx <- 1:S
        if (DST[i.d]) {
            Xout[i.d, 1:length(idx), ] <- X[k + 1:length(idx), ]
        } else {
            if (DST.SPRING[i.d]) {
                idx <- 1:(S - S / 24)
                ## MARCH
                for (i.S in seq_along(idx)) {
                    if (i.S <= Shift * S / 24) Xout[i.d, i.S, ] <- X[k + i.S, ]
                    if (i.S == Shift * S / 24) Xout[i.d, i.S + 1:(S / 24), ] <- t(X[k + i.S, ] + t(tcrossprod((1:(S / 24)) / (length(1:(S / 24)) + 1), (X[k + i.S + 1, ] - X[k + i.S, ]))))
                    if (i.S > Shift * S / 24) Xout[i.d, i.S + S / 24, ] <- X[k + i.S, ]
                } # i.S
            } else {
                idx <- 1:(S + S / 24)
                ## October
                for (i.S in seq_along(idx)) {
                    if (i.S <= Shift * S / 24) Xout[i.d, i.S, ] <- X[k + i.S, ]
                    if (i.S %in% (Shift * S / 24 + 1:(S / 24))) Xout[i.d, i.S, ] <- 0.5 * (X[k + i.S, ] + X[k + i.S + S / 24, ])
                    if (i.S > (Shift + 2) * S / 24) Xout[i.d, i.S - S / 24, ] <- X[k + i.S, ]
                } # i.S
            }
        } #
        k <- k + length(idx)
    }
    ## last
    i.d <- length(DLf)
    idx <- grep(DLf[i.d], xxf.end)
    if (DST[i.d]) {
        Xout[i.d, 1:length(idx), ] <- X[k + 1:length(idx), ]
    } else {
        if (DST.SPRING[i.d]) {
            ## MARCH
            for (i.S in seq_along(idx)) {
                if (i.S <= Shift * S / 24) Xout[i.d, i.S, ] <- X[k + i.S, ]
                if (i.S == Shift * S / 24) Xout[i.d, i.S + 1:(S / 24), ] <- t(X[k + i.S, ] + t(tcrossprod((1:(S / 24)) / (length(1:(S / 24)) + 1), (X[k + i.S + 1, ] - X[k + i.S, ]))))
                if (i.S > Shift * S / 24) Xout[i.d, i.S + S / 24, ] <- X[k + i.S, ]
            } # i.S
        } else {
            ## October
            for (i.S in seq_along(idx)) {
                if (i.S <= Shift * S / 24) Xout[i.d, i.S, ] <- X[k + i.S, ]
                if (i.S %in% (Shift * S / 24 + 1:(S / 24))) Xout[i.d, i.S, ] <- 0.5 * (X[k + i.S, ] + X[k + i.S + S / 24, ])
                if (i.S > (Shift + 2) * S / 24) Xout[i.d, i.S - S / 24, ] <- X[k + i.S, ]
            } # i.S
        }
    } #
    k <- k + length(idx)


    Xout
}



# get.tempets: get smoothed temperature as fit of ANN ets model for specified values of alpha ------------------------------

get.tempets <- function(xt, ialpha) {
    xtt <- na.approx(xt, rule = 2, na.rm = FALSE)
    get.ets <- function(z) ets(xtt, "ANN", alpha = 1 / z)$fitted
    TETS <- as.matrix(sapply(ialpha, get.ets))
    dimnames(TETS) <- list(NULL, paste0("temp_", ialpha))
    TETS
}


# forecast_temp: gam + AR model to forecast ets-smoothed temperature -------------------------------------------------

forecast_temp <- function(S = 24, talpha, temp, HoD_utc, HoY_utc, time_utc, H, HoD, knots) {
    # note: temp input in UTC (meteorological data should not be clock change adjusted but modelled in UTC; humans adjust to local time the planet not)
    # fitted and forecasted temperature are adjusted after modelling to local time

    TETS <- get.tempets(temp, talpha)
    TEMPETS <- rbind(TETS, matrix(nrow = H, ncol = dim(TETS)[2]))

    df.temp <- data.frame(TEMPETS, HoD_utc, HoY_utc)

    tempterm <- paste0(
        "+ s(HoD_utc, bs='cp', k=S) + s(HoY_utc, bs='cp', k=6) + ",
        "ti(HoY_utc, HoD_utc, bs='cp', k=c(6,8))"
    )

    for (i.t in 1:dim(TEMPETS)[2]) {
        formdet <- as.formula(paste0(dimnames(TEMPETS)[[2]][i.t], " ~ ", tempterm))

        system.time(modtemp <- bam(
            formdet,
            select = TRUE, data = df.temp[1:length(temp), ],
            discrete = TRUE, gamma = 1, nthreads = 4, knots = knots
        ))

        pp <- predict(modtemp, newdata = tail(df.temp, H), type = "response")

        temp_ar <- ar(residuals(modtemp), order.max = 4 * 168)

        df.temp[(1:H + length(temp)), dimnames(TEMPETS)[[2]][i.t]] <-
            pp + predict(temp_ar, n.ahead = H, se.fit = FALSE)
    }

    temp_lt <- matrix(NA, nrow = length(HoD), ncol = (dim(TEMPETS)[2]))
    colnames(temp_lt) <- colnames(TEMPETS)
    for (i.t in 1:dim(TEMPETS)[2]) {
        temp_s <- as.numeric(t(DST.trafo(
            X = df.temp[, dimnames(TEMPETS)[[2]][i.t]],
            Xtime = time_utc,
            Xtz = "CET"
        )[, , 1]))

        temp_ts <- tail(temp_s, length(temp_s) - (HoD[1]))
        if (tail(HoD, 1) != 23) {
            temp_lt[, i.t] <- head(temp_ts, -(23 - tail(HoD, 1)))
        } else {
            temp_lt[, i.t] <- temp_ts
        }
    }
    for (i.t in 1:dim(TEMPETS)[2]) {
        if (any(is.na(temp_lt[, i.t]) == TRUE)) {
            print("note dst shift")
            print(which(is.na(temp_lt[, i.t])))
        }
    }

    return(list("temp_lt" = temp_lt, "temp_ar" = temp_ar))
}


# get.impact: calculate an hourly impact factor for fixed date holidays -------------------------------------------

get.impact <- function(xtime, y, idy) {
    # xtime: time vector as.POSIXct (output impact will have the same length)
    # y: dependent variable with weekly pattern (not necessarily aligned with xtime)
    # idy: index set in xtime aligning y with xtime

    HoD <- lubridate::hour(xtime[1:max(idy)]) + 1
    DoW <- lubridate::wday(xtime[1:max(idy)], week_start = 1) # DoW starts at 1 for Monday, goes to 7 for Sunday
    HoW <- HoD + 24 * (DoW - 1) # HoW starts with 1 at Monday 0:00

    # index of first Mo 0:00 in y
    fi <- which(HoW[idy][1:(7 * S)] == 1)

    # index of first Mo 0:00 in xtime
    fi_xtime <- which(HoW[1:(7 * S)] == 1)


    # calculate median per HoW
    if (length(y) %% (7 * S) == 0) {
        WK <- apply(matrix(y, 7 * S), 1, median, na.rm = TRUE)
    } else {
        WK <- apply(matrix(head(y, -(length(y) %% (7 * S))), 7 * S), 1, median, na.rm = TRUE)
    }

    if (fi != 1) WK <- WK[c(fi:(S * 7), 1:(fi - 1))]
    WKmax <- apply(matrix(WK[S + 1:(3 * S)], S), 1, mean) # mean of Tue-Thu per HoD
    WKmin <- tail(WK, S) ## Sun
    himpact <- pmin(pmax((WK - WKmin) / (WKmax - WKmin), 0), 1)
    if (fi != 1) {
        impact <- c(himpact[HoW[1:(fi_xtime - 1)]], rep_len(himpact, length(xtime) - fi + 1))
    } else {
        impact <- rep_len(himpact, length(xtime))
    }
    return(impact)
}



# get.statesmulti: calculate socio-economic states for avergae of given day and hours  ------------------------------------------

get.statesmulti <- function(rr, freq = 168, ltimestamps = list(1), nobs = length(rr) + freq, lagsADAM, hADAM) {
    # rr: input ts
    # freq: frequency of substructure, e.g. for hourly data with daily slicing freq=24, with weekly slicing freq=168, for have hourly data freq=48 and 336 etc.
    # ltimestamps: A list of timestamps per week in 0,...,freq, to specify the aggregation(~mean) on the input and output level:
    # aggregation on the input level specified by list elements (~ vectors); aggregation on the output level specified by length of list (~ number of list elements);
    # - Example 1: Average input time series for HoD 11, 12, 13 on Mon, Tue, Wed.
    #   -> ltimestamps = list(c(11:13, 11:13+24, 11:13+48))
    # **Note**: The time series should start on a Monday at 00:00 ~ HoW = 0.
    # Otherwise, add monid ~ index of first occurance of Mon 0:00 in ts to ltimestamps.
    # nobs: number of output observations, usaully length(rr)+ forecasting horizon
    # lagsADAM: 1 for ANN model, period of seasonal component otherwise, e.g. lagsADAM=52, does not allow multiple frequencies
    # hADAM: in-sample TMAE-loss horizon

    get.seqOut <- function(vtimestamp) sort(do.call(c, lapply(vtimestamp, function(x) seq(x, nobs, by = freq))))

    get.inrr <- function(vtimestamp) {
        seqOut <- sort(do.call(c, lapply(vtimestamp, function(x) seq(x, nobs, by = freq))))
        seqAv <- seqOut[seqOut < length(rr)]
        mat <- matrix(c(rr[seqAv], rep(NA, (ceiling(length(seqAv) / length(vtimestamp)) * length(vtimestamp) - length(seqAv)))), byrow = TRUE, nrow = ceiling(length(seqAv) / length(vtimestamp)), ncol = length(vtimestamp))
        rowMeans(mat, na.rm = TRUE)
    }

    seqOut <- lapply(ltimestamps, function(x) get.seqOut(x))
    inrr <- lapply(ltimestamps, function(x) get.inrr(x))
    inlen <- lapply(ltimestamps, function(x) length(x))

    modelADAM <- numeric(9)
    get.statesid <- function(inrr, seqOut, inlen) {
        qmod <- smooth::adam(inrr,
            model = "ANA", lags = lagsADAM, orders = list(ma = c(0, 0), ar = c(0, 0)), initial = "backcasting",
            distribution = "dlaplace", loss = "TMAE", h = hADAM
        )
        modelADAM <<- qmod
        print(qmod)
        if (length(seqOut) %% inlen != 0) {
            seqOut <- head(seqOut, -(length(seqOut) %% inlen))
        }
        htmp <- length(seqOut) / inlen - length(inrr)
        states <- matrix(, nrow = nobs, ncol = 3)
        states[seqOut, 1] <- c(rep(fitted(qmod), each = inlen), rep(forecast(qmod, h = htmp)$mean, each = inlen))
        states[seqOut, 2] <- c(rep(tail(qmod$states[, 1], length(inrr)), each = inlen), rep(tail(qmod$states[, 1], 1), (htmp) * inlen))

        if (dim(qmod$states)[2] >= 2) {
            if (lagsADAM < htmp) {
                id.fill <- htmp - lagsADAM
                shat <- tail(qmod$states[, 2], lagsADAM)
                shat <- c(shat, head(shat, id.fill))
            }
            if (lagsADAM == htmp) {
                shat <- tail(qmod$states[, 2], lagsADAM)
            }
            if (lagsADAM > htmp) {
                shat <- head(tail(qmod$states[, 2], lagsADAM), htmp)
            }
            states[seqOut, 3] <- c(rep(tail(qmod$states[, 2], length(inrr)), each = inlen), rep(shat, each = inlen))
        } else {
            states[, 3] <- 0
        }

        get.state <- function(z) na.approx(states[, z], na.rm = FALSE, yleft = states[head(which(!is.na(states[, 1])), 1), z], yright = states[tail(which(!is.na(states[, z])), 1), z])
        tmp <- sapply(2:3, get.state)
        dimnames(tmp) <- list(NULL, paste(c("level", "season"), sep = "_"))
        print(qmod$B["alpha"])
        return(list(tmp = tmp, alp = qmod$B))
    }

    STATES <- list()
    ALP <- list()
    for (i.r in seq_along(inlen)) {
        STATES[[i.r]] <- get.statesid(inrr = inrr[[i.r]], seqOut = seqOut[[i.r]], inlen = inlen[[i.r]])[[1]]
        ALP[[i.r]] <- get.statesid(inrr = inrr[[i.r]], seqOut = seqOut[[i.r]], inlen = inlen[[i.r]])[[2]]
    }

    return(list(as.data.frame(Reduce(`+`, STATES) / length(inlen)), ALP))
}



# forecast_AR_GAM: gam + AR model to forecast load -------------------------------------

forecast_AR_GAM <- function(Y, HoD, HoW, HoY, HLDP, HLDfixed, HLDweekday,
                            rho, gamma, order.max, H, S, impact, temp_lt, ltimestamps,
                            freq, ADAM, lagsADAM, hADAM, khldfixed, khldweekday,
                            khldp, ktemp, kHoW, kHoD, kHoY, kHoYHoD, kHoYHoW,
                            kstates, knots, temp_in) {
    # extend Y with NA values for forecasting horizon
    Yext <- c(Y, rep(NA, H))

    # handle newly introduced holidays
    inclHLDfixed <- apply(HLDfixed, 2, sum) > 10
    inclHLDweekday <- apply(HLDweekday, 2, sum) > 10

    # construct deterministic term formula
    detterm <- paste0(
        paste0("s(", dimnames(HLDfixed)[[2]][inclHLDfixed], ", bs='cp', k=", khldfixed, ", by=impact, pc=0)", collapse = "+"), "+",
        paste0("s(", dimnames(HLDweekday)[[2]], ", bs='cp', k=", khldweekday, ", pc=0)", collapse = "+"), "+",
        paste0("s(", dimnames(HLDP)[[2]], ", bs='cp', k=", khldp, ", pc=0)", collapse = "+"),
        " + s(HoD, bs='cp', k=", kHoD, ")",
        "+ s(HoW, bs='cp', k=", kHoW, ")",
        "+ s(HoY, bs='cp', k=", kHoY, ")",
        "+ ti(HoY, HoD, bs='cp', k=c(", kHoYHoD[1], ",", kHoYHoD[2], "))",
        "+ ti(HoY, HoW, bs='cp', k=c(", kHoYHoW[1], ",", kHoYHoW[2], "))"
    )

    # construct data matrices for GAM
    if (temp_in) {
        Mgam <- cbind(
            data.frame(Y = Yext, temp_lt, impact, HoD, HoW, HoY),
            HLDP, HLDfixed, HLDweekday
        )
        dettempterm <- paste0(
            detterm,
            paste0(" + s(", dimnames(temp_lt)[[2]], ", k=", ktemp, ", bs='ps')", collapse = "")
        )
        form <- as.formula(paste0("Y ~ ", dettempterm))
    } else {
        Mgam <- cbind(
            data.frame(Y = Yext, impact, HoD, HoW, HoY),
            HLDP, HLDfixed, HLDweekday
        )
        form <- as.formula(paste0("Y ~ ", detterm))
    }

    # split training and test data
    Mgam_train <- Mgam[1:length(Y), ]
    Mgam_test <- tail(Mgam, H)

    # fit the first-stage GAM model

    moddet <- bam(
        form,
        select = TRUE, data = Mgam_train, discrete = TRUE, gamma = gamma,
        nthreads = 4, control = gam.control(trace = FALSE), rho = rho, knots = knots
    )

    # if ADAM model is used
    if (ADAM) {
        rr <- residuals(moddet)

        # ddjust timestamps for state estimation
        id.Mo <- which(HoW[1:(7 * S)] == 0)
        vltimestamps <- lapply(ltimestamps, function(x) x + id.Mo)

        # compute states
        STAT <- get.statesmulti(
            vltimestamps,
            freq = freq, rr = rr, nobs = dim(Mgam)[1],
            lagsADAM = lagsADAM, hADAM = hADAM
        )

        states <- STAT[[1]]
        stat_names <- dimnames(states)[[2]]


        # if only level is present
        if (all(states$season == 0)) {
            states <- states[, 1]
            stat_names <- stat_names[1]
        }
        alphas <- STAT[[2]]


        alphas <- STAT[[2]]

        # add state variables to GAM dataset
        Mgam_train[, stat_names] <- head(states, length(Y))
        Mgam_test[, stat_names] <- tail(states, H)

        # construct new GAM formula with state effects
        if (temp_in) {
            form <- as.formula(
                paste0(
                    "Y ~ ", dettempterm, " + ",
                    paste("s(HoD, by=", stat_names, " ,bs='cp', k=", kstates, ")",
                        collapse = "+"
                    )
                )
            )
        } else {
            form <- as.formula(
                paste0(
                    "Y ~ ", detterm, " + ",
                    paste("s(HoD, by=", stat_names, " ,bs='cp', k=", kstates, ")",
                        collapse = "+"
                    )
                )
            )
        }


        # fit final GAM model with state effects
        mod <- bam(
            form,
            select = TRUE, data = Mgam_train, discrete = TRUE,
            nthreads = 4, gamma = gamma, rho = rho, control = gam.control(trace = FALSE),
            knots = knots, coef = coef(moddet)
        )
    } else {
        mod <- moddet
        states <- NA
    }

    # extract intercept
    intercept <- coef(mod)[1]

    # estimate AR model on residuals

    modres <- ar(x = residuals(mod), aic = TRUE, order.max = order.max, method = "burg")

    # forecast using AR + GAM model

    forecasts <- predict(modres, n.ahead = H, se.fit = FALSE) +
        predict.bam(mod, newdata = Mgam_test, type = "response", se.fit = FALSE)

    # forecasts_GAM <- predict.bam(mod, newdata = Mgam_test, type = "response", se.fit = FALSE)

    # calculate fit
    # fit <- fitted(modres) + fitted(mod)

    output <- list(
        # "fit" = fit,
        "forecasts" = forecasts,
        # "forecastGAM" = forecasts_GAM,
        "modres" = modres,
        "GAMI" = moddet,
        "mod" = mod,
        "states" = states,
        "alphas" = alphas,
        "intercept" = intercept
    )
    return(output)
}

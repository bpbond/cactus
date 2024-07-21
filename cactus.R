
library(dplyr)
srdb_raw <- read.csv("srdb_data/srdb-data.csv")

srdb_raw %>%
    filter(Meas_method %in% c("IRGA", "Gas chromatography")) %>%
    select(Record_number, Study_midyear, Latitude, Longitude, Manipulation, Rs_annual, Rs_growingseason) %>%
    mutate(Latitude = round(Latitude, 2),
           Longitude = round(Longitude, 2),
           Year = floor(Study_midyear)) %>%
    filter(Manipulation == "None", Latitude > 0) %>%
    filter(!is.na(Rs_growingseason) | !is.na(Rs_annual)) %>%
    filter(!is.na(Latitude), !is.na(Longitude), !is.na(Study_midyear)) %>%
    as_tibble() ->
    srdb

coords_year <- srdb[c("Year", "Longitude", "Latitude")] %>% distinct()
coords <- srdb[c("Longitude", "Latitude")] %>% distinct()
message(nrow(coords), " distinct coordinate pairs")

# ---------- WorldClim climatology

message("Getting WorldClim data...")
library(geodata)
# geodata::worldclim_global is smart enough to *cache* the data--
# it will only download if the data don't exist at the path we give it
tavg <- worldclim_global("tavg", "5", "worldclim_data/")
prec <- worldclim_global("prec", "5", "worldclim_data/")

# Extract our points of interest. terra::extract() will handle making sure
# the coordinates get mapped to the correct grid cell(s) in the data
library(terra)
tavg_dat <- terra::extract(tavg, coords)
prec_dat <- terra::extract(prec, coords)

coords_wc <- coords
coords_wc$wcMAT <- rowMeans(tavg_dat[-1]) # monthly average
coords_wc$wcMAP <- rowSums(prec_dat[-1]) # monthly average

srdb <- left_join(srdb, coords_wc, by = c("Longitude", "Latitude"))

# ---------- ERA5 climate data

# This is expensive, so we cache the extracted ERA5 data
era5_dat_file <- paste("era5", nrow(coords_year), digest::digest(coords_year), sep = "_")
if(file.exists(era5_dat_file)) {
    message("Loading saved data ", era5_dat_file)
    era5_results <- readRDS(era5_dat_file)
} else {
    message("Getting ERA5 data; this is slow...")
    library(tidyr)
    get_era5 <- function(lon, lat, year) {
        f <- list.files("era5_data/", pattern = paste0("^", year), full.names = TRUE)
        if(length(f) == 0) return(NULL)

        era5 <- rast(f)
        dat <- terra::extract(era5, data.frame(lon, lat))
        dat <- pivot_longer(dat, -ID)
        dat$month <- rep(1:12, each = 4)
        dat$name[grep("2 metre temperature", dat$name)] <- "Tair"
        dat$name[grep("undefined", dat$name)] <- "LAI"
        dat$name[grep("Volumetric soil water", dat$name)] <- "VWC"
        dat %>%
            group_by(ID, month, name) %>%
            # sum the two LAI (high and low veg)
            summarise(value = sum(value), .groups = "drop")
    }

    # Do this row by row just to keep things manageable
    era5_results <- list()
    for(i in seq_len(nrow(coords_year))) {
        if(i %% 50 == 0) message(i, " of ", nrow(coords_year))
        era5_results[[i]] <- get_era5(coords_year$Longitude[i], coords_year$Latitude[i], coords_year$Year[i])
    }
    era5_results <- bind_rows(era5_results, .id = "row")
    saveRDS(era5_results, file = era5_dat_file)
}


# ---------- SPEI drought dataset

#coords <- coords[1:324,] # about 20%

# This is expensive, so we cache the extracted SPEI data
spei_dat_file <- paste("spei", nrow(coords), digest::digest(coords), sep = "_")
if(file.exists(spei_dat_file)) {
    message("Loading saved data ", spei_dat_file)
    spei_dat <- readRDS(spei_dat_file)
} else {
    message("Extracting SPEI data; this is slow...")
    # Documentation: https://spei.csic.es/database.html
    # Data downloaded 2024-07-17
    library(terra)
    spei <- rast("spei_data/spei01.nc")

    # Extract our points of interest. terra::extract() will handle making sure
    # the coordinates get mapped to the correct grid cell(s) in the data
    spei_dat <- terra::extract(spei, coords)
    saveRDS(spei_dat, file = spei_dat_file)
}

# Reshape data into a more manageable form
library(tidyr)
spei_monthly <- pivot_longer(spei_dat, -ID)
spei_monthly <- separate(spei_monthly, name, into = c("spei", "entry"), convert = TRUE)
# The SPEI data don't seem to provide 'time' explicitly in the netcdf, so
# compute it from the entries
spei_monthly$year <- ceiling(spei_monthly$entry / 12) + 1900
spei_monthly$month <- (spei_monthly$entry - 1) %% 12 + 1
spei_monthly$time <- with(spei_monthly, year + (month-1) / 12)

# Compute gsd - growing season drought
spei_monthly %>%
    arrange(ID, year, month) %>%
    group_by(ID, year) %>%
    summarise(gsd = mean(value[6:8]), .groups = "drop") ->
    gsd

# Merge back with coords
coords %>%
    mutate(ID = seq_len(nrow(coords))) %>%
    left_join(gsd, by = "ID") ->
    gsd

library(ggplot2)
p <- ggplot(gsd, aes(year, gsd, color = ID)) + geom_point()
print(p)

# For each lon/lat pair in the SRDB data, look up gsd0 (current year SPEI),
# gsd1 (last year), and gsd2 (etc)
get_spei <- function(lon, lat, yr) {
    gsd %>%
        filter(Longitude == lon, Latitude == lat) %>%
        arrange(year) -> x

    if(nrow(x) == 0) return(c(NA, NA, NA)) # not found

    y <- which(x$year == yr)
    x$gsd[c(y, y-1, y-2)]
}

srdb %>%
    select(Study_midyear, Longitude, Latitude) %>%
    mutate(year = floor(Study_midyear)) %>%
    distinct() ->
    srdb_spei

message("Looking up growing season SPEI (yr 0, -1, -2) for SRDB data...")
srdb_spei$gsd2 <- srdb_spei$gsd1 <- srdb_spei$gsd0 <- NA_real_
for(i in seq_len(nrow(srdb_spei))) {
    spei_i <- get_spei(srdb_spei$Longitude[i], srdb_spei$Latitude[i], srdb_spei$year[i])
    srdb_spei$gsd0[i] <- spei_i[1]
    srdb_spei$gsd1[i] <- spei_i[2]
    srdb_spei$gsd2[i] <- spei_i[3]
}

drought_cutoff <- -1

message("Merging back into SRDB data...")
srdb %>%
    left_join(srdb_spei, by = c("Study_midyear", "Longitude", "Latitude")) %>%

    mutate(
        # I am not sure about having an arbitrary 'drought cutoff'
        # What about a birch variable that is
        #   (gsd0 - gsd1) if gsd0 >=0, and 0 otherwise
        # In other words, the gsd difference from year 0 to year -1, but only
        # when gsd0 (the current year) is at least 0, i.e. normal
        birch = if_else(gsd0 > 0 & gsd1 < 0, gsd0 - gsd1, 0),

        # Logical variant
        birch_lgl = gsd0 > 0 & gsd1 < drought_cutoff,

        drought = gsd0 <= drought_cutoff,
        year1_drought = gsd1 <= drought_cutoff,
        year2_drought = gsd2 <= drought_cutoff,
        out_1year = !drought & year1_drought,
        out_2year = !drought & !year1_drought & year2_drought) ->
    srdb_final

# ================= Growing season tests

library(car)

srdb_final %>%
    filter(!is.na(Rs_growingseason),
           Rs_growingseason < 20) ->
    srdb_gs

# Null model - MAT, MAP, drought index as predictors
m0_gs <- lm( sqrt(Rs_growingseason) ~ wcMAT  + wcMAP + drought, data = srdb_gs)
print(summary(m0_gs))
car::Anova(m0_gs, type = "III")

# Birch effect model - does emergence from drought (no drought in the current
# year, but drought the previous year) have a significant effect?
m1_gs <- lm( sqrt(Rs_growingseason) ~ wcMAT  + wcMAP + drought + birch_lgl, data = srdb_gs)
print(summary(m1_gs))
car::Anova(m1_gs, type = "III")


# ================= Annual tests

srdb_final %>%
    filter(!is.na(Rs_annual),
           Rs_annual > 0,
           Rs_annual < 2000) ->
    srdb_an

# Null model - MAT, MAP, drought index as predictors
m0_an <- lm( sqrt(Rs_annual) ~ wcMAT  + wcMAP + drought, data = srdb_an)
print(summary(m0_an))
shapiro.test(residuals(m0_an))
car::Anova(m0_an, type = "III")

# Birch effect model - does emergence from drought (no drought in the current
# year, but drought the previous year) have a significant effect?
m1_an <- lm( sqrt(Rs_annual) ~ wcMAT  + wcMAP + drought + birch_lgl, data = srdb_an)
print(summary(m1_an))
car::Anova(m1_an, type = "III")


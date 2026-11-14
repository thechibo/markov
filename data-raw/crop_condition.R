# Run from the package root: Rscript data-raw/crop_condition.R
# To reuse a downloaded USDA archive, pass its path as the sole argument.
# Source: USDA NASS Quick Stats, https://www.nass.usda.gov/datasets
# The public crops archive is about 1 GB compressed. Read it in chunks.

years <- 1995:2022
crops <- c("CORN", "SOYBEANS")
states <- c("IA", "NE")
ratings <- c("very poor", "poor", "fair", "good", "excellent")
series <- expand.grid(rating = toupper(ratings), crop = crops,
                      stringsAsFactors = FALSE)
items <- paste0(series$crop, " - CONDITION, MEASURED IN PCT ", series$rating)

download_archive <- function() {
  site <- "https://www.nass.usda.gov"
  page <- readLines(paste0(site, "/datasets"), warn = FALSE)
  pattern <- "/sites/default/files/datasets/qs\\.crops_[0-9]{8}\\.txt\\.gz"
  paths <- regmatches(page, gregexpr(pattern, page))
  archive <- tail(sort(unique(unlist(paths))), 1L)
  if (!length(archive)) stop("USDA crops archive link not found.")

  url <- paste0(site, archive)
  message("Downloading ", url)
  destination <- file.path(tempdir(), basename(archive))
  old_options <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old_options))
  if (download.file(url, destination, mode = "wb") != 0L) {
    stop("USDA crops archive download failed.")
  }
  destination
}

read_ratings <- function(archive) {
  connection <- gzfile(archive, "rt")
  on.exit(close(connection))
  header <- readLines(connection, n = 1L)
  pattern <- paste0("\\t(", paste(items, collapse = "|"), ")\\t")
  matches <- list()
  repeat {
    # The bulk file contains NUL padding in some records.
    lines <- readLines(connection, n = 100000L, skipNul = TRUE)
    if (!length(lines)) break
    condition <- grep(" - CONDITION, MEASURED IN PCT ", lines,
                      value = TRUE, fixed = TRUE, useBytes = TRUE)
    hits <- grep(pattern, condition, value = TRUE, perl = TRUE, useBytes = TRUE)
    if (length(hits)) matches[[length(matches) + 1L]] <- hits
  }
  raw <- read.delim(text = paste(c(header, unlist(matches)), collapse = "\n"),
                    colClasses = "character", quote = "", check.names = FALSE)
  names(raw) <- tolower(names(raw))
  columns <- c("commodity_desc", "state_alpha", "year", "reference_period_desc",
               "agg_level_desc", "source_desc", "domain_desc", "short_desc",
               "value")
  if (!all(columns %in% names(raw))) stop("Unexpected Quick Stats schema.")
  raw <- subset(raw, agg_level_desc == "STATE" & source_desc == "SURVEY" &
                  domain_desc == "TOTAL" & state_alpha %in% states &
                  year %in% years & grepl("^WEEK #[0-9]+$", reference_period_desc))
  unique(data.frame(
    crop = raw$commodity_desc,
    state = raw$state_alpha,
    year = as.integer(raw$year),
    week = as.integer(sub("^WEEK #", "", raw$reference_period_desc)),
    rating = match(series$rating[match(raw$short_desc, items)], toupper(ratings)),
    percent = suppressWarnings(as.numeric(raw$value))
  ))
}

make_paths <- function(data, crop, state) {
  data <- data[data$crop == crop & data$state == state, ]
  if (anyDuplicated(data[c("year", "week", "rating")])) {
    stop("Conflicting ratings for ", crop, " in ", state, ".")
  }
  paths <- matrix(NA_integer_, length(years), 53L,
                  dimnames = list(years, sprintf("week_%02d", 1:53)))
  for (week in split(data, list(data$year, data$week), drop = TRUE)) {
    shares <- week$percent[match(seq_along(ratings), week$rating)]
    if (anyNA(shares) || any(shares < 0 | shares > 100) || sum(shares) == 0) next
    # Use half the reported total to allow for rounding. At an exact tie,
    # use the lower category.
    paths[as.character(week$year[1L]), week$week[1L]] <-
      which(cumsum(shares) >= sum(shares) / 2)[1L]
  }
  common <- which(colSums(is.na(paths)) == 0L)
  # A one-step transition must span exactly one week. Use the longest
  # consecutive run of common weeks, taking the earliest run if tied.
  runs <- split(common, common - seq_along(common))
  weeks <- if (length(runs)) runs[[which.max(lengths(runs))]] else integer()
  list(paths = paths[, weeks, drop = FALSE], common = length(common))
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
archive <- if (length(args)) args[[1L]] else download_archive()
message("Reading ", basename(archive))
data <- read_ratings(archive)

# Break a tie in length alphabetically by crop, then by state abbreviation.
comparison <- expand.grid(state = states, crop = crops, stringsAsFactors = FALSE)
candidates <- Map(function(crop, state) make_paths(data, crop, state),
                  comparison$crop, comparison$state)
comparison$common_weeks <- vapply(candidates, `[[`, integer(1L), "common")
comparison$trajectory_length <- vapply(candidates, function(x) ncol(x$paths),
                                       integer(1L))
print(comparison, row.names = FALSE)
winner <- which.max(comparison$trajectory_length)
crop_condition <- candidates[[winner]]$paths
if (ncol(crop_condition) < 2L) stop("No candidate has two consecutive common weeks.")

attr(crop_condition, "crop") <- comparison$crop[winner]
attr(crop_condition, "state") <- comparison$state[winner]
attr(crop_condition, "levels") <- ratings
attr(crop_condition, "source_archive") <- basename(archive)
message("Selected ", comparison$crop[winner], " / ", comparison$state[winner],
        ": ", nrow(crop_condition), " years by ", ncol(crop_condition), " weeks (",
        colnames(crop_condition)[1L], " to ", tail(colnames(crop_condition), 1L), ").")
saveRDS(crop_condition, file.path("data-raw", "crop_condition.rds"))



sky_snq_age_str_df <- function(connection) {
  #load libraries 
  library(DBI)
  library(odbc)
  
  # Hardcoded WRIA stream ID groups
  sky_ids <- c(41,42,38,30,43,21,39,37,34,22,36,35,19,31,32,29)
  snq_ids <- c(23, 24, 25, 26, 27, 28)
  
  query <- "
  SELECT
    f.FishSampleID,
    f.ScaleCardPageID,
    f.FinMarkAdSampleStatusID,
    fs.TotalAge,
    fs.FWAge,
    sc.SampleDate,
    ss.WRIAStreamID,
    cwtr.CWTReleaseInfoID,
    otolib.HatcheryMark AS OtolithLibrary_HatcheryMark
  FROM FishSample f
  JOIN FishScaleResult fs ON f.FishSampleID = fs.FishSampleID
  JOIN ScaleCardPage scp ON f.ScaleCardPageID = scp.ScaleCardPageID
  JOIN ScaleCard sc ON scp.ScaleCardID = sc.ScaleCardID
  JOIN Species s ON sc.SpeciesID = s.SpeciesID
  LEFT JOIN ScaleCardStreamSample ss ON sc.ScaleCardID = ss.ScaleCardID
  LEFT JOIN FishCWTResult cwtr ON f.FishSampleID = cwtr.FishSampleID
  LEFT JOIN FishOtolithResult fot ON f.FishSampleID = fot.FishSampleID
  LEFT JOIN OtolithLibrary otolib ON fot.OtolithLibraryID = otolib.OtolithLibraryID
  WHERE s.Name = 'Chinook'
  "
  
  chinook <- DBI::dbGetQuery(connection, query)
  
  # Shared filtering logic (now only uses `chinook` from outer scope)
  filter_group <- function(stream_ids) {
    fish <- chinook[
      !is.na(chinook$SampleDate) &
        !is.na(chinook$TotalAge) &
        chinook$FinMarkAdSampleStatusID == 2 &
        is.na(chinook$CWTReleaseInfoID) &
        is.na(chinook$OtolithLibrary_HatcheryMark) &
        chinook$WRIAStreamID %in% stream_ids &
        chinook$FWAge %in% c(1, 2),
    ]
    
    fish$ReturnYear <- as.integer(format(as.Date(fish$SampleDate), "%Y"))
    years <- sort(unique(fish$ReturnYear))
    
    summary_df <- data.frame(
      ReturnYear = years,
      TotalSamples = integer(length(years)),
      Age2_FingerlingType = integer(length(years)),
      Age3_FingerlingType = integer(length(years)),
      Age4_FingerlingType = integer(length(years)),
      Age5_FingerlingType = integer(length(years)),
      Age3_YearlingType = integer(length(years)),
      Age4_YearlingType = integer(length(years)),
      Age5_YearlingType = integer(length(years)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(years)) {
      yr <- years[i]
      fingerlings <- fish[fish$FWAge == 1 & fish$ReturnYear == yr, ]
      yearlings   <- fish[fish$FWAge == 2 & fish$ReturnYear == yr, ]
      
      summary_df$TotalSamples[i]         <- nrow(fingerlings) + nrow(yearlings)
      summary_df$Age2_FingerlingType[i]  <- sum(fingerlings$TotalAge == 2, na.rm = TRUE)
      summary_df$Age3_FingerlingType[i]  <- sum(fingerlings$TotalAge == 3, na.rm = TRUE)
      summary_df$Age4_FingerlingType[i]  <- sum(fingerlings$TotalAge == 4, na.rm = TRUE)
      summary_df$Age5_FingerlingType[i]  <- sum(fingerlings$TotalAge == 5, na.rm = TRUE)
      summary_df$Age3_YearlingType[i]    <- sum(yearlings$TotalAge == 3, na.rm = TRUE)
      summary_df$Age4_YearlingType[i]    <- sum(yearlings$TotalAge == 4, na.rm = TRUE)
      summary_df$Age5_YearlingType[i]    <- sum(yearlings$TotalAge == 5, na.rm = TRUE)
    }
    
    return(summary_df)
  }
  
  # Return both summaries
  list(
    sky_summary = filter_group(sky_ids),
    snq_summary = filter_group(snq_ids)
  )
}


# con <- dbConnect(odbc(),
#                  Driver = "SQL Server",
#                  Server = "SVR-WEBSQLCLUSTER",
#                  Database = "nrfss",
#                  Trusted_Connection = "Yes")
# source("functions/age_structure.R")
# 
# age_structure <- sky_snq_age_str_df(con)
# sky_age_structure <- age_structure$sky_summary
# snq_age_structure <- age_structure$snq_summary
# 
# dbDisconnect(con)
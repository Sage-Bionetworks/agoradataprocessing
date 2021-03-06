#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(synapser))
suppressPackageStartupMessages(library(assertr))
suppressPackageStartupMessages(library(agoradataprocessing))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-c", "--config"), type="character",
              help="Configuration file.", dest="config",
              metavar="config", default = "config-staging.json"),
  make_option(c("--store"), action="store_true", default=FALSE,
              dest="store", help="Store in Synapse [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

synLogin()

store <- FALSE

config <- jsonlite::fromJSON(opt$config)

processed_data <- agoradataprocessing::process_data(config = config)

#########################################
# Write out all data and store in Synapse
#########################################

processed_data$omics_scores %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$omicsOutputFileJSON)

processed_data$genetics_scores %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$geneticsOutputFileJSON)

processed_data$overall_scores %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$overallOutputFileJSON)

processed_data$teamInfo %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$teamInfoFileJSON)

processed_data$geneInfo %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$geneInfoFileJSON)

processed_data$diffExprData %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$diffExprFileJSON)

processed_data$network %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$networkOutputFileJSON)

processed_data$proteomics %>%
  jsonlite::toJSON(pretty=2, na=NULL) %>%
  readr::write_lines(config$proteomicsFileJSON)

processed_data$metabolomics %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$metabolomicsFileJSON)

processed_data$srm_data %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$srmDataOutputFileJSON)

processed_data$neuropath_data %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$neuropathOutputFileJSON)

processed_data$target_exp_validation_harmonized_data %>%
  jsonlite::toJSON(pretty=2, digits=NA) %>%
  readr::write_lines(config$targetExpValidationHarmonizedOutputFileJSON)

if (opt$store) {
  teamInfoJSON <- synStore(File(config$teamInfoFileJSON,
                                parent=config$outputFolderId),
                           used=c(config$teamInfoId,
                                  config$teamMemberInfoId),
                           forceVersion=FALSE)

  geneInfoFinalJSON <- synStore(File(config$geneInfoFileJSON,
                                     parent=config$outputFolderId),
                                used=c(config$diffExprDataId,
                                       config$igapDataId,
                                       config$eqtlDataId,
                                       config$medianExprDataId,
                                       config$brainExpressionDataId,
                                       config$targetListOrigId,
                                       config$druggabilityDataId),
                                forceVersion=FALSE)

  diffExprDataJSON <- synStore(File(config$diffExprFileJSON,
                                    parent=config$outputFolderId),
                               used=c(config$diffExprDataId,
                                      config$tissuesTableId,
                                      config$studiesTableId),
                               forceVersion=FALSE)

  networkDataJSON <- synStore(File(config$networkOutputFileJSON,
                                   parent=config$outputFolderId),
                              used=c(config$diffExprDataId,
                                     config$networkDataId),
                              forceVersion=FALSE)

  proteomicsDataJSON <- synStore(File(config$proteomicsFileJSON,
                                      parent=config$outputFolderId),
                                used=c(config$proteomicsDataId),
                                forceVersion=FALSE)

  metabolomicsDataJSON <- synStore(File(config$metabolomicsFileJSON,
                                        parent=config$outputFolderId),
                                   used=c(config$metabolomicsDataId),
                                   forceVersion=FALSE)

  omicsScoresDataJSON <- synStore(File(config$omicsOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$omicsScoresTableId),
                                        forceVersion=FALSE)

  geneticsScoresDataJSON <- synStore(File(config$geneticsOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$geneticsScoresTableId),
                                        forceVersion=FALSE)

  overallScoresDataJSON <- synStore(File(config$overallOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$geneticsScoresTableId,
                                               config$omicsScoresTableId),
                                        forceVersion=FALSE)

  srmDataJSON <- synStore(File(config$srmDataOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$srmDataId),
                                        forceVersion=FALSE)

  neuroPathDataJSON <- synStore(File(config$neuropathOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$neuropathDataId),
                                        forceVersion=FALSE)

  targetExpValidationHarmonizedDataJSON <- synStore(File(config$targetExpValidationHarmonizedOutputFileJSON,
                                        parent=config$outputFolderId),
                                        used=c(config$targetExpressionValidationHarmonizedId),
                                        forceVersion=FALSE)

  dataFiles <- c(diffExprDataJSON,
                 geneInfoFinalJSON,
                 teamInfoJSON,
                 networkDataJSON,
                 proteomicsDataJSON,
                 metabolomicsDataJSON,
                 omicsScoresDataJSON,
                 geneticsScoresDataJSON,
                 overallScoresDataJSON,
                 srmDataJSON,
                 neuroPathDataJSON,
                 targetExpValidationHarmonizedDataJSON)

  dataManifest <- purrr::map_df(.x=dataFiles,
                                .f=function(x) data.frame(id=x$properties$id,
                                                          version=x$properties$versionNumber))

  dataManifest %>% readr::write_csv(config$manifestFileCSV)

  dataManifestCsv <- synStore(File(config$manifestFileCSV,
                                   parent=config$outputFolderId),
                              used=dataFiles,
                              forceVersion=FALSE)
}

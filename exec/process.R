#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(synapser))
suppressPackageStartupMessages(library(assertr))
suppressPackageStartupMessages(library(agoradataprocessing))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-c", "--config"), type="character",
              help="Configuration file.", dest="config",
              metavar="config"),
  make_option(c("--store"), action="store_true", default=FALSE,
              dest="store", help="Store in Synapse [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

synLogin()

store <- FALSE

config <- jsonlite::fromJSON(opt$config)

studiesTable <- agoradataprocessing::get_studies_table(config$studiesTableId)

tissuesTable <- agoradataprocessing::get_tissues_table(config$tissuesTableId)

teamInfo <- agoradataprocessing::get_team_info(config$teamInfoId)
teamMemberInfo <- agoradataprocessing::get_team_member_info(config$teamMemberInfoId)
teamInfo <- agoradataprocessing::process_team(teamInfo, teamMemberInfo)

# Do things relying on mygene before loading synapser
# by running getMyGeneInfo.R
# now load them as processed here
geneInfoFinal <- agoradataprocessing::get_gene_info_table(config$mygeneInfoFileId)
orig_size <- nrow(geneInfoFinal)

####### Process target list
targetListOrig <- agoradataprocessing::get_target_list(config$targetListOrigId) %>%
  assertr::verify(Team %in% teamInfo$team)

geneInfoFinal <- agoradataprocessing::process_target_list(targetListOrig, geneInfoFinal) %>%
  assertr::verify(nrow(.) == orig_size)

# Add to gene info
igap <- agoradataprocessing::get_igap(config$igapDataId)
geneInfoFinal <- agoradataprocessing::process_igap(igap, geneInfoFinal) %>%
  assertr::verify(nrow(.) == orig_size)

eqtl <- agoradataprocessing::get_eqtl(config$eqtlDataId)
geneInfoFinal <- agoradataprocessing::process_eqtl(eqtl, geneInfoFinal) %>%
  assertr::verify(nrow(.) == orig_size)

brain_expression <- agoradataprocessing::get_brain_expression_data(config$brainExpressionDataId)

geneInfoFinal <- agoradataprocessing::process_brain_expression_data(brain_expression,
                                                                    geneInfoFinal,
                                                                    fdr_random_threshold=config$isChangedInADBrainThreshold) %>%
  assertr::verify(nrow(.) == orig_size)

median_expr_data <- agoradataprocessing::get_median_expression_data(config$medianExprDataId)
geneInfoFinal <- agoradataprocessing::process_median_expression_data(median_expr_data,
                                                                     geneInfoFinal) %>%
  assertr::verify(nrow(.) == orig_size)

# Druggability data
druggabilityData <- agoradataprocessing::get_druggability_data(config$druggabilityDataId)

geneInfoFinal <- agoradataprocessing::process_druggability_data(druggabilityData,
                                                                geneInfoFinal) %>%
  assertr::verify(nrow(.) == orig_size)


# Process gene expression (logfc and CI) data
# Drop existing gene symbol and add them later.
diffExprData <- agoradataprocessing::get_rnaseq_diff_expr_data(config$diffExprDataId,
                                                               models_to_keep = config$modelsToKeep) %>%
  assertr::verify(study %in% studiesTable$studyname)

diffExprData <- process_rnaseq_diff_expr_data(diffExprData, gene_info = geneInfoFinal,
                                              adj_p_value_threshold = config$adjPValThreshold,
                                              target_list = targetListOrig)


network <- agoradataprocessing::get_network_data(config$networkDataId)
network <- agoradataprocessing::process_network_data(network, geneInfoFinal)

# Data tests
stopifnot(config$geneInfoColumns %in% colnames(geneInfoFinal))
stopifnot(config$diffExprCols %in% colnames(diffExprData))
stopifnot(config$networkCols %in% colnames(network))

#########################################
# Write out all data and store in Synapse
#########################################

teamInfo %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$teamInfoFileJSON)

geneInfoFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$geneInfoFileJSON)

diffExprData %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$diffExprFileJSON)

network %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$networkOutputFileJSON)

if (opt$store) {
  teamInfoJSON <- synStore(File(config$teamInfoFileJSON,
                                parent=config$outputFolderId),
                           used=c(config$teamInfoId, config$teamMemberInfoId),
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

  dataFiles <- c(diffExprDataJSON,
                 geneInfoFinalJSON,
                 teamInfoJSON,
                 networkDataJSON)

  dataManifest <- purrr::map_df(.x=dataFiles,
                                .f=function(x) data.frame(id=x$properties$id,
                                                          version=x$properties$versionNumber))

  dataManifest %>% readr::write_csv(config$manifestFileCSV)

  dataManifestCsv <- synStore(File(config$manifestFileCSV,
                                   parent=config$outputFolderId),
                              used=dataFiles,
                              forceVersion=FALSE)
}

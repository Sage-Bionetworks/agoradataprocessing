#!/usr/bin/env Rscript

library(tidyverse)
library(synapser)
library(assertr)
library(agoradataprocessing)

synLogin()

store <- FALSE

config <- jsonlite::fromJSON("./config.json")

studiesTable <- agoradataprocessing::get_studies_table(config$studiesTableId)

tissuesTable <- agoradataprocessing::get_tissues_table(config$tissuesTableId)

teamInfo <- agoradataprocessing::get_team_info(config$teamInfoId)
teamMemberInfo <- agoradataprocessing::get_team_member_info(config$teamMemberInfoId)
teamInfo <- agoradataprocessing::process_team(teamInfo, teamMemberInfo)

####### Process target list
targetListOrig <- agoradataprocessing::get_target_list(config$targetListOrigId) %>%
  assertr::verify(Team %in% teamInfo$team)

nestedTargetListOrig <- agoradataprocessing::process_target_list(targetListOrig)

# Do things relying on mygene before loading synapser
# by running getMyGeneInfo.R
# now load them as processed here
geneInfoFinal <- agoradataprocessing::get_gene_info_table(config$mygeneInfoFileId)

# Add to gene info
igap <- agoradataprocessing::get_igap(config$igapDataId)
geneInfoFinal <- agoradataprocessing::process_igap(igap, geneInfoFinal)

eqtl <- agoradataprocessing::get_eqtl(config$eqtlDataId)
geneInfoFinal <- agoradataprocessing::process_eqtl(eqtl, geneInfoFinal)

brainExpression <-

geneInfoFinal <- geneInfoFinal %>% left_join(brainExpression)
geneInfoFinal$isChangedInADBrain[is.na(geneInfoFinal$isChangedInADBrain)] <- FALSE

medianExprData <- readr::read_csv(synGet(config$medianExprDataId)$path) %>%
  filter(ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
  dplyr::rename_all(tolower)

nestedMedianExprData <- medianExprData %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='medianexpression')

geneInfoFinal <- left_join(geneInfoFinal, nestedMedianExprData)

# Druggability data
druggabilityDataObj <- synapser::synGet(config$druggabilityDataId)
druggabilityData <- druggabilityDataObj$path %>%
  readr::read_csv() %>%
  select(-`HGNC Name`) %>%
  dplyr::rename(ensembl_gene_id=GeneID) %>%
  assertr::chain_start() %>%
  assertr::verify(ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
  assertr::verify(assertr::is_uniq(ensembl_gene_id)) %>%
  assertr::chain_end()

# Druggability data schemas
colnames(druggabilityData) <- gsub(" ", "_", tolower(colnames(druggabilityData)))

nestedDruggabilityData <- druggabilityData %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='druggability')

geneInfoFinal <- left_join(geneInfoFinal, nestedDruggabilityData)

targetListOrig <- targetListOrig %>%
  assertr::verify(ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id)

# Process gene expression (logfc and CI) data
# Drop existing gene symbol and add them later.
diffExprDataObj <- synapser::synGet(config$diffExprDataId)
diffExprData <- diffExprDataObj$path %>%
  readr::read_tsv() %>%
  dplyr::mutate(tmp=paste(Model, Comparison, Sex, sep=" ")) %>%
  dplyr::filter(tmp %in% modelsToKeep) %>%
  dplyr::select(-tmp) %>%
  dplyr::mutate(Study=stringr::str_replace(Study, "MAYO", "MayoRNAseq"),
                Study=stringr::str_replace(Study, "MSSM", "MSBB"),
                Sex=stringr::str_replace_all(Sex,
                                             c("ALL"="males and females",
                                               "FEMALE"="females only",
                                               "MALE"="males only")),
                Model=stringr::str_replace(Model, "\\.", " x "),
                Model=stringr::str_replace(Model, "Diagnosis", "AD Diagnosis"),
                logFC=round(logFC, digits = 3),
                fc=2**logFC
  ) %>%
  assertr::verify(Study %in% studiesTable$studyname) %>%
  dplyr::mutate(model=glue::glue("{model} ({sex})", model=Model, sex=Sex))

# filtering by p-value for a gene that is significant in at least one model or a
# nominated target
diffExprData <- diffExprData %>%
  dplyr::filter(adj.P.Val <= adjPValThreshold | (ensembl_gene_id %in% targetListOrig$ensembl_gene_id)) %>%
  dplyr::select(ensembl_gene_id) %>%
  distinct() %>%
  left_join(., diffExprData)

diffExprData <- diffExprData %>%
  dplyr::select(ensembl_gene_id, logFC, fc, CI.L, CI.R,
                adj.P.Val, Tissue, Study, model)

colnames(diffExprData) <- gsub("\\.", "_", tolower(colnames(diffExprData)))

diffExprDataFinal <- left_join(diffExprData,
                               geneInfoFinal %>%
                                 dplyr::select(ensembl_gene_id, hgnc_symbol)) %>%
  filter(!is.na(hgnc_symbol)) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, everything())

network <- readr::read_csv(synGet(config$networkDataId)$path)

network <- network %>%
  filter(geneA_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id,
         geneB_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id)

network <- network %>%
  dplyr::select(geneA_ensembl_gene_id, geneB_ensembl_gene_id, brainRegion) %>%
  left_join(., geneInfoFinal %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% distinct(),
            by=c("geneA_ensembl_gene_id"="ensembl_gene_id")) %>%
  dplyr::rename(geneA_external_gene_name=hgnc_symbol) %>%
  left_join(., geneInfoFinal %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% distinct(),
            by=c("geneB_ensembl_gene_id"="ensembl_gene_id")) %>%
  dplyr::rename(geneB_external_gene_name=hgnc_symbol) %>%
  assertr::chain_start() %>%
  assertr::verify(geneA_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
  assertr::verify(geneB_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
  assertr::verify(geneA_external_gene_name %in% geneInfoFinal$hgnc_symbol) %>%
  assertr::verify(geneB_external_gene_name %in% geneInfoFinal$hgnc_symbol) %>%
  assertr::chain_end()

# Data tests
stopifnot(geneInfoColumns %in% colnames(geneInfoFinal))
stopifnot(diffExprCols %in% colnames(diffExprDataFinal))
stopifnot(networkCols %in% colnames(network))

#########################################
# Write out all data and store in Synapse
#########################################

teamInfo %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$teamInfoFileJSON)

geneInfoFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$geneInfoFileJSON)

diffExprDataFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$diffExprFileJSON)

network %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(config$networkOutputFileJSON)

if (store) {
  teamInfoJSON <- synStore(File(config$teamInfoFileJSON,
                                parent=config$outputFolderId),
                           used=c(teamInfoId, teamMemberInfoId),
                           forceVersion=FALSE)

  geneInfoFinalJSON <- synStore(File(config$geneInfoFileJSON,
                                     parent=config$outputFolderId),
                                used=c(diffExprDataId,
                                       igapDataId,
                                       eqtlDataId,
                                       medianExprDataId,
                                       brainExpressionDataId,
                                       targetListOrigId,
                                       druggabilityDataId),
                                forceVersion=FALSE)

  diffExprDataJSON <- synStore(File(config$diffExprFileJSON,
                                    parent=config$outputFolderId),
                               used=c(diffExprDataId,
                                      tissuesTableId,
                                      studiesTableId),
                               forceVersion=FALSE)

  networkDataJSON <- synStore(File(config$networkOutputFileJSON,
                                   parent=config$outputFolderId),
                              used=c(diffExprDataId, networkDataId),
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

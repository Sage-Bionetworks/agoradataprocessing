library(tidyverse)
library(synapser)
library(assertr)

synLogin()

outputFolderId <- 'syn12177492'

studiesTableId <- "syn11363298"
tissuesTableId <- "syn12180244"

mygeneInfoFileId <- 'syn15666826'
diffExprDataId <- 'syn11180450'
targetListOrigId <- "syn12540368"
networkDataId <- "syn11685347"
igapDataId <- "syn12514826"
eqtlDataId <- "syn12514912"
medianExprDataId <- "syn12514804"
brainExpressionDataId <- "syn11914808"

# Team and team member info as csv files
teamInfoId <- "syn12615624"
teamMemberInfoId <- "syn12615633"

networkOutputFileJSON <- "./network.json"
diffExprFileJSON <- "./rnaseq_differential_expression.json"
teamInfoFileJSON <- "./team_info.json"
geneInfoFileJSON <- "./gene_info.json"

# Required columns for output testing
geneInfoColumns <- c("ensembl_gene_id", "name",
                     "summary",
                     "hgnc_symbol", "alias",
                     "go.MF", "nominatedtarget",
                     "nominations", "isIGAP",
                     "haseqtl", "isChangedInADBrain",
                     "medianexpression")

diffExprCols <- c("ensembl_gene_id", "logfc", "fc", "ci_l", "ci_r",
                  "adj_p_val", "tissue", "study", "model")

networkCols <- c("geneA_ensembl_gene_id", "geneB_ensembl_gene_id",
                 "geneA_external_gene_name", "geneB_external_gene_name",
                 "brainRegion")

# threshold for if fdr determines if changed in AD brain
isChangedInADBrainThreshold <- 0.05

# Threshold for keeping differentially expressed genes
adjPValThreshold <- 0.001

studiesTable <- synapser::synTableQuery(glue::glue("select * from {studiesTableId}")) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::select(-ROW_ID, -ROW_VERSION) %>%
  dplyr::rename(studyid=Study, studyname=StudyName)

tissuesTable <- synapser::synTableQuery(glue::glue("select * from {tissuesTableId}")) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::select(-ROW_ID, -ROW_VERSION)

tissuesLookup <- setNames(as.character(tissuesTable$name),
                          tissuesTable$shortname)

modelsToKeep <- c("Diagnosis AD-CONTROL ALL",
                  "Diagnosis.AOD AD-CONTROL ALL",
                  "Diagnosis.Sex AD-CONTROL FEMALE",
                  "Diagnosis.Sex AD-CONTROL MALE")

teamMemberInfo <- synGet(teamMemberInfoId)$path %>%
  readr::read_csv() %>%
  dplyr::arrange(name) %>%
  dplyr::group_by(team) %>%
  tidyr::nest(.key="members")

teamInfo <- synGet(teamInfoId)$path %>%
  readr::read_csv() %>%
  dplyr::left_join(teamMemberInfo)

####### Process target list
targetListOrig <- synGet(targetListOrigId)$path %>%
  readr::read_csv() %>%
  dplyr::select(-hgnc_symbol) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::mutate(data_synapseid=stringr::str_split(data_synapseid, ",")) %>%
  assertr::verify(team %in% teamInfo$team)

nestedTargetListOrig <- targetListOrig %>%
  dplyr::group_by(ensembl_gene_id) %>%
  tidyr::nest(.key='nominatedtarget') %>%
  dplyr::mutate(nominations=purrr::map_int(nominatedtarget, nrow))

# Do things relying on mygene before loading synapser
# by running getMyGeneInfo.R
# now load them as processed there.
mygeneInfoFileObj <- synGet(mygeneInfoFileId)
load(mygeneInfoFileObj$path)

# Merge in target info to gene info
geneInfoFinal <- geneTableMerged %>%
  dplyr::rename(hgnc_symbol=symbol) %>%
  dplyr::filter(!is.na(X_id)) %>%
  dplyr::select(-X_id, -X_score) %>%
  dplyr::full_join(., nestedTargetListOrig)

# Add to gene info
igap <- synGet(igapDataId)$path %>% readr::read_csv()

geneInfoFinal <- geneInfoFinal %>%
  dplyr::mutate(isIGAP=ensembl_gene_id %in% igap$ensembl_gene_id)

eqtl <- synGet(eqtlDataId)$path %>%
  readr::read_csv() %>%
  dplyr::select(ensembl_gene_id, haseqtl=hasEqtl)

geneInfoFinal <- geneInfoFinal %>%
  dplyr::left_join(eqtl)
geneInfoFinal$haseqtl[is.na(geneInfoFinal$haseqtl)] <- FALSE

brainExpression <- synapser::synGet(brainExpressionDataId)$path %>%
  readr::read_tsv() %>%
  dplyr::select(ensembl_gene_id, fdr.random) %>%
  dplyr::mutate(isChangedInADBrain = fdr.random <= isChangedInADBrainThreshold) %>%
  dplyr::select(-fdr.random)

geneInfoFinal <- geneInfoFinal %>% left_join(brainExpression)
geneInfoFinal$isChangedInADBrain[is.na(geneInfoFinal$isChangedInADBrain)] <- FALSE

medianExprData <- readr::read_csv(synGet(medianExprDataId)$path) %>%
  filter(ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
  dplyr::rename_all(tolower)

nestedMedianExprData <- medianExprData %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='medianexpression')

geneInfoFinal <- left_join(geneInfoFinal, nestedMedianExprData)

targetListOrig <- targetListOrig %>%
  assertr::verify(ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id)

# Process gene expression (logfc and CI) data
# Drop existing gene symbol and add them later.
diffExprDataObj <- synapser::synGet(diffExprDataId)
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
                # Tissue=stringr::str_replace_all(Tissue, tissuesLookup),
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

network <- readr::read_csv(synGet(networkDataId)$path)

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
  readr::write_lines(teamInfoFileJSON)

teamInfoJSON <- synStore(File(teamInfoFileJSON,
                              parent=outputFolderId),
                         used=c(teamInfoId, teamMemberInfoId),
                         forceVersion=FALSE)

geneInfoFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(geneInfoFileJSON)

geneInfoFinalJSON <- synStore(File(geneInfoFileJSON,
                                   parent=outputFolderId),
                              used=c(diffExprDataId, igapDataId,
                                     eqtlDataId, medianExprDataId,
                                     brainExpressionDataId,
                                     targetListOrigId),
                              forceVersion=FALSE)

diffExprDataFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(diffExprFileJSON)

diffExprDataJSON <- synStore(File(diffExprFileJSON,
                                  parent=outputFolderId),
                             used=c(diffExprDataId, tissuesTableId, studiesTableId),
                             forceVersion=FALSE)

network %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(networkOutputFileJSON)

networkDataJSON <- synStore(File(networkOutputFileJSON,
                                 parent=outputFolderId),
                            used=c(diffExprDataId, networkDataId),
                            forceVersion=FALSE)

dataFiles <- c(diffExprDataJSON, geneInfoFinalJSON,
               teamInfoJSON, networkDataJSON)

dataManifest <- purrr::map_df(.x=dataFiles,
                               .f=function(x) data.frame(id=x$properties$id,
                                                         version=x$properties$versionNumber))

dataManifest %>% readr::write_csv("./data_manifest.csv")

dataManifestCsv <- synStore(File("./data_manifest.csv",
                                 parent=outputFolderId),
                            used=dataFiles,
                            forceVersion=FALSE)

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
                  "adj_p_val", "tissue", "study", "model", "xxx")

networkCols <- c("geneA_ensembl_gene_id", "geneB_ensembl_gene_id",
                 "geneA_external_gene_name", "geneB_external_gene_name",
                 "brainRegion")

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
  arrange(name) %>%
  group_by(team) %>%
  nest(.key="members")

teamInfo <- synGet(teamInfoId)$path %>%
  readr::read_csv() %>%
  left_join(teamMemberInfo)

####### Process target list
targetListOrig <- synGet(targetListOrigId)$path %>%
  readr::read_csv() %>%
  dplyr::select(-hgnc_symbol) %>%
  dplyr::rename_all(tolower) %>%
  mutate(data_synapseid=stringr::str_split(data_synapseid, ",")) %>%
  assertr::verify(team %in% teamInfo$team)

nestedTargetListOrig <- targetListOrig %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='nominatedtarget') %>%
  mutate(nominations=purrr::map_int(nominatedtarget, nrow))

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
  full_join(., nestedTargetListOrig)

# Add to gene info
igap <- synGet(igapDataId)$path %>% readr::read_csv()

geneInfoFinal <- geneInfoFinal %>%
  mutate(isIGAP=ensembl_gene_id %in% igap$ensembl_gene_id)

eqtl <- synGet(eqtlDataId)$path %>% readr::read_csv() %>%
  dplyr::select(ensembl_gene_id, haseqtl=hasEqtl)

geneInfoFinal <- geneInfoFinal %>% left_join(eqtl)
geneInfoFinal$haseqtl[is.na(geneInfoFinal$haseqtl)] <- FALSE

brainExpression <- synapser::synGet(brainExpressionDataId)$path %>%
  readr::read_tsv() %>%
  dplyr::select(ensembl_gene_id, fdr.random) %>%
  dplyr::mutate(isChangedInADBrain = fdr.random <= 0.05) %>%
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
  filter(adj.P.Val <= 0.001 | (ensembl_gene_id %in% targetListOrig$ensembl_gene_id)) %>%
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

stopifnot(geneInfoColumns %in% colnames(geneInfoFinal))

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

stopifnot(diffExprCols %in% colnames(diffExprDataFinal))

diffExprDataFinal %>%
  jsonlite::toJSON(pretty=2) %>%
  readr::write_lines(diffExprFileJSON)

diffExprDataJSON <- synStore(File(diffExprFileJSON,
                                  parent=outputFolderId),
                             used=c(diffExprDataId, tissuesTableId, studiesTableId),
                             forceVersion=FALSE)

stopifnot(networkCols %in% colnames(network))

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

# networkNested <- network %>%
#   group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
#            geneA_external_gene_name, geneB_external_gene_name) %>%
#   nest(brainRegion) %>%
#   mutate(data=map(data, function(x) list(brainRegion=x$brainRegion))) %>%
#   assertr::verify(geneA_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
#   assertr::verify(geneB_ensembl_gene_id %in% geneInfoFinal$ensembl_gene_id) %>%
#   assertr::verify(geneA_external_gene_name %in% geneInfoFinal$hgnc_symbol) %>%
#   assertr::verify(geneB_external_gene_name %in% geneInfoFinal$hgnc_symbol)

# networkNested %>%
#   jsonlite::toJSON(pretty=2) %>%
#   readr::write_lines("network_nested.json")
#
# networkNestedDataJSON <- synStore(File("network_nested.json",
#                                        parent=outputFolderId),
#                                   used=c(diffExprDataId, networkDataId),
#                                   forceVersion=FALSE)

# medianExprData %>%
#   jsonlite::toJSON() %>%
#   readr::write_lines(paste0("./median_expression", ".json"))
#
# medianExprDataJSON <- synStore(File(paste0("./median_expression", ".json"),
#                                     parent=outputFolderId),
#                                used=c(diffExprDataId, medianExprDataId),
#                                forceVersion=FALSE)

# network2 <- network %>%
#   dplyr::group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
#                   geneA_external_gene_name, geneB_external_gene_name) %>%
#   dplyr::summarize(value=n_distinct(brainRegion)) %>%
#   #  regions=paste(unique(brainRegion), collapse=","))
#   dplyr::ungroup() %>%
#   dplyr::distinct()
#   # %>%
#   # dplyr::filter(value > 1)
#
# network2 <- network2 %>%
#   dplyr::left_join(network) %>%
#   dplyr::group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
#                   geneA_external_gene_name, geneB_external_gene_name, value) %>%
#   dplyr::summarize(title=paste(unique(brainRegion), collapse=",")) %>%
#   dplyr::ungroup()
#
# nodesForNetwork <- dplyr::bind_rows(network2 %>%
#                                       dplyr::select(gene=geneA_ensembl_gene_id,
#                                                     symbol=geneA_external_gene_name),
#                                     network2 %>%
#                                       dplyr::select(gene=geneB_ensembl_gene_id,
#                                                     symbol=geneB_external_gene_name)) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(label=ifelse(is.na(symbol), gene, symbol)) %>%
#   dplyr::select(id=gene, label) # %>%
#   # dplyr::left_join(scoreData %>% dplyr::select(id=ensembl.gene, Score)) %>%
#   # dplyr::mutate(title=sprintf("score: %0.1f", Score)) %>%
#   # dplyr::mutate(color=cut(Score, breaks=c(-Inf, 0, 2, 4, Inf),
#   #                         labels=c("red", "yellow", "orange", "green")))
#
# edgesForNetwork <- network2 %>%
#   dplyr::select(from=geneA_ensembl_gene_id, to=geneB_ensembl_gene_id, value, title)

# networkNodesOutputFile <- "./network_nodes.csv"
# networkEdgesOutputFile <- "./network_edges.csv"

# IMSRId <- 'syn11149859'
# IMSROutputFilePrefix <- "./IMSR_processed"
# IMSROutputFile <- "./IMSR_processed.feather"
# IMSR <- synGet(IMSRId)$path %>%
#   readr::read_csv() %>%
#   dplyr::select(`Strain ID`, `Strain/Stock`, Repository, `Gene Symbol`, URL) %>%
#   dplyr::mutate(`Gene Symbol`=toupper(`Gene Symbol`)) %>%
#   dplyr::mutate(`Strain/Stock`=stringr::str_c("<a href='", URL, "'>", `Strain/Stock`, "</a>"))
#
# feather::write_feather(IMSR, IMSROutputFile)
# fIMSR <- synStore(File(IMSROutputFile,
#                        parent=outputFolderId),
#                   used=c(IMSRId),
#                   forceVersion=FALSE)
#
#
#
# scoreDataId <- 'syn11688680'
# scoreData <- synGet(scoreDataId)$path %>%
#   readr::read_csv() %>%
#   dplyr::select(ensembl_gene_id=gene, hgnc_symbol=external_gene_name,
#                 ad_driver_score=adDriverScore)

# druggabilityDataId <- 'syn11420932'
# druggabilityDataOutputFile <- "./druggabilityData_IGAP.csv"

# druggabilityData <- synGet(druggabilityDataId)$path %>%
#   read_csv() %>%
#   select(-Center, -DDI_Interest, -some_number, -GENE_DESCRIPTION) %>%
#   rename(ensembl.gene=ENSG) %>%
#   mutate(status_assays=forcats::fct_recode(status_assays, unknown=""),
#          status_crystal_structure=forcats::fct_recode(status_crystal_structure, unknown=""),
#          status_pocket=forcats::fct_recode(status_pocket, unknown=""),
#          status_in_vivo_work=forcats::fct_recode(status_in_vivo_work, unknown=""),
#          status_known_ligands=forcats::fct_recode(status_known_ligands, unknown=""),
#          Lilly_DrugEBIlity_Consensus_Score=forcats::fct_explicit_na(as.character(Lilly_DrugEBIlity_Consensus_Score),
#                                                            na_level = "unk"),
#          `Lilly_GW_Druggability_Structure-based`=forcats::fct_explicit_na(as.character(`Lilly_GW_Druggability_Structure-based`),
#                                                                  na_level = "unk")) %>%
#   group_by(GENE_SYMBOL) %>%
#   slice(1) %>%
#   ungroup()
#
# druggabilityData <- druggabilityData %>%
#   select(GENE_SYMBOL, starts_with('status')) %>%
#   distinct() %>%
#   tidyr::gather(category, status, starts_with('status')) %>%
#   mutate(status_numeric=forcats::fct_recode(status, `0`="unknown", `0`='bad',
#                                    `1`="medium", `2`="good"),
#          status_numeric=levels(status_numeric)[as.numeric(status_numeric)]) %>%
#   select(GENE_SYMBOL, status_numeric) %>%
#   group_by(GENE_SYMBOL) %>%
#   summarize(sum_status=sum(as.numeric(status_numeric))) %>%
#   ungroup() %>%
#   right_join(druggabilityData, by=c('GENE_SYMBOL')) %>%
#   arrange(GENE_SYMBOL)
#
# write_csv(druggabilityData, druggabilityDataOutputFile)
#
# fDrugData <- synStore(File(druggabilityDataOutputFile,
#                            parent=outputFolderId),
#                       used=druggabilityDataId, forceVersion=FALSE)

####### Process target list
# targetListOrig <- synGet(targetListOrigId)$path %>%
#   read_csv()
#
# targetList <- targetListOrig %>%
#   select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)
#
# targetListDistinct <- targetList %>%
#   group_by(Gene, ensembl.gene) %>%
#   mutate(Centers=paste(Center, collapse=", "),
#          nominations=length(unique(Center))) %>%
#   ungroup() %>%
#   select(-Center) %>%
#   distinct() %>%
#   arrange(-nominations)
#
# targetManifest <- targetListDistinct %>%
#   left_join(druggabilityData, by=c('Gene'='GENE_SYMBOL',
#                                    'ensembl.gene'='ensembl.gene')) %>%
#   select(Gene,
#          Centers,
#          nominations) %>%
#   distinct() %>%
#   arrange(-nominations)
#
# gtexObj <- synGet('syn7542283')
#
# gtex <- fread(getFileLocation(gtexObj), data.table=FALSE) %>%
#   mutate(ensembl.gene=str_replace(Name, "\\..*", "")) %>%
#   # dplyr::filter(ensembl.gene %in% c(genesForNetwork$gene,
#   #                                   targetList$ensembl.gene)) %>%
#   select(ensembl.gene, hgnc_symbol=Description, starts_with('Brain'))
#
# gtex <- gtex %>%
#   tidyr::gather(tissue, medianFPKM, 3:ncol(gtex)) %>%
#   mutate(tissue=str_replace(tissue, "Brain - ", ""))
#
# medianGTEx <- median(gtex$medianFPKM)
# #
# nTargets <- targetListOrig %>% count(group)
#
# save(ISMR, nTargets, network, targetList, targetListOrig,
#      druggabilityData, targetManifest, targetListDistinct, genesForNetwork,
#      gg, geneFPKMLong, gtex, medianGTEx, diffExprData, geneDF, dForFilter,
#      file="./data2load.RData")
# synStore(File("./data2load.RData", parent="syn7525089"))
# diffExprData <- synapseClient::synGet("syn11180450") %>%
#   synapseClient::getFileLocation() %>%
#   readr::read_csv()

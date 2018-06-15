library(tidyverse)
library(mygene)
library(synapser)
library(biomaRt)

# library(assertr)

synLogin()

wotFolderId <- 'syn12177492'

geneExprDataId <- 'syn11180450'
IMSRId <- 'syn11149859'
geneFPKMId <- "syn5581268"
geneCovariatesId <- "syn5581227"
studiesTableId <- "syn11363298"
tissuesTableId <- "syn12180244"

scoreDataId <- 'syn11688680'
targetListOrigId <- "syn12540368"
networkDataId <- "syn11685347"

targetListOutputFile <- "./targetList.csv"
targetListDistinctOutputFile <- "./targetListDistinct.csv"
targetManifestOutputFile <- "./targetManifest.csv"
networkOutputFile <- "./network.csv"
networkNodesOutputFile <- "./network_nodes.csv"
networkEdgesOutputFile <- "./network_edges.csv"
networkOutputFileJson <- "./network.json"
networkOutputFilePrefix <- "./network"

targetListOutputFile <- "./targetList_IGAP.csv"
targetListDistinctOutputFile <- "./targetListDistinct_IGAP.csv"
targetManifestOutputFile <- "./targetManifest_GAP.csv"

fGeneExprDataOutputFilePrefix <- "./geneExprData"
IMSROutputFilePrefix <- "./IMSR_processed"
geneFPKMLongOutputFilePrefix <- "./geneFPKMLong"

fGeneExprDataOutputFile <- "./geneExprData.feather"
IMSROutputFile <- "./IMSR_processed.feather"
geneFPKMLongOutputFile <- "./geneFPKMLong.feather"

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

keep <- c("Diagnosis AD-CONTROL ALL",
          "Diagnosis.AOD AD-CONTROL ALL",
          "Diagnosis.Sex AD-CONTROL FEMALE",
          "Diagnosis.Sex AD-CONTROL MALE")

# Process gene expression (logfc and CI) data
# Drop existing gene symbol and add them later.
geneExprData <- synapser::synGet(geneExprDataId)$path %>%
  readr::read_tsv() %>%
  dplyr::mutate(tmp=paste(Model, Comparison, Sex, sep=" ")) %>%
  dplyr::filter(tmp %in% keep) %>%
  dplyr::select(-tmp) %>%
  dplyr::mutate(Study=stringr::str_replace(Study, "MAYO", "MayoRNAseq"),
                Study=stringr::str_replace(Study, "MSSM", "MSBB"),
                Sex=stringr::str_replace_all(Sex,
                                             c("ALL"="Males and females",
                                               "FEMALE"="Females only",
                                               "MALE"="Males only")),
                # Tissue=stringr::str_replace_all(Tissue, tissuesLookup),
                Model=stringr::str_replace(Model, "\\.", " x "),
                Model=stringr::str_replace(Model, "SourceDiagnosis", "Study-specific Diagnosis"),
                logFC=round(logFC, digits = 3)
  ) %>%
  assertr::verify(Study %in% studiesTable$studyname) %>%
  dplyr::mutate(model=glue::glue("{model} ({sex})", model=Model, sex=Sex))


# # Maybe for filtering by p-value
# geneExprData <- geneExprData %>%
#   filter(adj_p_val <= 0.001) %>%
#   dplyr::select(ensembl_gene_id) %>%
#   distinct() %>%
#   left_join(., geneExprData)

geneExprData <- geneExprData %>%
  dplyr::select(ensembl_gene_id, logFC, CI.L, CI.R, adj.P.Val, Tissue, Study, model=model)

colnames(geneExprData) <- gsub("\\.", "_", tolower(colnames(geneExprData)))

# For unique list of genes
geneTable <- geneExprData %>%
  dplyr::select(ensembl_gene_id) %>%
  distinct()

# Get current hgnc symbol from biomart
ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
fromBiomart <- biomaRt::getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"),
                              filters="ensembl_gene_id",
                              values=geneTable$ensembl_gene_id,
                              mart=ensembl)

# Use mygene.info to get name and summary
geneInfoRes <- mygene::getGenes(fromBiomart$ensembl_gene_id,
                                fields=c("symbol", "name", "summary", "type_of_gene"))

geneInfo <- geneInfoRes %>%
  as.data.frame() %>%
  as.tibble() %>%
  group_by(query) %>%
  top_n(1, X_score) %>%
  ungroup()

geneInfo <- left_join(fromBiomart, geneInfo,
                      by=c('ensembl_gene_id'='query',
                           'hgnc_symbol'='symbol')) %>%
  as.tibble() %>%
  distinct()

geneInfoFinal <- geneInfo %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, name, summary) %>%
  mutate(symbol=ifelse(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol)) %>%
  group_by(ensembl_gene_id) %>%
  summarize(n=n_distinct(hgnc_symbol)) %>%
  left_join(geneInfo, .) %>%
  arrange(desc(n), ensembl_gene_id) %>%
  filter(!is.na(X_id)) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, name, summary)

geneExprData <- left_join(geneExprData,
                          geneInfoFinal %>% dplyr::select(ensembl_gene_id, hgnc_symbol))

geneExprData %>%
  jsonlite::toJSON() %>%
  readr::write_lines(paste0(fGeneExprDataOutputFilePrefix, ".json"))

fGeneExprDataJSON <- synStore(File(paste0(fGeneExprDataOutputFilePrefix, ".json"),
                                   parent=wotFolderId),
                              used=c(geneExprDataId, tissuesTableId, studiesTableId),
                              forceVersion=FALSE)

fGeneExprDataCsv <- synStore(File(paste0(fGeneExprDataOutputFilePrefix, ".csv"),
                                  parent=wotFolderId),
                             used=c(geneExprDataId, tissuesTableId, studiesTableId),
                             forceVersion=FALSE)

# Add to gene info
igapDataId <- "syn12514826"
igap <- synGet(igapDataId)$path %>% readr::read_csv()

geneInfoFinal <- geneInfoFinal %>%
  mutate(isIGAP=ensembl_gene_id %in% igap$ensembl_gene_id)

eqtlDataId <- "syn12514912"
eqtl <- synGet(eqtlDataId)$path %>% readr::read_csv() %>%
  dplyr::select(ensembl_gene_id, haseqtl=hasEqtl)

geneInfoFinal <- geneInfoFinal %>% left_join(eqtl)
geneInfoFinal$haseqtl[is.na(geneInfoFinal$haseqtl)] <- FALSE

medianExprDataId <- "syn12514804"
medianExprData <- readr::read_csv(synGet(medianExprDataId)$path) %>%
  filter(ensembl_gene_id %in% geneExprData$ensembl_gene_id) %>%
  dplyr::rename_all(tolower)

nestedMedianExprData <- medianExprData %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='medianexpression')

geneInfoFinal <- left_join(geneInfoFinal, nestedMedianExprData)

####### Process target list
targetListOrig <- synGet(targetListOrigId)$path %>%
  readr::read_csv() %>%
  dplyr::select(-hgnc_symbol) %>%
  dplyr::rename_all(tolower) %>%
  mutate(data_synapseid=stringr::str_split(data_synapseid, ","))

nestedTargetListOrig <- targetListOrig %>%
  group_by(ensembl_gene_id) %>%
  nest(.key='nominatedtarget')

geneInfoFinal <- left_join(geneInfoFinal, nestedTargetListOrig)

geneInfoFinal %>%
  jsonlite::toJSON() %>%
  readr::write_lines(paste0("./gene_info", ".json"))

geneInfoFinalJSON <- synStore(File(paste0("./gene_info", ".json"),
                                   parent=wotFolderId),
                              used=c(geneExprDataId, igapDataId, eqtlDataId, medianExprDataId,
                                     targetListOrigId),
                              forceVersion=FALSE)


network <- readr::read_csv(synGet(networkDataId)$path)

network <- network %>%
  dplyr::select(geneA_ensembl_gene_id, geneB_ensembl_gene_id, brainRegion) %>%
  left_join(., geneExprData %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% distinct(),
            by=c("geneA_ensembl_gene_id"="ensembl_gene_id")) %>%
  dplyr::rename(geneA_external_gene_name=hgnc_symbol) %>%
  left_join(., geneExprData %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% distinct(),
            by=c("geneB_ensembl_gene_id"="ensembl_gene_id")) %>%
  dplyr::rename(geneB_external_gene_name=hgnc_symbol)

network <- network %>%
  mutate(geneA_external_gene_name=ifelse(is.na(geneA_external_gene_name),
                                         geneA_ensembl_gene_id, geneA_external_gene_name),
         geneB_external_gene_name=ifelse(is.na(geneB_external_gene_name),
                                         geneB_ensembl_gene_id, geneB_external_gene_name))


network <- network %>%
  filter(geneA_ensembl_gene_id %in% geneExprData$ensembl_gene_id,
         geneB_ensembl_gene_id %in% geneExprData$ensembl_gene_id)

network %>%
  jsonlite::toJSON() %>%
  readr::write_lines(paste0(networkOutputFilePrefix, ".json"))

networkDataJSON <- synStore(File(paste0(networkOutputFilePrefix, ".json"),
                                 parent=wotFolderId),
                            used=c(geneExprDataId, networkDataId),
                            forceVersion=FALSE)

networkDataCsv <- synStore(File(paste0(networkOutputFilePrefix, ".csv"),
                                parent=wotFolderId),
                           used=c(geneExprDataId, networkDataId),
                           forceVersion=FALSE)


medianExprData %>%
  jsonlite::toJSON() %>%
  readr::write_lines(paste0("./median_expression", ".json"))

medianExprDataJSON <- synStore(File(paste0("./median_expression", ".json"),
                                    parent=wotFolderId),
                               used=c(geneExprDataId, medianExprDataId),
                               forceVersion=FALSE)

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

# IMSR <- synGet(IMSRId)$path %>%
#   readr::read_csv() %>%
#   dplyr::select(`Strain ID`, `Strain/Stock`, Repository, `Gene Symbol`, URL) %>%
#   dplyr::mutate(`Gene Symbol`=toupper(`Gene Symbol`)) %>%
#   dplyr::mutate(`Strain/Stock`=stringr::str_c("<a href='", URL, "'>", `Strain/Stock`, "</a>"))
#
# feather::write_feather(IMSR, IMSROutputFile)
# fIMSR <- synStore(File(IMSROutputFile,
#                        parent=wotFolderId),
#                   used=c(IMSRId),
#                   forceVersion=FALSE)
#
#
#
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
#                            parent=wotFolderId),
#                       used=druggabilityDataId, forceVersion=FALSE)

####### Process target list
# targetListOrig <- synGet(targetListOrigId)$path %>%
#   read_csv()
#
# targetList <- targetListOrig %>%
#   select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)
#
# write_csv(targetList, targetListOutputFile)
#
# fTargetList <- synStore(File(targetListOutputFile,
#                              parent=wotFolderId),
#                         used=targetListOrigId, forceVersion=FALSE)
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
# write_csv(targetList, targetListDistinctOutputFile)
#
# fTargetListDistinct <- synStore(File(targetListDistinctOutputFile,
#                              parent=wotFolderId),
#                         used=fTargetList, forceVersion=FALSE)
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
# write_csv(targetManifest, targetManifestOutputFile)
#
# fTargetManifest <- synStore(File(targetManifestOutputFile,
#                                      parent=wotFolderId),
#                                 used=c(fTargetListDistinct, fDrugData),
#                             forceVersion=FALSE)


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
#      gg, geneFPKMLong, gtex, medianGTEx, geneExprData, geneDF, dForFilter,
#      file="./data2load.RData")
# synStore(File("./data2load.RData", parent="syn7525089"))
# diffExprData <- synapseClient::synGet("syn11180450") %>%
#   synapseClient::getFileLocation() %>%
#   readr::read_csv()

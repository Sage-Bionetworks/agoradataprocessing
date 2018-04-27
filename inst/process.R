


library(tidyverse)
library(synapser)
library(feather)

synLogin()

wotFolderId <- 'syn7525089'

geneExprDataId <- 'syn11180450'
IMSRId <- 'syn11149859'
geneFPKMId <- "syn5581268"
geneCovariatesId <- "syn5581227"


scoreDataId <- 'syn11688680'
targetListOrigId <- "syn8656625"
networkDataId <- "syn11685347"

targetListOutputFile <- "./targetList.csv"
targetListDistinctOutputFile <- "./targetListDistinct.csv"
targetManifestOutputFile <- "./targetManifest.csv"
networkOutputFile <- "./network.csv"
networkNodesOutputFile <- "./network_nodes.csv"
networkEdgesOutputFile <- "./network_edges.csv"
networkOutputFileJson <- "./network.json"

targetListOutputFile <- "./targetList_IGAP.csv"
targetListDistinctOutputFile <- "./targetListDistinct_IGAP.csv"
targetManifestOutputFile <- "./targetManifest_GAP.csv"

fGeneExprDataOutputFilePrefix <- "./geneExprData"
IMSROutputFilePrefix <- "./IMSR_processed"
geneFPKMLongOutputFilePrefix <- "./geneFPKMLong"

fGeneExprDataOutputFile <- "./geneExprData.feather"
IMSROutputFile <- "./IMSR_processed.feather"
geneFPKMLongOutputFile <- "./geneFPKMLong.feather"

### Process gene expression (logfc and CI) data
geneExprData <- synGet(geneExprDataId)$path %>%
  read_tsv() %>%
  filter(stringr::str_detect(Model, "Diagnosis")) %>%
  mutate(hgnc_symbol=ifelse(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol),
         Model=stringr::str_replace(Model, "\\.", " x "),
         Model=stringr::str_replace(Model, "SourceDiagnosis", "Study-specific Diagnosis"),
         neg.log10.adj.P.Val=-log10(adj.P.Val)) %>%
  tidyr::unite("tissue_study", Tissue, Study, sep=", ", remove=FALSE) %>%
  tidyr::unite("comparison_model_sex", Comparison, Model, Sex, sep=", ", remove=FALSE) %>%
  mutate(tissue_study_pretty=glue::glue("{tissue} ({study})", tissue=Tissue, study=Study),
         comparison_model_sex_pretty=glue::glue("{comparison} {model} ({sex})", comparison=Comparison, model=Model, sex=Sex))

colnames(geneExprData) <- gsub("\\.", "_", tolower(colnames(geneExprData)))

feather::write_feather(geneExprData, paste0(fGeneExprDataOutputFilePrefix, ".feather"))
readr::write_csv(geneExprData, paste0(fGeneExprDataOutputFilePrefix, ".csv"))

fGeneExprDataFeather <- synStore(File(paste0(fGeneExprDataOutputFilePrefix, ".feather"),
                                      parent=wotFolderId),
                                 used=c(geneExprDataId),
                                 forceVersion=FALSE)

fGeneExprDataCsv <- synStore(File(paste0(fGeneExprDataOutputFilePrefix, ".csv"),
                                  parent=wotFolderId),
                             used=c(geneExprDataId),
                             forceVersion=FALSE)

IMSR <- synGet(IMSRId)$path %>%
  readr::read_csv() %>%
  dplyr::select(`Strain ID`, `Strain/Stock`, Repository, `Gene Symbol`, URL) %>%
  mutate(`Gene Symbol`=toupper(`Gene Symbol`)) %>%
  mutate(`Strain/Stock`=stringr::str_c("<a href='", URL, "'>", `Strain/Stock`, "</a>"))

write_feather(IMSR, IMSROutputFile)
fIMSR <- synStore(File(IMSROutputFile,
                       parent=wotFolderId),
                  used=c(IMSRId),
                  forceVersion=FALSE)



scoreData <- synGet(scoreDataId)$path %>%
  read_csv() %>%
  rename(ensembl.gene=gene, Score=adDriverScore, Gene=external_gene_name)

####### Process target list
targetListOrig <- synGet(targetListOrigId)$path %>%
  read_csv()

targetList <- targetListOrig %>%
  select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)

write_csv(targetList, targetListOutputFile)

fTargetList <- synStore(File(targetListOutputFile,
                             parent=wotFolderId),
                        used=targetListOrigId, forceVersion=FALSE)

targetListDistinct <- targetList %>%
  group_by(Gene, ensembl.gene) %>%
  mutate(Centers=paste(Center, collapse=", "),
         nominations=length(unique(Center))) %>%
  ungroup() %>%
  select(-Center) %>%
  distinct() %>%
  arrange(-nominations)

write_csv(targetList, targetListDistinctOutputFile)

fTargetListDistinct <- synStore(File(targetListDistinctOutputFile,
                                     parent=wotFolderId),
                                used=fTargetList, forceVersion=FALSE)

targetManifest <- scoreData %>%
  filter(Gene %in% targetListDistinct$Gene) %>%
  select(Gene, Score) %>%
  mutate(Score=round(Score, 1)) %>%
  arrange(-Score)

write_csv(targetManifest, targetManifestOutputFile)

fTargetManifest <- synStore(File(targetManifestOutputFile,
                                 parent=wotFolderId),
                            used=c(fTargetListDistinct, scoreDataId),
                            forceVersion=FALSE)


network <- readr::read_csv(synGet(networkDataId)$path)

network2 <- network %>%
  group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
           geneA_external_gene_name, geneB_external_gene_name) %>%
  summarize(value=n_distinct(brainRegion)) %>%
  #  regions=paste(unique(brainRegion), collapse=","))
  ungroup() %>%
  distinct() %>%
  filter(value > 1)

network2 <- network2 %>% left_join(network) %>%
  group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
           geneA_external_gene_name, geneB_external_gene_name, value) %>%
  summarize(title=paste(unique(brainRegion), collapse=",")) %>%
  ungroup()

nodesForNetwork <- dplyr::bind_rows(network2 %>%
                                      select(gene=geneA_ensembl_gene_id,
                                             symbol=geneA_external_gene_name),
                                    network2 %>%
                                      select(gene=geneB_ensembl_gene_id,
                                             symbol=geneB_external_gene_name)) %>%
  distinct() %>%
  mutate(label=ifelse(is.na(symbol), gene, symbol)) %>%
  select(id=gene, label) %>%
  left_join(scoreData %>% select(id=ensembl.gene, Score)) %>%
  mutate(title=sprintf("score: %0.1f", Score)) %>%
  mutate(color=cut(Score, breaks=c(-Inf, 0, 2, 4, Inf),
                   labels=c("red", "yellow", "orange", "green")))

edgesForNetwork <- network2 %>%
  select(from=geneA_ensembl_gene_id, to=geneB_ensembl_gene_id, value, title)

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

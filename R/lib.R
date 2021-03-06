#' Load study information from Synapse
#'
#' @export
get_studies_table <- function(id) {
  synapser::synTableQuery(glue::glue("select * from {id}"))$asDataFrame() %>%
    tibble::as_tibble() %>%
    dplyr::select(-ROW_ID, -ROW_VERSION) %>%
    dplyr::rename(studyid=Study, studyname=StudyName)
}

#' Load tissue information from Synapse
#'
#' @export
get_tissues_table <- function(id) {
  synapser::synTableQuery(glue::glue("select * from {id}"))$asDataFrame() %>%
  tibble::as_tibble() %>%
  dplyr::select(-ROW_ID, -ROW_VERSION)
}

#' Load team info from Synapse
#'
#' @export
get_team_info <- function(id){
  synGet(id)$path %>% readr::read_csv()
}

#' Load team member info from Synapse
#'
#' @export
get_team_member_info <- function(id) {
    synGet(id)$path %>%
    readr::read_csv() %>%
    dplyr::arrange(name) %>%
    dplyr::group_by(team) %>%
    tidyr::nest(.key="members")
}

#' Process team info by joining with team member info.
#'
#' @export
process_team <- function(team_info, team_member_info) {
  team_info %>% dplyr::left_join(team_member_info)
}

#' Load gene info data from Synapse
#'
#' @export
get_gene_info_table <- function(id) {
  mygeneInfoFileObj <- synGet(id)
  load(mygeneInfoFileObj$path)

  # This rename/filter/select should be moved to getMyGeneInfo.R
  geneTableMerged %>%
    dplyr::rename(hgnc_symbol=symbol) %>%
    dplyr::rename(go_MF=go.MF) %>%
    dplyr::filter(!is.na(X_id)) %>%
    dplyr::select(-X_id, -X_score)
}

#' Load iGAP data from Synapse
#'
#' @export
get_igap <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_csv()
}

#' Add iGAP boolean to gene info table
#'
#' @export
process_igap <- function(data, gene_info) {
  gene_info %>%
    dplyr::mutate(isIGAP=ensembl_gene_id %in% data$ensembl_gene_id)
}

#' Load eQTL data from Synapse
#'
#' @export
get_eqtl <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_csv() %>%
    dplyr::rename_all(tolower) %>%
    dplyr::select(ensembl_gene_id, haseqtl)
}

#' Process eQTL data
#' @export
process_eqtl <- function(data, gene_info) {
  data %>%
    assertr::verify(assertr::has_all_names("ensembl_gene_id", "haseqtl")) %>%
    dplyr::left_join(gene_info, .) %>%
    dplyr::mutate(haseqtl=ifelse(is.na(haseqtl), FALSE, haseqtl))
}

#' Load target list from Synapse
#'
#' @export
get_target_list <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_csv()
}

#' Process the target list by getting specific columns,
#' splitting Synapse IDs for data used, and grouping
#' and nesting by Ensembl gene id into nominatedtarget column.
#'
#' @export
process_target_list <- function(data, gene_info) {
  data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename_all(.funs=stringr::str_replace_all, pattern=" ", replacement="_") %>%
    assertr::verify(assertr::has_all_names("ensembl_gene_id", "data_synapseid")) %>%
    dplyr::select(-hgnc_symbol) %>%
    dplyr::mutate(data_synapseid=stringr::str_split(data_synapseid, ",")) %>%
    dplyr::group_by(ensembl_gene_id) %>%
    tidyr::nest(.key='nominatedtarget') %>%
    dplyr::mutate(nominations=purrr::map_int(nominatedtarget, nrow)) %>%
    dplyr::left_join(gene_info, .)
}

#' Get brain expression data
#'
#' @export
get_brain_expression_data <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_tsv() %>%
    dplyr::select(ensembl_gene_id, fdr.random)
}

#' Process brain expression data
#' @export
process_brain_expression_data <- function(data, gene_info, fdr_random_threshold) {
  data %>%
    assertr::verify(assertr::has_all_names("ensembl_gene_id", "fdr.random")) %>%
    dplyr::mutate(isChangedInADBrain = fdr.random <= fdr_random_threshold) %>%
    dplyr::select(-fdr.random) %>%
    left_join(gene_info, .)
}

process_rna_change_in_the_brain <- function(gene_info, diff_exp) {

  rna_changed_data <- diff_exp %>%
      dplyr::select(ensembl_gene_id, adj_p_val) %>%
      assertr::verify(not_na(ensembl_gene_id)) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(minimum_adj_p_value = min(adj_p_val)) %>%
      dplyr::mutate(isAnyRNAChangedInADBrain = minimum_adj_p_value <= 0.05) %>%
      dplyr::select(- minimum_adj_p_value)
      
  data <- gene_info %>%
      assertr::verify(not_na(ensembl_gene_id)) %>%
      left_join(., rna_changed_data, by = c("ensembl_gene_id"))
}

#' Combine Proteomics Data with Gene Info
#' Relationship between esemble_gene_id and Cor_PVal is 1:n
#' 
#' @export 
process_proteomics_data <- function(gene_info, proteomics) {
proteomics %>% 
  dplyr::select(ensembl_gene_id, cor_pval) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(min_cor_pval = min(cor_pval)) %>%
  dplyr::mutate(isAnyProteinChangedInADBrain = min_cor_pval <= 0.05) %>%
  dplyr::select(-min_cor_pval) %>%
  dplyr::left_join(gene_info, ., by=c("ensembl_gene_id"))
}

#' Get median expression data
#'
#' @export
get_median_expression_data <- function(id) {
  synGet(id)$path %>%
    readr::read_csv() %>%
    dplyr::rename_all(tolower) %>%
    assertr::verify(assertr::has_all_names("ensembl_gene_id", "medianlogcpm", "tissue"))
}

#' Process median expression data
#'
#' @export
process_median_expression_data <- function(data, gene_info) {
  data %>%
    assertr::verify(assertr::has_all_names("ensembl_gene_id")) %>%
    group_by(ensembl_gene_id) %>%
    nest(.key='medianexpression') %>%
    left_join(gene_info, .)
}

#' Get druggability data from Synapse
#'
#' @export
get_druggability_data <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_csv() %>%
    dplyr::rename(ensembl_gene_id=GeneID) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename_all(.funs=stringr::str_replace_all, pattern=" ", replacement="_") %>%
    assertr::verify(assertr::is_uniq(ensembl_gene_id))
}

#' Process druggability data
#'
#' @export
process_druggability_data <- function(data, gene_info) {
  data %>%
    group_by(ensembl_gene_id) %>%
    nest(.key='druggability') %>%
    left_join(gene_info, .)
}

#' Get rnaseq differential expression data from Synapse
#'
#' @export
get_rnaseq_diff_expr_data <- function(id, models_to_keep) {
  synapser::synGet(id)$path %>%
    readr::read_tsv() %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename_all(.funs=stringr::str_replace_all, pattern=" ", replacement="_") %>%
    dplyr::rename_all(.funs=stringr::str_replace_all, pattern="\\.", replacement="_") %>%
    dplyr::mutate(tmp=paste(model, comparison, sex, sep=" ")) %>%
    dplyr::filter(tmp %in% models_to_keep) %>%
    dplyr::select(-tmp) %>%
    dplyr::mutate(study=stringr::str_replace(study, "MAYO", "MayoRNAseq"),
                  study=stringr::str_replace(study, "MSSM", "MSBB"),
                  sex=stringr::str_replace_all(sex,
                                               c("ALL"="males and females",
                                                 "FEMALE"="females only",
                                                 "MALE"="males only")),
                  model=stringr::str_replace(model, "\\.", " x "),
                  model=stringr::str_replace(model, "Diagnosis", "AD Diagnosis"),
                  logfc=round(logfc, digits = 3),
                  fc=2**logfc
    ) %>%
    dplyr::mutate(model=glue::glue("{model} ({sex})", model=model, sex=sex))
  # %>%
    # dplyr::select(-model) %>%
    # dplyr::rename(model=tmp)
}

#' Process rnaseq diff expr data
#'
#' @export
process_rnaseq_diff_expr_data <- function(data, gene_info, adj_p_value_threshold, target_list) {
  keep <- data %>%
    dplyr::filter(adj_p_val <= adj_p_value_threshold | (ensembl_gene_id %in% target_list$ensembl_gene_id),
                  ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    dplyr::select(ensembl_gene_id) %>%
    dplyr::distinct()

  data %>%
    dplyr::filter(ensembl_gene_id %in% keep$ensembl_gene_id) %>%
    assertr::chain_start() %>%
    assertr::verify(ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    assertr::chain_end() %>%
    dplyr::select(ensembl_gene_id, logfc, fc, ci_l, ci_r,
                  adj_p_val, tissue, study, model) %>%
    dplyr::left_join(.,
                     gene_info %>% dplyr::select(ensembl_gene_id, hgnc_symbol),
                     by="ensembl_gene_id") %>%
    filter(!is.na(hgnc_symbol)) %>%
    dplyr::select(ensembl_gene_id, hgnc_symbol, everything())
}

#' Get network data from Synapse
#'
#' @export
get_network_data <- function(id) {
  synapser::synGet(id)$path %>%
    readr::read_csv()
}

#' Process rnaseq diff expr data
#'
#' @export
process_network_data <- function(data, gene_info) {
  data %>%
    filter(geneA_ensembl_gene_id %in% gene_info$ensembl_gene_id,
           geneB_ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    dplyr::select(geneA_ensembl_gene_id,
                  geneB_ensembl_gene_id, brainRegion) %>%
    left_join(.,
              gene_info %>%
                dplyr::select(ensembl_gene_id,
                              hgnc_symbol) %>%
                distinct(),
              by=c("geneA_ensembl_gene_id"="ensembl_gene_id")) %>%
    dplyr::rename(geneA_external_gene_name=hgnc_symbol) %>%
    left_join(.,
              gene_info %>%
                dplyr::select(ensembl_gene_id,
                              hgnc_symbol) %>%
                distinct(),
              by=c("geneB_ensembl_gene_id"="ensembl_gene_id")) %>%
    dplyr::rename(geneB_external_gene_name=hgnc_symbol) %>%
    assertr::chain_start() %>%
    assertr::verify(geneA_ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    assertr::verify(geneB_ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    assertr::verify(geneA_external_gene_name %in% gene_info$hgnc_symbol) %>%
    assertr::verify(geneB_external_gene_name %in% gene_info$hgnc_symbol) %>%
    assertr::chain_end()
}

#' Get proteomics data
#'
#' @export
get_proteomics_data <- function(id) {
  synGet(id)$path %>%
    readr::read_csv() %>%
    dplyr::rename(ensembl_gene_id=ENSG) %>%
    dplyr::rename(hgnc_symbol=GeneName) %>%
    dplyr::rename_all(tolower)
}

#' Get metabolomics data
#'
#' @export
get_metabolomics_data <- function(id, env) {
  synGet(id)$path %>% load(envir = env)
}

#' Get entire Synapse table by id (works with both large and small tables)
#'
#' @export
get_whole_table <- function (id) {
    synapser::synTableQuery(glue::glue("select * from {id}"), resultsAs='csv')$filepath %>%
    readr::read_csv() %>%
    dplyr::select(-ROW_ID, -ROW_VERSION)
}

#' Load srm data from Synapse
#'
#' @export
get_srm_data <- function(id){
  synGet(id)$path %>% readr::read_csv()
}

#' Load neuropath data from Synapse
#'
#' @export
get_neuropath_data <- function(id){
  synGet(id)$path %>% readr::read_csv() %>% dplyr::rename_all(tolower)
}

#' Load target expression validation harmonized data from Synapse
#'
#' @export
get_target_exp_validation_harmonized_data <- function(id){
  synGet(id)$path %>% readr::read_csv()
}

#' Process metabolomics. For stats,
#' - convert list into dataframe
#' - make id column into factor
#' - transpose boxplot stats
#' - trim extraneous columns
#' Then, join processed stats and gene associations.
#'
#' @export
process_metabolomics <- function(raw_metabolomics_gene_associations,
                                 raw_metabolomics_stats) {
  standardize_columns <- function(df) {
    df %>%
      dplyr::rename_all(tolower) %>%
      dplyr::rename_all(.funs=stringr::str_replace_all, pattern="\\.", replacement="_")
  }

  check_shape <- function(x) { identical(dim(x),c(2,5)) }

  has_expected_shape <- function(data) {
    data %>%
      dplyr::select(transposed.boxplot.stats) %>%
      map(dim) %>%
      lapply(check_shape) %>%
      purrr::keep(.,FALSE) %>%
      length %>%
      all.equal(., 0)
  }

  stats_with_transposed <- raw_metabolomics_stats %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    assertr::verify(assertr::is_uniq(metabolite.id)) %>%
    dplyr::mutate(metabolite.id=factor(unlist(metabolite.id))) %>%
    dplyr::mutate(boxplot.stats=map(boxplot.stats,unlist)) %>%
    dplyr::mutate(transposed.boxplot.stats=lapply(boxplot.stats, t))

  assertthat::assert_that(has_expected_shape(stats_with_transposed))

  metabolomics_stats <- stats_with_transposed %>%
    dplyr::select(-metabolite.full.name,-boxplot.stats) %>%
    standardize_columns

  raw_metabolomics_gene_associations %>%
    standardize_columns %>%
    dplyr::left_join(metabolomics_stats, by="metabolite_id") %>%
    assertr::verify(nrow(.) == nrow(raw_metabolomics_gene_associations))
}

#' @export
process_data <- function(config) {
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

  geneInfoFinal <- process_rna_change_in_the_brain(geneInfoFinal, diffExprData)
  
  network <- agoradataprocessing::get_network_data(config$networkDataId)
  network <- agoradataprocessing::process_network_data(network, geneInfoFinal)

  proteomics <- get_proteomics_data(config$proteomicsDataId)

  omics_scores <- get_whole_table(config$omicsScoresTableId)
  genetics_scores <- get_whole_table(config$geneticsScoresTableId)
  overall_scores <- get_whole_table(config$overallScoresTableId) %>% dplyr::select(1:6)

  srm_data <- get_srm_data(config$srmDataId)
  target_exp_validation_harmonized_data <- get_target_exp_validation_harmonized_data(config$targetExpressionValidationHarmonizedId)

  geneInfoFinal <- process_proteomics_data(geneInfoFinal, proteomics)

  neuropath_data <- get_neuropath_data(config$neuropathDataId)

  env <- environment()
  get_metabolomics_data(config$metabolomicsDataId, env)
  metabolomics <- process_metabolomics(agora.metabolite.gene.associations,
                                       agora.metabolite.stats)

  # Data tests
  stopifnot(config$geneInfoColumns %in% colnames(geneInfoFinal))
  stopifnot(config$diffExprCols %in% colnames(diffExprData))
  stopifnot(config$networkCols %in% colnames(network))

   return(list(teamInfo=teamInfo, geneInfo=geneInfoFinal,
              diffExprData=diffExprData, network=network,
              proteomics=proteomics, metabolomics=metabolomics,
              omics_scores=omics_scores, 
              genetics_scores=genetics_scores,
              overall_scores=overall_scores,
              srm_data=srm_data,
              neuropath_data=neuropath_data,
              target_exp_validation_harmonized_data=target_exp_validation_harmonized_data
              ))
}

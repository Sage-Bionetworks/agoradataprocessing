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
    select(-`HGNC Name`) %>%
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
    assertr::chain_start() %>%
    assertr::verify(ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    assertr::chain_end() %>%
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
    dplyr::filter(adj_p_val <= adj_p_value_threshold | (ensembl_gene_id %in% target_list$ensembl_gene_id)) %>%
    dplyr::select(ensembl_gene_id) %>%
    distinct()

  data %>%
    dplyr::filter(ensembl_gene_id %in% keep$ensembl_gene_id) %>%
    assertr::chain_start() %>%
    assertr::verify(ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
    assertr::chain_end() %>%
    dplyr::select(ensembl_gene_id, logfc, fc, ci_l, ci_r,
                  adj_p_val, tissue, study, model) %>%
    left_join(.,
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
# Process gene expression (logfc and CI) data
# Drop existing gene symbol and add them later.

# filtering by p-value for a gene that is significant in at least one model or a
# nominated target
# colnames(diffExprData) <- gsub("\\.", "_", tolower(colnames(diffExprData)))
#
# diffExprDataFinal <-

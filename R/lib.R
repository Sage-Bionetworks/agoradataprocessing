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
process_igap <- function(igap, gene_info) {
  gene_info %>%
    dplyr::mutate(isIGAP=ensembl_gene_id %in% igap$ensembl_gene_id)
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
  gene_info %>%
    dplyr::left_join(data) %>%
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
process_target_list <- function(data) {
  data %>%
    dplyr::select(-hgnc_symbol) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::mutate(data_synapseid=stringr::str_split(data_synapseid, ",")) %>%
    dplyr::group_by(ensembl_gene_id) %>%
    tidyr::nest(.key='nominatedtarget') %>%
    dplyr::mutate(nominations=purrr::map_int(nominatedtarget, nrow))
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
process_brain_expression_data <- function(data, gene_info, fdr.random) {
  data %>%
    dplyr::mutate(isChangedInADBrain = fdr.random <= isChangedInADBrainThreshold) %>%
    dplyr::select(-fdr.random) %>%
    left_join(gene_info)
}

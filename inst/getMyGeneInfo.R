library(tidyverse)
library(mygene)

## Do things relying on mygene before loading synapser
# Get set of genes and their Ensembl Gene IDs as the universe of genes
refFlatFileUrl <- "https://gist.github.com/kdaily/2ed85e0dd3048fea8424b40243ddfa1c/raw/420086bd941962df66992667972c13462e504cc6/gencode.v24.primary_assembly.refFlat.txt"

geneTable <- readr::read_tsv(refFlatFileUrl, col_names=FALSE) %>%
  dplyr::select(ensembl_gene_id=X1) %>%
  dplyr::mutate(ensembl_gene_id=stringr::str_remove(ensembl_gene_id, "\\..*")) %>%
  dplyr::distinct()

# Use mygene.info to get name and summary
geneInfoRes <- mygene::getGenes(geneTable$ensembl_gene_id,
                                fields=c("symbol", "name", "summary",
                                         "type_of_gene", "alias", "go.MF"))

geneInfo <- geneInfoRes %>%
  as.data.frame() %>%
  as.tibble() %>%
  group_by(query) %>%
  top_n(1, X_score) %>%
  ungroup() %>%
  mutate(go.MF=purrr::map(go.MF, tibble::as_tibble))

geneTableMerged <- left_join(geneTable, geneInfo,
                             by=c("ensembl_gene_id"="query"))

save(geneTable, geneInfoRes, geneInfo, geneTableMerged, file="/tmp/geneInfo.RData")

synapser::synLogin()
f <- synapser::File("/tmp/geneInfo.RData", parentId="syn7525089")
f <- synapser::synStore(f, used=refFlatFileUrl)

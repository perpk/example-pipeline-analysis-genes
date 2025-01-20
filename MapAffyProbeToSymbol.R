require(biomaRt)

MapAffyProbeToSymbol <- function(probeids) {
  mart<-useMart("ENSEMBL_MART_ENSEMBL")
  mart<-useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart = mart,
    attributes = c("affy_hg_u133_plus_2", "external_gene_name"),
    filter = "affy_hg_u133_plus_2",
    values = probeids,
    uniqueRows = TRUE
  )
  annotLookup <- annotLookup[annotLookup$external_gene_name != "", ]
  unique_mapping <- annotLookup[!duplicated(annotLookup$affy_hg_u133_plus_2), ]
  probe_to_gene <- setNames(unique_mapping$external_gene_name, unique_mapping$affy_hg_u133_plus_2)
  gene_symbols <- sapply(probeids, function(probe) {
    if (probe %in% names(probe_to_gene)) {
      probe_to_gene[[probe]]
    } else {
      probe
    }
  })
  return(data.frame(geneSymbol = gene_symbols))
}

require(biomaRt)

GetAffyAnnotationLookupForProbe<-function(probeids) {
  mart<-useMart("ENSEMBL_MART_ENSEMBL")
  mart<-useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart = mart,
    attributes = c("affy_hg_u133_plus_2", "external_gene_name"),
    filter = "affy_hg_u133_plus_2",
    values = probeids,
    uniqueRows = TRUE
  )
  return(annotLookup);  
}

GetAffyAnnotationLookupForProbe(c("243934_at"))

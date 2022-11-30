library(STRINGdb)
library(igraph)
library(linkcomm)
library(data.table)

# wd <- readline("Insert the code directory: ")
# setwd(wd)

string_db <- STRINGdb$new(version="11", species=9606, score_threshold=700)

string.network <- string_db$get_graph()

genes <- read.csv("genes_for_HP_0030880.csv", sep=";")

genes_entrez <- string_db$map(my_data_frame = genes, my_data_frame_id_col_names = "GENE_SYMBOL")

hits.network <- string_db$get_subnetwork(genes_entrez$STRING_id)

# Obtain neighbours

first.neigh <- (neighbors(graph = string.network, v = V(hits.network)$name, mode = "all"))$name

hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, first.neigh)))

hits.df <- igraph::as_data_frame(hits.network, what="edges") 

lc <- linkcomm::getLinkCommunities(hits.df, hcmethod = "average")

# Saving plots

png(file="../results/Raynaud_genes-graph.png", width=850, height=850)
plot(
  hits.network,
  vertex.label = genes_entrez$GENE_SYMBOL,
  vertex.frame.color = NA,
  vertex.size = 10,
  edge.weight = edge.betweenness(hits.network),
  edge.color = "grey",
  reescale=TRUE,
  vertex.label.dist = 1.8,
  vertex.label.cex = 1,
  vertex.label.color = "black",
  layout = layout.auto(hits.network)
)
dev.off()

png(file="../results/Raynaud_genes-dendrogram.png", width=500, height=500)
plot(lc, type = "summary")
dev.off()


png(file="../results/Raynaud_genes-comunity_members_matrix.png", width=500, height=500)
plot(lc, type = "members")
dev.off()

png(file="../results/Raynaud_genes-comunities_graph.png", width=1080, height=1080)
plot(lc, type = "graph", layout = layout.fruchterman.reingold, ewidth = 2, vlabel = FALSE)
dev.off()


################# Nested communities ##################

# Nested communities: Son aquellas comunidades que son independientes respecto a otras
getAllNestedComm(lc)
png("../results/Raynaud-nested_communities.png")
plot(lc, type = "graph", vlabel=FALSE)
dev.off()

################# Estudio funcionalidad ##################

enriquecimiento <- function(cluster) {
  # Extraemos sus datos
  comunity_genes <- getNodesIn(lc, clusterids = cluster, type = "names")
  expr_c <- genes_entrez[genes_entrez$STRING_id %in% comunity_genes, ]
  
  # Enriquecimiento funcional con GO
  enrichmentGO <- string_db$get_enrichment(expr_c$STRING_id, category = "Process", iea = TRUE, )
  enrichmentGO$ontology <- rep("GO")
  
  # Enriquecimiento funcional con KEGG
  enrichmentKEGG <- string_db$get_enrichment(expr_c$STRING_id, category = "KEGG", iea = TRUE)
  enrichmentKEGG$ontology <- rep("KEGG")
  
  # Enriquecimiento funcional con Pfam
  # enrichmentPfam <- string_db$get_enrichment(expr_c$STRING_id, category = "Pfam", iea = TRUE)
  # enrichmentPfam$ontology <- rep("Pfam")
  
  # Juntamos las dos
  enrichment <- rbind(enrichmentGO)
  data.table::setcolorder(enrichment, c(11, c(1:10)))
  
  return(enrichment)
}

for (i in c(1:6)) {
  tryCatch(
    {
      comunidad <- enriquecimiento(i)
      dir <- paste("../results/Raynaud_Enriquecimiento_funcional_Comunidades/Comunity", i, ".csv", sep = "")
      write.csv(comunidad, dir)
    },
    error=function(err) {
      message(paste("Error in comunity", i, ":", sep = " "))
      message(err)
      message("\n")
      return(NA)
    }
  )
}

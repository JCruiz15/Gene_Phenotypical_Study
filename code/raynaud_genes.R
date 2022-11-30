library(STRINGdb)
library(igraph)
library(linkcomm)

wd <- readline("Insert the code directory: ")
setwd(wd)

string_db <- STRINGdb$new(version="11", species=9606, score_threshold=700)

string.network <- string_db$get_graph()

genes <- read.csv("genes_for_HP_0030880.csv", sep=";")

genes_entrez <- string_db$map(my_data_frame = genes, my_data_frame_id_col_names = "GENE_SYMBOL")


# g <- igraph::graph_from_data_frame(vertices = nodes$V2, d = links, directed = FALSE)

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
  vertex.size = 10,
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
# png("14_lc_nested_communities.png")
plot(lc, type = "graph")
# dev.off()

################# Estudio funcionalidad ##################

enriquecimiento <- function(cluster) {
  # Extraemos sus datos
  comunity_genes <- getNodesIn(lc, clusterids = i, type = "names")
  expr_c <- genes_entrez[genes_entrez$STRING_id %in% comunity_genes, ]
  
  # Enriquecimiento funcional con GO
  enrichmentGO <- string_db$get_enrichment(expr_c$STRING_id, category = "Process", iea = TRUE, )
  enrichmentGO$ontology <- rep("GO")
  
  # Enriquecimiento funcional con KEGG
  enrichmentKEGG <- string_db$get_enrichment(expr_c$STRING_id, category = "KEGG", iea = TRUE)
  enrichmentKEGG$ontology <- rep("KEGG")
  
  # Enriquecimiento funcional con Pfam
  #enrichmentPfam <- string_db$get_enrichment(expr_c$STRING_id, category = "Pfam", iea = TRUE)
  #enrichmentPfam$ontology <- rep("Pfam")
  
  # Juntamos las dos
  enrichment <- rbind(enrichmentGO, enrichmentKEGG)
  data.table::setcolorder(enrichment, c(11, c(1:10)))

  # Agrupamos filas repetidas indicando las ontolog?as de origen
  for (i in enrichment$term[!duplicated(enrichment$term)]) {
    num_filas <- which(enrichment$term == i)
    ontologias <- paste(enrichment[num_filas, 1], collapse = ", ")
    nueva_fila <- enrichment[num_filas[1],]
    nueva_fila$ontology <- ontologias
    enrichment <- enrichment[-num_filas, ]
    enrichment <- rbind(enrichment, nueva_fila)
  }
  enrichment
  return(enrichment)
}

for (i in c(1,2,3,4,5,6)) {
  comunidad <- enriquecimiento(i)
  View(comunidad)
  dir <- paste("../results/Raynaud_Enriquecimiento_funcional_Comunidades", i, ".csv", sep = "")
  write.csv(comunidad, dir)
}

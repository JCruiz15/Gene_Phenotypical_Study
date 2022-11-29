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

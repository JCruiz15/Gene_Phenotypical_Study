library(STRINGdb)
library(igraph)
library(linkcomm)

string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="")

nodes <- read.table("/home/jesusaldanamartin/Documentos/SistBio/raynaud/Gene_Phenotipycal_Study/code/string_node_degrees.tsv", header = FALSE)

df <- read.table("/home/jesusaldanamartin/Documentos/SistBio/raynaud/Gene_Phenotipycal_Study/code/string_interactions.tsv")

links <- data.frame(df$V3, df$V4)

# links
# nodes

g <- igraph::graph_from_data_frame(vertices = nodes$V2, d = links, directed = FALSE)

g

hits.network <- string_db$get_subnetwork(nodes$V2)

first.neigh <- (neighbors(graph = g, v = V(hits.network)$name, mode = "all"))$name

hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, first.neigh)))

hits.df <- igraph::as_data_frame(hits.network, what="edges") 

lc <- getLinkCommunities(hits.df, hcmethod = "average")

plot(lc, type = "members")

plot(
  hits.network,
  vertex.label = nodes$V1,
  vertex.size = 10,
  vertex.label.color = "black",
  layout = layout.auto(g)
)

plot(
  g,
  vertex.label = nodes$V1,
  vertex.size = 10,
  vertex.label.color = "black",
  layout = layout.auto(g)
     )


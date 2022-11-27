library(STRINGdb)
library(igraph)
library(linkcomm)

nodes <- read.table("string_node_degrees.tsv", header = FALSE)

df <- read.table("string_interactions.tsv")

links <- data.frame(df$V3, df$V4)

# links
# nodes

g <- igraph::graph_from_data_frame(vertices = nodes$V2, d = links, directed = FALSE)

plot(
  g,
  vertex.label = nodes$V1,
  vertex.size = 10,
  vertex.label.color = "black",
  layout = layout.auto(g)
     )

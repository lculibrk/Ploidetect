compressdata <- function(compress, x, segmentation_threshold = segmentation_threshold){
  time1 <- Sys.time()
  # Arrange data by position and remove tibble-ness
  compress <- compress %>% arrange(pos) %>% as.data.frame()
  # If this is the first iteration, set npoints to one (npoints records how many original vertices were compressed into one vertex)
  if(is.null(compress$npoints)){
    compress$npoints <- 1
  }
  # Get differences between neighbouring points
  diffs <- abs(diff(compress$corrected_depth/compress$npoints))
  # Initialize graph from data
  graph <- graph(edges = c(row.names(compress)[1], rep(row.names(compress[-c(1, nrow(compress)),]), each = 2), row.names(compress)[nrow(compress)]), directed = F)
  # Set edges to have the diffs attribute
  graph <- set_edge_attr(graph, name = "diff", value = diffs)
  # Give vertices appropriate attributes
  graph <- set_vertex_attr(graph, name = "corrected_depth", value = compress$corrected_depth)
  graph <- set_vertex_attr(graph, name = "npoints", value = compress$npoints)
  graph <- set_vertex_attr(graph, name = "from", value = compress$pos)
  graph <- set_vertex_attr(graph, name = "to", value = compress$end)
  # loop over all vertices
  sort_by_diff <- data.frame("vertex" = V(graph)$name)
  sort_by_diff$diff <- 0
  sort_by_diff$e
  time2 <- Sys.time()
  #for(row in 1:nrow(sort_by_diff)){
  #  sort_by_diff$diff[row] <- max(edge_attr(graph, "diff", incident(graph, sort_by_diff$vertex[row])))
  #  sort_by_diff$e[row] <- which.max(edge_attr(graph, "diff", incident(graph, sort_by_diff$vertex[row])))
  #}
  edges <- edge_attr(graph, "diff", E(graph))
  e1 <- c(NA, edges)
  e2 <- c(edges, NA)
  
  criteria <- segmentation_threshold*x
  ## Diff check
  both <- which(e1 > criteria & e2 > criteria)
  graph <- delete_edges(graph, c(both, na.omit(unique((1:length(e1)) - as.numeric(e1 > e2)))))
  #na.omit(unique((1:length(e1)) - as.numeric(e1 > e2)))
  
  
  time3 <- Sys.time()
  #sort_by_diff <- sort_by_diff %>% arrange(diff)
  #time26 <- Sys.time()
  #delete_edges(graph, E(graph)[edge_attr(graph, "diff") > segmentation_threshold*x])
  #todel <- c()
  #for(vertex in sort_by_diff$vertex){
  #  if(length(incident(graph, vertex)) == 0){
  #    next
  #  }
  #  # If vertex is an outlier (diffs are over some threshold fraction of what we expect for a copy change) then break all edges
  #  if(all(edge_attr(graph, "diff", incident(graph, vertex)) > segmentation_threshold*x)){
  #    todel <- c(todel, incident(graph, vertex))
  #    #graph <- delete_edges(graph, incident(graph, vertex))
  #    next
  #  }
  #  # If vertex has two edges, break the one with larger "diff"
  #  #if(length(incident(graph, vertex)) == 2){
  #  #  graph <- delete_edges(graph, incident(graph, vertex)[which.max(get.edge.attribute(graph, "diff", incident(graph, vertex)))])
  #  #}
  #  if(length(incident(graph, vertex)) == 2){
  #    todel <- c(todel, incident(graph, vertex)[which.max(get.edge.attribute(graph, "diff", incident(graph, vertex)))])
  #  }
  #}
  #time27 <- Sys.time()
  #todel <- E(graph)[edge_attr(graph, "diff") > segmentation_threshold*x]
  #graph <- delete_edges(graph, todel)
  #delete_edges(graph, todel)
  #time3 <- Sys.time()
  #time3-time26
  #time27-time26
  #print(get.vertex.attribute(toy, "npoints", V(toy)))
  #print(toy)
  # Get list of all vertex pairs to merge
  tomerge <- ends(graph, E(graph))
  # Get all vertices
  vertices <- V(graph)
  # Coerce vertices into a format where value is the vertex value and name is vertex name
  vertnames <- names(vertices)
  vertices <- as.numeric(vertices)
  names(vertices) <- vertnames
  # Change "tomerge" from names to values
  tomerge[,2] <- vertices[which(names(vertices) %in% tomerge[,2])]
  tomerge[,1] <- vertices[which(names(vertices) %in% tomerge[,1])]
  # Not needed I think
  #todelete <- vertices[which(vertices %in% tomerge[,2])]
  # Change pairs of vertices to repeat the same vertex twice (used in contract.vertices() to map which vertices to contract)
  vertices[which(vertices %in% tomerge[,2])] <- tomerge[,1]
  mode(vertices) <- "integer"
  lint <- vertices[1]
  time4 <- Sys.time()
  for(i in 2:length(vertices)){
    if(vertices[i] - 1 > lint){
      vertices[vertices == vertices[i]] <- lint + 1
    }
    lint <- vertices[i]
  }
  time5 <- Sys.time()
  
  # Merge connected vertices
  graph <- contract.vertices(graph, mapping = vertices, vertex.attr.comb = list("corrected_depth" = "sum", "npoints" = "sum", "from" = "min", "to" = "max", "name" = "first"))
  # Delete all old vertices
  #toy <- delete.vertices(toy, which(names(V(toy)) == "character(0)"))
  # Reconstruct the data.frame we began with
  dat <- data.frame("corrected_depth" = get.vertex.attribute(graph, "corrected_depth"), 
                    "npoints" = get.vertex.attribute(graph, "npoints"), 
                    "pos" = get.vertex.attribute(graph, "from"),
                    "end" = get.vertex.attribute(graph, "to"))
  #print(dat[which(dat$npoints == 1),])
  if(T){
    print(paste0("Iteration complete, segment count = ", nrow(dat)))
  }
  time6 <- Sys.time()
  
  #print(c(time2 - time1, time2-time3, time4-time3, time5-time4, time6-time5))
  #print(dat)
  return(dat)
}

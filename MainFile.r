library("GNAR")
library("igraph")

# Generate data from a random graph and GNAR
graph <- make_graph(edges = c(1, 2, 2, 3, 2, 4, 4, 5), n = 5, directed = FALSE) 
matrix <- as_adjacency_matrix(graph, type = "upper")
data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5,0.3)))

errors <- c(0)
for (i in 0:(2^10-1)){
    # Convert i to binary
    binary <- as.integer(intToBits(i))
    binary <- binary[1:10]
    # Convert binary to a matrix
    matrix <- matrix(0,5,5)
    matrix[upper.tri(matrix)] <- binary
    # Convert matrix to a graph
    net <- matrixtoGNAR(matrix)
    # Fit GNAR
    fit <- GNARfit(data[], net, alphaOrder = )
    # Calculate error
    errors <- c(errors, fit$error)
}
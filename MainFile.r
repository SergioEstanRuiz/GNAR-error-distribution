library("GNAR")
library("igraph")

# Generate data from a random graph and GNAR
graph <- make_graph(edges = c(1, 2, 2, 3, 2, 4, 4, 5), n = 5, directed = FALSE) 
realMatrix <- as_adjacency_matrix(graph, type = "upper")
data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
data.ts <- ts(data)

generateError <- function(size = 5){
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
    data.ts <- ts(data)
    errors <- c(0)
    numberofGraphs <- 2^(size*(size-1)/2)-1
    for (i in 1:numberofGraphs){
        # Convert i to binary
        binary <- as.integer(intToBits(i))
        binary <- binary[1:10]
        # Convert binary to a matrix
        matrix <- matrix(0,size,size)
        matrix[upper.tri(matrix)] <- binary
        # Convert matrix to a graph
        net <- matrixtoGNAR(matrix)
        # Fit GNAR
        fit <- GNARfit(vts = data.ts[1:199,], net = net, alphaOrder = 1, betaOrder = c(1))
        # Calculate error
        err <- sum((data[200] - predict(fit))^2)
        errors <- c(errors, err)
    }
        return(errors)
}

errors <- generateError(size=5)

pdf(file="./testFile.pdf")
# plot histogram of errors 
hist(errors, breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors")
# add description
mtext("Histogram representing distribution of errors for data \n simulated from a GNAR(1,[1]) for a 5-node network", side=3)
dev.off()


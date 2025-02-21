library("GNAR")
library("igraph")

int_to_mat <- function(i, size){
    # Convert i to binary
    binary <- as.integer(intToBits(i))
    length <- size*(size-1)/2
    binary <- binary[1:length]
    # Convert binary to a matrix
    matrix <- matrix(0,size,size)
    matrix[upper.tri(matrix)] <- binary
    return(matrix)
}

mat_to_int <- function(matrix){
    return(sum(matrix[upper.tri(matrix)]*2^(0:(length(matrix[upper.tri(matrix)])-1))))
}

generateError <- function(size = 5, k = 1, n=500, seed = NA, order = 1, globalAlpha = TRUE){
    if (!is.na(seed)){
        set.seed(seed)
    }
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    # Generate list of length size with random values between -0.5 and 0.5
    alphaParams <- lapply(1:order, function(x) runif(size, min=-0.5, max=0.5))
    betaParams <- lapply(1:order, function(x) runif(1, min=-0.5, max=0.5))
    data <- GNARsim(n=n, net = igraphtoGNAR(graph), alphaParams = alphaParams, betaParams = betaParams)

    data.ts <- ts(data)
    errors <- c()
    numberofGraphs <- 2^(size*(size-1)/2)-1
    for (i in 1:numberofGraphs){
        # Convert i to binary
        if (i %% 100 == 0) {
            print(i)
        }
        matrix <- int_to_mat(i, size)
        err <- 0
        net <- matrixtoGNAR(matrix)
        for (j in 1:k){
            fit <- GNARfit(vts = data.ts[j:(n-1-k+j),], net = net, alphaOrder = 2, betaOrder = c(1,1), globalalpha = globalAlpha)
            err <- err + sum((data[(n-k+j),] - predict(fit))^2)
        }
        errors <- c(errors, err/k)
    }
    return(errors)
}
getRandomSample <- function(size = 5, errors, sample = 100, erdos_prob=0.3){
    err <- c()
    for (i in 1:sample){
        graph <- erdos.renyi.game(size, p=erdos_prob, type = "gnp")
        matrix <- as_adjacency_matrix(graph, type = "upper")
        i <- mat_to_int(matrix)
        err <- c(err, errors[i])
    }
    return(err)
}
hamming_distance <- function(i, j){
    return(sum(as.integer(intToBits(i)) != as.integer(intToBits(j))))
}
edge_to_int <- function(edge,size){
    # Parameters
    # edge: list of vectors (a,b) where 1 < a < b, eg. list(c(4,5),c(2,3))
    # size: number of nodes
    
    index <- function(edge){
        a <- edge[1]
        b <- edge[2]
        return(a + (b-1)*(b-2)/2)
    }
    indices <- sapply(edge, index)
    result <- c()
    for (i in 1:(2^(size*(size-1)/2) - 1)) {
        if (all(as.integer(intToBits(i))[indices] == 1)) {
            result <- c(result, i)
        }
    }
    return(result)
}

error <- generateError(size=6, k=40, n=500, seed = 1000, order =2, globalAlpha = FALSE)
write.table(error, file = "./error007.txt", sep = "\t")

# Plot histogram of error and sample_errors
pdf(file="./hist007.pdf")
# plot histogram of errors
hist(as.numeric(error), breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors", xlim = range(c(error)), ylim = c(0, 10000))
legend("topright", legend = c("Errors"), fill = c("lightblue"))
dev.off()

error <- as.matrix(read.table("./error007.txt", sep = "\t"))
set.seed(1000)
true_graph <- erdos.renyi.game(6, p=0.3, type = "gnp")
true_int <- mat_to_int(as_adjacency_matrix(true_graph))
distances <- lapply(1:(2^15-1), function(x) hamming_distance(true_int, x))
pdf(file="./hHammingDistances007.pdf")
plot(error, distances, main = "Error vs Hamming Distance", xlab = "Error", ylab = "Hamming Distance")
dev.off()

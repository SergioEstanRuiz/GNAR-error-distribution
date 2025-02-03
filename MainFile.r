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

generateError <- function(size = 5, k = 1, seed = NA){
    if (!is.na(seed)){
        set.seed(seed)
    }
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,size))), betaParams = list(c(0.5)))
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
            fit <- GNARfit(vts = data.ts[j:(199-k+j),], net = net, alphaOrder = 1, betaOrder = c(1))
            err <- err + sum((data[(200-k+j),] - predict(fit))^2)
        }
        errors <- c(errors, err/k)
    }
    return(errors)
}

generateErrorAll <- function(size = 5, k = 10, seed = NA){
    if (!is.na(seed)){
        set.seed(seed)
    }
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    data <- GNARsim(n=250, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
    data.ts <- ts(data)
    numberofGraphs <- 2^(size*(size-1)/2)-1
    errors <- matrix(0, k, numberofGraphs)
    for (i in 1:numberofGraphs){
        # Convert i to binary
        binary <- as.integer(intToBits(i))
        binary <- binary[1:10]
        # Convert binary to a matrix
        matrix <- matrix(0,size,size)
        matrix[upper.tri(matrix)] <- binary
        # Convert matrix to a graph
        err <- 0
        net <- matrixtoGNAR(matrix)
        for (j in 1:k){
            fit <- GNARfit(vts = data.ts[j:(249-k+j),], net = net, alphaOrder = 1, betaOrder = c(1))
            err <- err + sum((data[(250-k+j),] - predict(fit))^2)
            errors[j,i] <- err
        }
    }
    for (j in 1:k){
        errors[j,] <- errors[j,]/j
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

error <- generateError(size=6, k=40, seed = 1234)
write.table(error, file = "./error_n=6_k=40.txt", sep = "\t")
error <- read.table("./error_n=6_k=40.txt", sep = "\t")
error <- as.matrix(error)
sample_error <- getRandomSample(size=6, errors = error, sample = 1000, erdos_prob = 0.1)

# Plot histogram of error and sample_errors
pdf(file="./hist(6,40).pdf")
# plot histogram of errors
hist(as.numeric(error), breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors", xlim = range(c(error, sample_error)), ylim = c(0, 4000))
hist(as.numeric(sample_error), breaks = 100, col = "#187534", border = "black", add = TRUE)
legend("topright", legend = c("Errors", "Sample Errors"), fill = c("lightblue", "#187534"))
dev.off()

summary(error)
summary(sample_error)
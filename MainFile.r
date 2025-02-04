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

generateError <- function(size = 5, k = 1, seed = NA, order = 1){
    if (!is.na(seed)){
        set.seed(seed)
    }
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    # Generate list of length size with random values between -0.5 and 0.5
    # 
    alphaParams <- list(c(runif(size*order, min=-0.5, max=0.5)))
    alphaParams <- lapply(1:order, function(x) runif(size, min=-0.5, max=0.5))
    betaParams <- lapply(1:order, function(x) runif(1, min=-0.5, max=0.5))
    data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = alphaParams, betaParams = betaParams)
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
            fit <- GNARfit(vts = data.ts[j:(199-k+j),], net = net, alphaOrder = 2, betaOrder = c(1,1))
            err <- err + sum((data[(200-k+j),] - predict(fit))^2)
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

error <- generateError(size=6, k=40, seed = 1234, order =2)
write.table(error, file = "./error_n=6_k=40_true=2_fit=2.txt", sep = "\t")
error <- read.table("./error_n=6_k=40_true=2_fit=2.txt", sep = "\t")
error <- as.matrix(error)
sample_error <- getRandomSample(size=6, errors = error, sample = 1000, erdos_prob = 0.1)

# Plot histogram of error and sample_errors
pdf(file="./hist(n=6,k=40,t=2,f=2).pdf")
# plot histogram of errors
hist(as.numeric(error), breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors", xlim = range(c(error, sample_error)), ylim = c(0, 10000))
hist(as.numeric(sample_error), breaks = 100, col = "#187534", border = "black", add = TRUE)
legend("topright", legend = c("Errors", "Sample Errors"), fill = c("lightblue", "#187534"))
dev.off()

# min_error <- 0
# for (i in 1:100){
#     sample_error <- getRandomSample(size=6, errors = error, sample = 5, erdos_prob = 0.1)
#     min_error <- min_error + min(sample_error)
# }
# print(min_error/100)
# summary(error)
# summary(sample_error)
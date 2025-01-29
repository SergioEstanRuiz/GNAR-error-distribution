library("GNAR")
library("igraph")

generateError <- function(size = 5, k = 1, seed = NA){
    if (!is.na(seed)){
        set.seed(seed)
        print("Seed is set")
    }
    graph <- erdos.renyi.game(size, p=0.3, type = "gnp")
    data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
    data.ts <- ts(data)
    errors <- c()
    numberofGraphs <- 2^(size*(size-1)/2)-1
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

# set.seed(1234)
# graph <- erdos.renyi.game(5, p=0.3, type = "gnp")
# data <- GNARsim(n=250, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
# data.ts <- ts(data)
# net <- igraphtoGNAR(graph)
# err <- 0
# for (j in 1:40){
#     fit <- GNARfit(vts = data.ts[j:(249-40+j),], net = net, alphaOrder = 1, betaOrder = c(1))
#     err <- err + sum((data[(250-40+j),] - predict(fit))^2)
# }
# err <- err/40
# summary(as.numeric(errorsAll[40,]))

error <- generateError(size=5, k=40, seed = 1234)
write.table(error, file = "./error_n=5_k=40.txt", sep = "\t")
error <- read.table("./error_n=5_k=40.txt", sep = "\t")



pdf(file="./histogram_n=5_k=40_part2.pdf")
# plot histogram of errors 
hist(as.numeric(error), breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors")
mtext("Histogram representing distribution of errors for data \n simulated from a GNAR(1,[1]) for a 5-node network, k=40", side=3)
dev.off()


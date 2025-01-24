library("GNAR")
library("igraph")

# Generate data from a random graph and GNAR
graph <- make_graph(edges = c(1, 2, 2, 3, 2, 4, 4, 5), n = 5, directed = FALSE) 
realMatrix <- as_adjacency_matrix(graph, type = "upper")
data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
data.ts <- ts(data)

generateError <- function(size = 5, k = 1, seed = NA){
    if (!is.na(seed)){
        set.seed(seed)
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
    data <- GNARsim(n=200, net = igraphtoGNAR(graph), alphaParams = list(c(rep(0.2,5))), betaParams = list(c(0.5)))
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
            fit <- GNARfit(vts = data.ts[j:(199-k+j),], net = net, alphaOrder = 1, betaOrder = c(1))
            err <- err + sum((data[(200-k+j),] - predict(fit))^2)
            errors[j,i] <- err
        }
    }
    for (j in 1:k){
        errors[j,] <- errors[j,]/j
    }
        return(errors)
}
plotting <- function(x, mean, std, xlab = "x-axis", ylab = "y-axis", title = "Title of Plot"){
    
    set.seed(1234)
    df <- data.frame(x =x, F = mean, L = mean - std, U = mean + std)

    upperlim <- max(mean) + 1.5*max(std)
    lowerlim <- min(mean) - 1.5*max(std)

    plot(df$x, df$F, ylim = c(lowerlim, upperlim), type = "l", xlab = xlab, ylab = ylab, main = title)
    #make polygon where coordinates start with lower limit and 
    # then upper limit in reverse order
    polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
    lines(df$x, df$F, lwd = 2)
    #add red lines on borders of polygon
    lines(df$x, df$U, col="red",lty=2)
    lines(df$x, df$L, col="red",lty=2)
}

errors1 <- generateError(size=5, k=2, seed= 1234)
errors2 <- generateError(size=5, k=2, seed = 1234) # good one
print(errors1-errors2)
print(errors2- errorsAll[2,])
errorsAll <- generateErrorAll(size=5, k=10, seed = 1234)
x <- 1:10
mean <- apply(errorsAll, 1, mean)
std <- apply(errorsAll, 1, sd)

pdf(file="./histogram_n=5_k=5to10.pdf")
colors <- rainbow(5)
plot(NULL, xlim = c(4, 5.2), ylim = c(0, 10), xlab = "Error", ylab = "Density", main = "Density Plots of Errors for Each k")
for (i in 1:5) {
    lines(density(errorsAll[(i+5),]), col = colors[i], lwd = 2)
}
legend("topright", legend = paste("k =", 6:10), col = colors, lwd = 2)
dev.off()

pdf(file="./histogram_n=5_k=2.pdf")
# plot histogram of errors 
hist(errors, breaks = 100, col = "lightblue", border = "black", xlab = "Error", main = "Histogram of Errors")
# add description
mtext("Histogram representing distribution of errors for data \n simulated from a GNAR(1,[1]) for a 5-node network, k=1", side=3)
dev.off()


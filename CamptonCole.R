
# The purpose of this project is to record the relationship of several factors, the average number of
# paths between 2 nodes (AvPaths), the % of time a graph has a Euler Path, The average largest eigenvalue
# and the average diameter of a graph with respect to the number of nodes and number of paths. 
# The result is a matrix M that has columns with the 4 properties for each number of nodes, and with
# rows representing the number of paths.

#update.packages()
library(igraph)
library(igraphdata)
library(gtools)

#Control Parameters
NODES = 10
EDGES = 10
LOOPS = 100
set.seed(37251)
## Plot for euler Paths alone - because of forced axis scaling to 10x10.
#plot(0:10, 0:10,xlab="Nodes",ylab="% Euler Paths", yaxt='n')
plot(0:10, 0:10,xlab="Nodes",ylab="")
par(xpd=TRUE)
colors = heat.colors(NODES)

## Create random Erdos-Renyi graph of n nodes and e edges
## Range of nodes
n = 2:NODES
## Range of of edges for each number of nodes
e = 1:EDGES
# Creates a matrix for data which has columns of node number rows of edge numbers
M = matrix(nrow=length(e), ncol=4*length(n))
colnames = list()
# Creates labels for the matrix
for( i in n){
  nodenum = paste(as.character(i),'nodes AvPaths')
  cname = c(nodenum,'%Euler','AvEVCent','AvDiameter')
  colnames = c(colnames,cname)
}
colnames(M) = colnames

for(i in n){
  x = matrix(nrow=(i), ncol=2)
for(j in e){
  
# Because the graph is simple the graph can have at most the complete number of edges
# The complete graph of n nodes edges given by the triangular number n(n-1)/2
if(j <= (i*(i-1)/2)){

  
##Loop through loops times for average values for each combination of i nodes and j edges.
loops = 1:LOOPS
avPath = 0
ep = 0
totsubs = 0
centerDeg = 0
diameter = 0
for(k in loops){
  
## Create random Euler-Renyi graph of i(1-n) nodes and j(1-e) edges
G = erdos.renyi.game(i,j,type="gnm")


## Calculates number of paths between Nodes ni and nj and finds mean path length.
nodes = 1:i
# Finds the combinations (without replacement) of nodes
combs = combinations(length(nodes),2,nodes)

# Finds number of edges between 2 edges
f <- function(x) {
  e = edge.disjoint.paths(G, x[1], x[2])
}
# Finds average number of paths between two nodes by applying f over the rows (combinations of nodes)
avPath = avPath + mean(apply(combs,1,f))


## Calculates if an Euler path exists. If 0 or 2 nodes have odd degree one exists
deg = degree(G)
#odds = sum(odd(deg))
subgraphs = clusters(G)
oddDegrees = array(0,dim=subgraphs$no)
for (n in 1:length(subgraphs$membership)){
  if (odd(deg[n])){
  oddDegrees[subgraphs$membership[n]] = oddDegrees[subgraphs$membership[n]] + 1
  }
}
for( odds in oddDegrees){
  if( odds == 0 || odds == 2){
   ep = ep + 1.0
  }
}
totsubs = totsubs + length(oddDegrees)

## Largest centrality measure from Eigenvector corresponding to largest eigenvalue
ev = evcent(G)$vector
maxi = which(ev == max(ev))
centerDeg = centerDeg + degree(G,maxi[1])

##Diameter
diameter = diameter + diameter(G)

}
## Matrix of N by E
M[j,4*(i-1)-3] = signif(avPath/length(loops),digits=4)
M[j,4*(i-1)-2] = ep/totsubs*10
M[j,4*(i-1)-1] = centerDeg/length(loops)
M[j,4*(i-1)] = diameter/length(loops)

}
}
########Plots!

# Plot 1 - Average Diameter
## Y is the output that is being correlated, X is the number of nodes. 
y <- M[,4*(i-1)]
x <- 1:length(y)
# For more than 2 nodes (where a linear model is relevant), fits a model to all data for a given number of nodes
## Can uncomment for closer fitting for 10th degree polynomial 'lm' -> linear model
if(i>2){
  #lo <- loess(y~x)
  lm <- lm(y ~ poly(x, 10, raw=TRUE))
  p1 <- points(x,y, pch=20, col=colors[i])
  lines(predict(lm), col=colors[i], lwd=2)
  #xspline(x,y)
}
legend(-3,-3, legend = 2:NODES, col=colors, pch=20, horiz=TRUE, cex=.65, title = "Average Graph Diameter vs Number of Paths")

# # Plot 2 - Average Maximum Eigenvalue Centrality
# ## Y is the output that is being correlated, X is the number of nodes. 
# y <- M[,4*(i-1)-1]
# x <- 1:length(y)
# # For more than 2 nodes (where a linear model is relevant), fits a model to all data for a given number of nodes
# ## Can uncomment for closer fitting for 10th degree polynomial 'lm' -> linear model
# if(i>2){
#   lo <- loess(y~x)
#   #lm <- lm(y ~ poly(x, 10, raw=TRUE))
#   p2 <- points(x,y, pch=20, col=colors[i])
#   lines(predict(lo), col=colors[i], lwd=2)
#   #xspline(x,y)
# }
# legend(-3,-3, legend = 2:NODES, col=colors, pch=20, horiz=TRUE, cex=.65, title = "Average Maximum Eigenvalue vs Number of Paths")


#if(i>2){ ggplot2.multiplot(p1, p2, cols=2)}


# # Plot 3 - Pecentage of the time there are Euler Paths
# ## Y is the output that is being correlated, X is the number of nodes. 
# y <- M[,4*(i-1)-2]
# x <- 1:length(y)
# # For more than 2 nodes (where a linear model is relevant), fits a model to all data for a given number of nodes
# ## Can uncomment for closer fitting for 10th degree polynomial 'lm' -> linear model
# if(i>2){
#   lo <- loess(y~x)
#   #lm <- lm(y ~ poly(x, 10, raw=TRUE))
#   points(x,y, pch=20, col=colors[i])
#   lines(predict(lo), col=colors[i], lwd=2)
#   #xspline(x,y)
# }
# legend(-3,-3, legend = 2:NODES, col=colors, pch=20, horiz=TRUE, cex=.65, title = "Percentage Euler Path vs Number of Paths")
 

# # Plot 4 - Average Path Length
# ## Y is the output that is being correlated, X is the number of nodes. 
# y <- M[,4*(i-1)-3]
# x <- 1:length(y)
# # For more than 2 nodes (where a linear model is relevant), fits a model to all data for a given number of nodes
# ## Can uncomment for closer fitting for 10th degree polynomial 'lm' -> linear model
# if(i>2){
#   lo <- loess(y~x)
#   #lm <- lm(y ~ poly(x, 10, raw=TRUE))
#   points(x,y, pch=20, col=colors[i])
#   lines(predict(lo), col=colors[i], lwd=2)
#   #xspline(x,y)
# }
# legend(-3,-3, legend = 2:NODES, col=colors, pch=20, horiz=TRUE, cex=.65, title = "Average Path Length vs Number of Paths")

}
print(M)
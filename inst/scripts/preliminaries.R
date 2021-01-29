library(POP)
library(sp)
library(nabor)
library(Matrix)
 
data(POPhexagons)

#set data directory of package as working directory
setwd('/home/jay/data/2019_papers/POP/POP_package/POP/data')

# create a list IDs of nearest neighbors within 30000 meters
# maximum will be 6 nearest neighbors for hexagonal grids
xy = coordinates(POPhexagons)
knnout = knn(xy, k = 7)
nb_list = vector("list", dim(knnout$nn.idx)[1]) 
for(i in 1:dim(knnout$nn.idx)[1])
  nb_list[[i]] = knnout$nn.idx[i,2:7][knnout$nn.dist[i,2:7] < 30000]

# Let W be indicator matrix of neighbors
# Create W matrix as sparse matrix
sparsei = NULL
sparsej = NULL
for(i in 1:length(nb_list)) {
	sparsei = c(sparsei,rep(i, times = length(nb_list[[i]])))
	sparsej = c(sparsej, nb_list[[i]])
}
W = sparseMatrix(x = 1, i = sparsei, j = sparsej)

# Let D be diagonal matrix with number of neighbors
# create D matrix as sparse matrix
Dmat = sparseMatrix(i = 1:length(nb_list), j = 1:length(nb_list),
	x = unlist(lapply(nb_list, length)))

#pre-compute determinant of I-rho*W, and store
# both determinants and sparse matrices of I-rho*W
# for grid of rho values to use as look-up during MCMC
logit_rho_grid = (-40:40)/5
rhogrid = exp(logit_rho_grid)/(1 + exp(logit_rho_grid))
detsig = NULL
SigiList = NULL
for(rhoi in rhogrid) {
		SigiList = c(SigiList, Dmat - rhoi*W)
		detsig = c(detsig, determinant(Dmat - rhoi*W, 
			log = TRUE)$modulus)
}
    
# save these preliminaries as data objects that are available from the package
save(rhogrid, file = 'rhogrid.rda')
save(SigiList, file = 'SigiList.rda')
save(detsig, file = 'detsig.rda')
save(nb_list, file = 'nb_list.rda')
save(W, file = 'W.rda')
save(Dmat, file = 'Dmat.rda')


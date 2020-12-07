# hp_functions.R

fit_gamma <- function(com, margin){

	if (margin==1){
		m <- mean(rowSums(com))
		v <- var(rowSums(com))
	}

	if (margin==2){
		m <- mean(colSums(com))
		v <- var(colSums(com))
	}
	
	scale <- v/m
	shape <- m*m/v
    return(c("scale"=scale, "shape"=shape))
}


# the following were taken (and in some cases modified) from melmasri/HP-prediction

plot_Z<-function(Z, xlab, ylab, ...){
    ## ploting interaction matrix as a binary image
    if(missing(ylab)) ylab = 'hosts'
	if(missing(xlab)) xlab = 'parasites'
    # par(mar = c(5,5,1,1)+0.1)
	image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1,]),
		col = c('white', 'black'), ylab=ylab, xlab=xlab,
          useRaster=TRUE,srt=45, axes=TRUE,cex.lab=1)
}

lof<-function(Z, indices = FALSE){
    ## Given a binary matrix Z. Where the rows is fixed
    ## A function that left orders the matrix sequentially from row 1 to n
    ## based on first appearance of columns.
    if(min(range(Z))<0) stop('Range is less that 0.')

    active_col <- apply(Z,1,function(r) which(as.vector(r)>0))
	bank = active_col[[1]]
	for(i in 1:nrow(Z)){
			a = setdiff(active_col[[i]], bank)
			bank=c(bank,a)
		}
    if(indices) bank else  Z[,bank]
}

plot_degree <- function(Z, Z_est, type='both', host.col='blue', parasite.col='red'){
    ## Plots the degree distribution per marginal on a bipartite biadjacency matrix
    ## input
    ## Z - interaction matrix
    ## Z_est = posterior interaction matrix (optional)
    ## type  = hosts, parasites, both
    ## extra input for colours 

    ## Optional estimated presence/absence matrix (Z_est) can be added to existing plot.
    para_degrees <- as.data.frame(table(colSums(Z)))
    para_degrees$Var1 <- as.numeric(para_degrees$Var1)
    ## para_degrees = para_degrees[-which(para_degrees$Var1<2),]
    host_degrees <- as.data.frame(table(rowSums(Z)))
    host_degrees$Var1 <- as.numeric(host_degrees$Var1)
    ## host_degrees = host_degrees[-which(host_degrees$Var1<2),]

    xlim = c(1, max(para_degrees$Var1,host_degrees$Var1)*1.5)
    ylim = c(1, max(para_degrees$Freq,host_degrees$Freq)*1.5)

    if (!missing(Z_est)){
        para_est <- as.data.frame(table(colSums(Z_est)))
        para_est$Var1 <- as.numeric(para_est$Var1)
        ## para_est = para_est[-which(para_est$Var1<2),]
        host_est <- as.data.frame(table(rowSums(Z_est)))
        host_est$Var1 <- as.numeric(host_est$Var1)
        ## host_est = host_est[-which(host_est$Var1<2),]
        xlim = c(1, max(para_degrees$Var1,host_degrees$Var1,
            para_est$Var1, host_est$Var1)*1.5)
        ylim = c(1, max(para_degrees$Freq,host_degrees$Freq,
            para_est$Freq, host_est$Freq)*1.5)
    }
    gpch = c('+', '*')
    if(type=='parasites'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        legend(xlim[2]*0.03, ylim[2]*1.4, c("Parasites"), col = parasite.col,
               pch = gpch[2], bty="n",pt.cex=1.5, cex=2)
        if(!missing(Z_est)){
            points((para_est), type="p", col=parasite.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.8, c("Est"), col = parasite.col,
               pch = 16, bty='n', pt.cex=1.5, cex=2)
        }
    }
    if(type=='hosts'){
        plot((host_degrees), type="p", col=host.col, pch=gpch[1], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        legend(xlim[2]*0.03, ylim[2]*1.4, c("Hosts"), col = host.col,
               pch = gpch[1], bty="n", pt.cex=1.5, cex=2)
        if(!missing(Z_est)){
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.8, c("Est"), col = host.col,
                   pch = 16, bty='n',pt.cex=1.5, cex=2)
        }
    }
    if(type=='both'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = 1.5,cex.axis = 1.5)
        points((host_degrees), type="p", col=host.col, pch=gpch[1])
    legend(xlim[2]*0.03, ylim[2]*1.4, c("Parasites", "Hosts"), col = c(parasite.col, host.col),
           pch = gpch[2:1], bty = 'n', pt.cex=1.5, cex=2)
        if (!missing(Z_est)) {
            points((para_est), type="p", col=parasite.col, pch=16)
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.03, ylim[2]*0.33, c("Est"), col = c("black"),
                   pch = 16, bty='n', pt.cex=1.5, cex=2)
        }
    }
}


# Fast EB scaling

eb.phylo<-function(phy, heights, a){
    ## modified form .eb.phylo
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R
	## #ht=heights.phylo(phy)
	## N=Ntip(phy)
	## Tmax=ht$start[N+1]
	## mm=match(1:nrow(ht), phy$edge[,2])
	## ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	## ht$t2=ht$start-ht$end+ht$t1
    if(a==0) return(phy)
    bl = (exp(a*heights$t2)-exp(a*heights$t1))/(a)
    phy$edge.length=bl[phy$edge[,2]]
    phy
}

heights.phylo<-function(x){
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res
}


generate_interactions<-function(r, c, eta = 0.3, aj=0.5, ai=0.5, D, mcmc_X=10){

    ## Generate a sample of Z, D, w and y
    i = floor(r*1.5)
    j = floor(c*1.5)
    # print(sprintf('Attempting an %ix%i matrix, with (r,c)=(%i,%i) truncation', i,j,r,c))

    distance_matrix<-function(n){
        d = matrix(0, ncol=n, nrow=n)
        a = rexp(sum(upper.tri(d)),5)
        d[upper.tri(d)]<-a
        d = t(d)
        d[upper.tri(d)]<-a
        d
    }
    if(missing(D)){
        D = distance_matrix(i)
    }
    else{
        # print(sprintf('Passed similarity matrix of dim %ix%i', nrow(D), ncol(D)))        
        if(i >= nrow(D)){
            i = nrow(D)
            # warning(paste0('The number of hosts is lowered to: ', nrow(D)))
        }else{
            minset = sample(1:nrow(D), i)
            D = D[minset, ]
            D = D[, minset]
        }
    }

    De = D^eta
    # This needs to be updated to new EB matrix transfomation

    w = rgamma(j ,aj, 1)
    y = rgamma(i, ai, 1)
    ## Setting of first interaction
    yw = outer(y,w)
    
    Z=matrix(0, i,j)
    # print(sprintf('Running an MCMC for %i iterations',mcmc_X*i))
    ##Sampling interactions
    for(s in 1:(mcmc_X*i)){
        sam = sample(1:i, j, replace=TRUE)
        Diag = cbind(sam,1:j)
        dd= (De%*%Z)[Diag]
        dd = 1*(dd==0) + dd*1*(dd>0)
        Z[Diag]<- 1*(runif(j)<= 1-exp(-yw[Diag]*dd))
    }

    # preserving rownames before pruning
    rownames(Z) <- rownames(D)

    # print('Removing parasites with no interactions...')
    ## Removing empty or extra columns
    aux = which(colSums(Z)==0)
    Z = Z[,-aux]
    w = w[-aux]
    dim(Z)

    # print('Removing hosts with no interactions...')
    ## Removing empty or extra rows
    aux = which(rowSums(Z)==0)
    Z = Z[-aux,]
    y = y[-aux]
    dim(Z)

    # if(ncol(Z)>c){
    #     aux = order(colSums(Z), decreasing=FALSE)
    #     aux = aux[1:(ncol(Z)-c)]
    #     Z = Z[,-aux]
    #     w = w[-aux]
    # }else{
    #     if(ncol(Z)!=c)
    #         stop("No. columns less than c")
    # }
    
    # # ## Removing row numbers
    # aux = which(rowSums(Z)<1)
    # if(length(aux)>0){
    #     Z = Z[-aux,]
    #     y = y[-aux]
    #     D = D[-aux,]
    #     D = D[,-aux]
    # }

    ## Re ordering
    # mohamad originally had "indeces" instead of "indices" which produces an error
    # bank = lof(Z, indices = TRUE)
    # w = w[bank]
    # Z = Z[,bank]
    list(w=w, y=y, eta=eta, Z=Z, D=D)
}

remove_noninteracting <-function(Z){
    ## Removing empty or extra columns
    aux = which(colSums(Z)==0)
    Z = Z[,-aux]
   
    ## Removing empty or extra rows
    aux = which(rowSums(Z)==0)
    Z = Z[-aux,]
   
    return(Z)
}  


# Plot trees with certain tips highlighted
# modified from Liam Revell's (phytools) 'droptipTree' 
# function within fancyTree.R 
# https://github.com/liamrevell/phytools/blob/master/R/fancyTree.R

droptipTree<-function(tree,tip,tip.labs){
	edges<-rep(0,nrow(tree$edge))
	names(edges)<-tree$edge[,2]
	keep<-setdiff(tree$tip.label,tip)
	ca<-findMRCA(tree,keep)
	root.node<-length(tree$tip)+1
	if(ca!=root.node){
		z<-setdiff(getDescendants(tree,root.node),getDescendants(tree,ca))
		edges[as.character(z)]<-1
	}
	z<-getDescendants(tree,ca)
	foo<-function(x,tree){
		n<-length(tree$tip.label)
		y<-getDescendants(tree,x)
		y<-y[y<=n]
		return(y)
	}
	y<-lapply(z,foo,tree=tree)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%keep)) edges[as.character(z[i])]<-1
	# par(mfrow=c(2,1))
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE, show.tip.label=tip.labs, label.offset=7)
	# dtree<-drop.tip(tree,tip); dtree$root.edge<-max(nodeHeights(tree))-max(nodeHeights(dtree))
	# plot.phylo(dtree,edge.width=2,no.margin=TRUE,root.edge=TRUE)
	# return(dtree)
}


sim_HP <- function(n_host=150, n_para=500, phy, method="elmasri", 
                    h_gamma, p_gamma, D) {

        if (method=="binomial"){
            hp <- matrix(rbinom(n_host*n_para, 1, 0.007), nrow=n_host, ncol=n_para)
            hp <- remove_noninteracting(hp)
            return(hp)
        }

        if (method=="elmasri"){
            n_host=nrow(D)
            hp <- generate_interactions(r=n_host, c=n_para, 
            eta = 1, ai=h_gamma["shape"], aj=p_gamma["shape"], D=D, 
            mcmc_X=1)
            return(hp$Z)
        }
}


# plot degree distribution across simulations

dd_overlay_plot <- function(Z, para=TRUE, col=rgb(0,0,0,0.1), main=""){

    if (para==TRUE){
        prich_max <- max(unlist(lapply(Z, rowSums)))
        max_nhosts <- max(unlist(lapply(Z, function(x) dim(x)[1])))

        Z[[1]] %>% rowSums() %>% sort() %>% rev() %>% log() %>%
                plot(., pch=20, main=main, ylab="parasite richness", xlab="Node rank", col=rgb(0,0,0,0.1), xlim=c(0,max_nhosts), ylim=c(0,log(prich_max)))

        for(i in 2:length(Z)){
            Z[[i]] %>% rowSums() %>% sort() %>% rev() %>%
                            log() %>% points(., pch=20, col=col)            
        }
    }   

    if (para==FALSE){           
        hrich_max <- max(unlist(lapply(Z, colSums)))
        max_nparas <- max(unlist(lapply(Z, function(x) dim(x)[2])))

        Z[[1]] %>% colSums() %>% sort() %>% rev() %>% log() %>%
                plot(., pch=20, main=main, ylab="host richness", xlab="Node rank", col=rgb(0,0,0,0.1), xlim=c(0,max_nparas), ylim=c(0,log(hrich_max)))

        for(i in 2:length(Z)){
            Z[[i]] %>% colSums() %>% sort() %>% rev() %>%
                            log() %>% points(., pch=20, col=col)
        }

    }

}




# ses.maxD
# code "modified" by Max Farrell (stolen) from Steve Kembel
require(picante)

maxD <- function (samp, dis, abundance.weighted=FALSE) {

    N <- dim(samp)[1]
    maxD <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            maxD[i] <- max(sample.dis)        
        }
        else {
            maxD[i] <- 0
        }
    }
    maxD
}


ses.maxD <- function (samp, dis, null.model = c("taxa.labels", "richness", 
    "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
    "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000) 
{
    dis <- as.matrix(dis)
    maxD.obs <- maxD(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    maxD.rand <- switch(null.model, taxa.labels = t(replicate(runs, 
        maxD(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted))), 
        richness = t(replicate(runs, maxD(randomizeMatrix(samp, 
            null.model = "richness"), dis, abundance.weighted))), 
        frequency = t(replicate(runs, maxD(randomizeMatrix(samp, 
            null.model = "frequency"), dis, abundance.weighted))), 
        sample.pool = t(replicate(runs, maxD(randomizeMatrix(samp, 
            null.model = "richness"), dis, abundance.weighted))), 
        phylogeny.pool = t(replicate(runs, maxD(randomizeMatrix(samp, 
            null.model = "richness"), taxaShuffle(dis), abundance.weighted))), 
        independentswap = t(replicate(runs, maxD(randomizeMatrix(samp, 
            null.model = "independentswap", iterations), dis, 
            abundance.weighted))), trialswap = t(replicate(runs, 
            maxD(randomizeMatrix(samp, null.model = "trialswap", 
                iterations), dis, abundance.weighted))))
    maxD.rand.mean <- apply(X = maxD.rand, MARGIN = 2, FUN = mean, 
        na.rm = TRUE)
    maxD.rand.sd <- apply(X = maxD.rand, MARGIN = 2, FUN = sd, 
        na.rm = TRUE)
    maxD.obs.z <- (maxD.obs - maxD.rand.mean)/maxD.rand.sd
    maxD.obs.rank <- apply(X = rbind(maxD.obs, maxD.rand), MARGIN = 2, 
        FUN = rank)[1, ]
    maxD.obs.rank <- ifelse(is.na(maxD.rand.mean), NA, maxD.obs.rank)
    data.frame(ntaxa = specnumber(samp), maxD.obs, maxD.rand.mean, 
        maxD.rand.sd, maxD.obs.rank, maxD.obs.z, maxD.obs.p = maxD.obs.rank/(runs + 
            1), runs = runs, row.names = row.names(samp))
}

# data(phylocom)
# test <- ses.maxD(phylocom$sample, cophenetic(phylocom$phylo), null.model="phylogeny.pool")



creepjump <- function (samp, dis, abundance.weighted=FALSE) {

    N <- dim(samp)[1]
    maxD <- numeric(N)
    mntd <- numeric(N)
    creep_jump <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            diag(sample.dis) <- NA
            maxD[i] <- max(sample.dis)
            mntd[i] <- mean(apply(sample.dis,2,min,na.rm=TRUE))
            creep_jump[i] <- mntd[i]/maxD[i]        
        }
        else {
            maxD[i] <- 0
            mntd[i] <- 0        
            creep_jump[i] <- 0        
        }
    }
    creep_jump
}




calc_sim_variables <- function(original_mat, ext_mat, n_trees, n_sims, ext_type,int_type){

# setting up dataframe
df <- data.frame(tree=rep(1:n_trees, each=n_sims), sim=rep(1:n_trees, times=n_trees), 
    extinction=ext_type, interactions=int_type)

## parasite responses
npara <- sapply(original_mat, ncol)
npara_ext <- sapply(ext_mat, function(x) sum(colSums(x)>0))

npara_lost <- npara-npara_ext

nSHP <- sapply(original_mat, function(x) sum(colSums(x)==1))
nSHP_ext <- sapply(ext_mat, function(x) sum(colSums(x)==1))
nSHP_lost <- nSHP - nSHP_ext

# Add to data frame
df$npara <- npara
df$npara_ext <- npara_ext
df$npara_lost <- npara_lost


# mean host richness
HR_mean <- sapply(sapply(original_mat, function(x) colSums(x)), mean)
HR_ext_mean <- sapply(sapply(ext_mat, function(x) colSums(x[,colSums(x)>0])), mean)
HR_mean_change <- HR_mean-HR_ext_mean


# Add to data frame
df$HR_mean <- HR_mean
df$HR_ext_mean <- HR_ext_mean
df$HR_mean_change <- HR_mean_change



# Phylo metrics

# get extant trees
trees_extant <- NULL

for (phy in 1:n_trees){

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim

            trees_extant[[cntr]] <- drop.tip(trees[[phy]], setdiff(trees[[phy]]$tip.label, rownames(ext_mat[[cntr]])))

    }
    
}

MPD <- NULL
MPD_ext <- NULL

for (phy in 1:n_trees){

    coph <- cophenetic(trees[[phy]])

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim
            MPD[[cntr]] <- mpd(t(original_mat[[cntr]]),coph)
            MPD[[cntr]][is.na(MPD[[cntr]])] <- 0

            coph_ext <- cophenetic(trees_extant[[cntr]])
            MPD_ext[[cntr]] <- mpd(t(ext_mat[[cntr]]),coph_ext)
            MPD_ext[[cntr]][is.na(MPD_ext[[cntr]])] <- 0

    }
    
}

mpd_mean <- sapply(MPD, mean)
mpd_ext_mean <- sapply(MPD_ext, mean)
mpd_mean_change <- mpd_mean - mpd_ext_mean


# Add to data frame
df$mpd_mean <- mpd_mean
df$mpd_ext_mean <- mpd_ext_mean
df$mpd_mean_change <- mpd_mean_change


# mntd
MNTD <- NULL
MNTD_ext <- NULL
MNTD_change <- NULL

for (phy in 1:n_trees){

    coph <- cophenetic(trees[[phy]])

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim
            MNTD[[cntr]] <- mntd(t(original_mat[[cntr]]),coph)
            MNTD[[cntr]][is.na(MNTD[[cntr]])] <- 0

            coph_ext <- cophenetic(trees_extant[[cntr]])
            MNTD_ext[[cntr]] <- mntd(t(ext_mat[[cntr]]),coph_ext)
            MNTD_ext[[cntr]][is.na(MNTD_ext[[cntr]])] <- 0

    }
    
}

mntd_mean <- sapply(MNTD, mean)
mntd_ext_mean <- sapply(MNTD_ext, mean)
mntd_mean_change <- mntd_mean - mntd_ext_mean


# Add to data frame
df$mntd_mean <- mntd_mean
df$mntd_ext_mean <- mntd_ext_mean
df$mntd_mean_change <- mntd_mean_change


# maxD
maxD <- NULL
maxD_ext <- NULL
maxD_change <- NULL

for (phy in 1:n_trees){

    coph <- cophenetic(trees[[phy]])

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim
            maxD[[cntr]] <- maxD(t(original_mat[[cntr]]),coph)
            maxD[[cntr]][is.na(maxD[[cntr]])] <- 0

            coph_ext <- cophenetic(trees_extant[[cntr]])
            maxD_ext[[cntr]] <- maxD(t(ext_mat[[cntr]]),coph_ext)
            maxD_ext[[cntr]][is.na(maxD_ext[[cntr]])] <- 0

    }
    
}

maxD_mean <- sapply(maxD, mean)
maxD_ext_mean <- sapply(maxD_ext, mean)
maxD_mean_change <- maxD_mean - maxD_ext_mean

# Add to data frame
df$maxD_mean <- maxD_mean
df$maxD_ext_mean <- maxD_ext_mean
df$maxD_mean_change <- maxD_mean_change



# # creepjump
# creepjump <- NULL
# creepjump_ext <- NULL
# creepjump_change <- NULL

# for (phy in 1:n_trees){

#     coph <- cophenetic(trees[[phy]])

#     for (sim in 1:n_sims){

#             cntr <- ((phy-1)*n_trees)+sim
#             creepjump[[cntr]] <- creepjump(t(original_mat[[cntr]]),coph)
#             creepjump[[cntr]][is.na(creepjump[[cntr]])] <- 0

#             coph_ext <- cophenetic(trees_extant[[cntr]])
#             creepjump_ext[[cntr]] <- creepjump(t(ext_mat[[cntr]]),coph_ext)
#             creepjump_ext[[cntr]][is.na(creepjump_ext[[cntr]])] <- 0

#     }
    
# }

# creepjump_mean <- sapply(creepjump, mean)
# creepjump_ext_mean <- sapply(creepjump_ext, mean)
# creepjump_mean_change <- creepjump_mean - creepjump_ext_mean


# # Add to data frame
# df$creepjump_mean <- creepjump_mean
# df$creepjump_ext_mean <- creepjump_ext_mean
# df$creepjump_mean_change <- creepjump_mean_change



# pd
PD <- NULL
PD_ext <- NULL
PD_change <- NULL

# mapply(myfun, arg1, arg2)
int_mats_t <- lapply(original_mat,t)
trees_rep <- rep(trees, each=n_sims)
PD <- mapply(pd, int_mats_t, trees_rep, include.root=FALSE, SIMPLIFY=FALSE)
# saveRDS(PD,"dummy_PD.rds")
PD_PD <- sapply(PD, function(x) x[[1]])
PD_PD <- lapply(PD_PD, function(x) { x[is.na(x)] <- 0; x})
PD_SR <- sapply(PD, function(x) x[[2]])

int_mats_t <- lapply(ext_mat,t)
PD_ext <- mapply(pd, int_mats_t, trees_extant, include.root=FALSE, SIMPLIFY=FALSE)
# saveRDS(PD_ext,"dummy_PD_ext.rds")
PD_ext_PD <- sapply(PD_ext, function(x) x[[1]])
PD_ext_PD <- lapply(PD_ext_PD, function(x) { x[is.na(x)] <- 0; x})
PD_ext_SR <- sapply(PD_ext, function(x) x[[2]])

# PD change for extant parasites (infect at least one host currently)
PD_change <- mapply(function(x_PD, y_PD, y_SR) 
                             x_PD[y_SR!=0]-y_PD[y_SR!=0], 
                            PD_PD, PD_ext_PD, PD_ext_SR)

# PD change for extant multi-host parasites (infect at least two hosts currently)
PD_change_ext_multi <- mapply(function(x_PD, y_PD, y_SR) 
                                   x_PD[y_SR<1]-y_PD[y_SR<1], 
                                   PD_PD, PD_ext_PD, PD_ext_SR)


pd_mean <- sapply(PD_PD, mean)
pd_ext_mean <- sapply(PD_ext_PD, mean)
pd_mean_change <- pd_mean - pd_ext_mean

# Add to data frame
df$pd_mean <- pd_mean
df$pd_ext_mean <- pd_ext_mean
df$pd_mean_change <- pd_mean_change



# host extinction measures

# number of hosts lost (matrix level)
nhost <- sapply(original_mat, nrow)
nhost_ext <- sapply(ext_mat, nrow)
nhost_lost <- nhost-nhost_ext
nhost_lost_prop <- nhost_lost/nhost

# Add to data frame
df$nhost <- nhost
df$nhost_ext <- nhost_ext
df$nhost_lost <- nhost_lost


# Evolutionary Distinctivness
ED <- lapply(trees, function(x) evol.distinct(x,type="equal.splits"))
# names(ED)[2] <- "ED"

ED_extant <- lapply(trees_extant, function(x) evol.distinct(x,type="equal.splits"))

ED_gain <- NULL

for (phy in 1:n_trees){

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim
            names(ED_extant[[cntr]])[2] <- "ED_extant"
            ED_gain[[cntr]] <- left_join(ED_extant[[cntr]], ED[[phy]])
            ED_gain[[cntr]]$ED_gain <- ED_gain[[cntr]]$ED_extant - ED_gain[[cntr]]$w
    }
    
}

ed <- sapply(ED, function(x) mean(x$w))
ed_ext <- sapply(ED_extant, function(x) mean(x$ED_extant))
ed_gain_mean <- sapply(ED_gain, function(x) mean(x$ED_gain))


# Add to data frame
df$ED <- rep(ed, each=10)
df$ED_ext <- ed_ext
df$ED_gain_mean <- ed_gain_mean



# Phylogenetic diversity (total for whole matrix)

# Calculating clade matrices per tree
clade_mats <- lapply(trees, clade.matrix)
clade_mats_ext <- lapply(trees_extant, clade.matrix)

pd_tot <- sapply(clade_mats, function(x) pd.calc(x)[1])
pd_ext <- sapply(clade_mats_ext, function(x) pd.calc(x)[1])

pd_loss <- NULL

for (phy in 1:n_trees){

    for (sim in 1:n_sims){

            cntr <- ((phy-1)*n_trees)+sim
            pd_loss[cntr] <- pd_tot[phy] - pd_ext[cntr]
    }
    
}

# Add to data frame
df$pd_tot <- rep(pd_tot, each=10)
df$pd_ext <- pd_ext
df$pd_loss <- pd_loss

return(df)


}














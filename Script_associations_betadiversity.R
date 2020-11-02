library(reshape2)
library(plyr)
library(dplyr)
library(lattice)
library(ggplot2)
library(compare)
library(nlme)
library(classInt)
library(maptools)
library(stringr)
library(rgeos)
library(rgdal)
library(sp)
library(mgcv)
library(data.table)
library(automap)
library(glmnet)
library(raster)
library(viridis)
library(see)
library(netassoc)
library(picante)
library(doParallel)
library(igraph)
library(ggraph)
library(spgwr)
library(entropart)

#### Load data

dataprp <- read.csv2("data_publi.csv")

square_centroid <- read.csv2("square_centroid.csv")

#### Calculating species associations ----

# Function to randomize following make_netassoc_network (Morueta-Holmes et al., 2016). The idea is that a species enters a 
# community with a probability that is proportional to its abundance in the regional pool (the whole dataset). 
randomizeM <- function(mat) { 
  # Compute total species frequencies
  total_freqs <- apply(mat, 2, sum)
  total_freqs <- total_freqs / sum(total_freqs)
  result <- apply(mat, 1, function(X) { 
    # We sample sum(X) individuals from the regional pool. With the probability 
    # of a species entering the community proportional to its abundance in the 
    # pool (the whole dataset)
    picked_species <- sample(seq.int(ncol(mat)), 
                             replace = TRUE, 
                             size = sum(X),
                             prob = total_freqs)
    tabulate(picked_species, nbins = ncol(mat))
  })
  t(result) 
}


# Function for result presentation
tabularize <- function(coocmat, name = "assoc") { 
  # Zero-size matrix
  if(is.null(dim(coocmat))){return(NA)}
  # Add numbers if species names are absent in the original matrix
  if(is.null(colnames(coocmat))){ 
    spnames <- seq.int(ncol(coocmat))
  } else { 
    spnames <- colnames(coocmat)
  }
  # Create new table
  tab <- data.frame(expand.grid(spi = spnames, spj = spnames), 
                    as.vector(coocmat))
  names(tab) <- c('spi', 'spj', name)
  return(tab)
}

# Function to calculate species associations
g<-function(mat, # observed data
            N.null, # number of null matrices to simulate
            occmin){ # minimum occurence number for a species to be selected
  
  mat[,c("code_point","year","habit","zonebio")]<-list(NULL)
  
  mat<-round(mat)
 
  # select species with more occurence than occmin
  matpa <- replace(mat, mat != 0, 1)
  mat <- mat[, colSums(matpa) > occmin]
  name <- names(mat)
  
  # log-transformation of abundance data as recommanded in Morueta-Holmes et al. 2016
  mat <- log(mat+1e-6) - log(1e-6)
  
  # observed association matrix
  mat <- t(as.matrix(mat))
  D.obs <- partial_correlation(mat,"shrinkage")
  
  # random association matrices
  D.null <- array(rep(NA,nrow(mat)*nrow(mat)*N.null),
                  dim = c(nrow(mat),nrow(mat),N.null))

  for(null in 1:N.null) {
    D.null[ , ,null] <- partial_correlation(apply(replicate(1000, t(randomizeM(t(mat)))), c(1, 2), mean),
                                            "shrinkage")
  }
  
  D.null.mean<-apply(D.null, c(1,2), mean)
  D.null.sd<-apply(D.null, c(1,2), sd)
  
  SES<-(D.obs-D.null.mean)/D.null.sd
  
  # pvalue matrix
  M <- array(c(D.obs, D.null), dim=c(nrow(mat),nrow(mat),N.null+1))
  p.val <- apply(M, c(1,2), function(x) {
    if ( (rank(x)[1]/(N.null+1)) < 0.5 ) { 
      rank(x)[1]/(N.null+1)
    } else { 
      ( 1 - rank(x)[1]/(N.null+1) )
    } 
  }) 
  
  p.val<-2*p.val
  diag(p.val)<-NA
  
  # adjusted pvalues
  p.val.adj<-matrix(p.adjust(p.val, method = "BH"), nrow=nrow(p.val), ncol=ncol(p.val))
  
  # select significant SES
  SES[p.val.adj>0.05] <- NA
  SES <- matrix(SES, nrow = nrow(p.val.adj), ncol = ncol(p.val.adj))
  
  # retrieve selected species names
  attr(SES,"dimnames")<-list(name)
  attr(SES,"dimnames")[[2]]<-attr(SES,"dimnames")[[1]]

  results <- data.frame(tabularize(D.obs, name = "obs"), 
                        ses  = tabularize(SES)[ ,3], 
                        pval = tabularize(p.val.adj)[ ,3], 
                        nullmean = tabularize(D.null.mean)[ ,3], 
                        nullsd = tabularize(D.null.sd)[ ,3])
  
  return(results)
}

registerDoParallel(cores = c(detectCores()-1)) # parallelise computation as this might take a while

dataperpoint <- dcast(dataprp,code_point+year+habit+zonebio~code_sp,fun.aggregate = sum,value.var="abond") # use data from all years

matrix_association <- ddply(dataperpoint, .(year, habit, zonebio), .fun=g, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)

matrix_association <- as.data.frame(matrix_association
                                    %>% group_by(habit,zonebio,spi,spj)
                                    %>% summarize(ses=mean(ses, na.rm=T)))

#### Computing community intensity and attractiveness ----

community_association <- function(x, nb.iter){
  hab <- names(which.max(table(droplevels(x)$habit)))
  zb <- names(which.max(table(droplevels(x)$zonebio)))
  aa <- droplevels(subset(matrix_association, habit==hab))
  aa <- droplevels(subset(aa, zonebio==zb))
  
  # calculate indices for the observed data
  b <- droplevels(x) %>% group_by(code_sp) %>% summarize(abond=sum(abond))
  
  datasso <- data.frame(spi=c(t(combn(levels(b$code_sp),2))[,1], t(combn(levels(b$code_sp),2))[,2]), spj=c(t(combn(levels(b$code_sp),2))[,2], t(combn(levels(b$code_sp),2))[,1]))
  aaa <- droplevels(subset(aa, spi %in% datasso$spi))
  aaa <- droplevels(subset(aaa, spj %in% datasso$spj))
  aaa <- aaa[aaa$spi!=aaa$spj,]
  aaa$abond <- b$abond[match(aaa$spi,b$code_sp)]*b$abond[match(aaa$spj,b$code_sp)]/2
  aaa$ses_abond <- aaa$abond*aaa$ses
  
  aaa <- data.frame(ab_ses_abond=sum(abs(aaa$ses_abond), na.rm = TRUE),
                  abond=sum(aaa$abond, na.rm = TRUE),
                  nb_p=length(which(sign(aaa$ses)==1)),
                  nb_n=length(which(sign(aaa$ses)==-1)))
  if((aaa$nb_p+aaa$nb_n)>0){
    result_intensity <- aaa$ab_ses_abond/aaa$abond
    result_attractiveness <- (aaa$nb_p-aaa$nb_n)/(aaa$nb_p+aaa$nb_n)
  }else{
    result_intensity <- 0
    result_attractiveness <- 0
  }
  
  # calculate indices for random data
  data.null <- droplevels(subset(dataprp, code_square==levels(droplevels(x$code_square))))

  for(i in 1: nb.iter){
    data.null2 <- droplevels(data.null) %>% group_by(code_sp) %>% summarize(abond=mean(abond))
    
    data.null2 <- data.null2[sample(nrow(data.null2), length(levels(droplevels(x$code_sp)))), ]
    
    b <- droplevels(data.null2) %>% group_by(code_sp) %>% summarize(abond=sum(abond))
    datasso <- data.frame(spi=c(t(combn(levels(b$code_sp),2))[,1], t(combn(levels(b$code_sp),2))[,2]), spj=c(t(combn(levels(b$code_sp),2))[,2], t(combn(levels(b$code_sp),2))[,1]))
    aaa <- droplevels(subset(aa, spi %in% datasso$spi))
    aaa <- droplevels(subset(aaa, spj %in% datasso$spj))
    aaa$abond <- b$abond[match(aaa$spi,b$code_sp)]*b$abond[match(aaa$spj,b$code_sp)]/2
    aaa$ses_abond <- aaa$abond*aaa$ses
    
    aaa <- data.frame(ab_ses_abond=sum(abs(aaa$ses_abond), na.rm = TRUE),
                    abond=sum(aaa$abond, na.rm = TRUE),
                    nb_p=length(which(sign(aaa$ses)==1)),
                    nb_n=length(which(sign(aaa$ses)==-1)))
    
    if((aaa$nb_p+aaa$nb_n)>0){
      result_intensity[i+1] <- aaa$ab_ses_abond/aaa$abond
      result_attractiveness[i+1] <- (aaa$nb_p-aaa$nb_n)/(aaa$nb_p+aaa$nb_n)
    }else{
      result_intensity[i+1] <- 0
      result_attractiveness[i+1] <- 0
    }
  }
  
  # compute the SES
  
  result <- data.frame(habit=hab,
                     zonebio=zb,
                     intensity=(result_intensity[1]-mean(result_intensity[2:(nb.iter+1)]))/sd(result_intensity[2:(nb.iter+1)]),
                     attractiveness=(result_attractiveness[1]-mean(result_attractiveness[2:(nb.iter+1)]))/sd(result_attractiveness[2:(nb.iter+1)]))
  
  
  return(result)
}

community_association2 <- function(x, nb.iter){tryCatch(community_association(x, nb.iter),
                                                    error=function(e) data.frame(habit="NA", zonebio="NA", intensity=NA,attractiveness=NA))}

com_association <- ddply(droplevels(dataprp), .(code_point, year), .fun=community_association2,.parallel =T, nb.iter=200)


#### Computing community structure ----

a  <-  matrix_association
a$pair <- as.factor(paste0(a$spi,sep="_",a$spj))
pairselect <- as.factor(paste0(levels(a$spi)[combn(1:109,2)[1,]],sep="_",levels(a$spj)[combn(1:109,2)[2,]]))
a <- droplevels(subset(a, pair %in% levels(pairselect)))

deg_comm <- function(x){
  links <- a[,c(3,4,5,1,2)]
  links <- na.omit(links)
  links$ses <- abs(links$ses)
  
  # select association from the right habitat and biogeographic region
  links <- subset(links, links$habit==names(which.max(table(x$habit))))
  links <- subset(links, links$zonebio==names(which.max(table(x$zonebio))))
  colnames(links)[3]  <-  "weight"
  rownames(links)  <-  NULL
  
  # specify the species
  nodes <- x[,c("code_sp","code_point","year")]
  nodes2 <- droplevels(nodes)
  nodes2 <- as.data.frame(nodes2 %>% group_by(code_sp) %>% summarize(count=n()))
  links2 <- subset(links, links$spi %in% levels(droplevels(nodes2$code_sp)))
  links2 <- droplevels(subset(links2, links2$spj %in% levels(droplevels(nodes2$code_sp))))
  
  # make the network
  net  <-  graph_from_data_frame(d=links2, vertices=nodes2, directed=F)
  
  if(nrow(nodes2)>=5){
    mean_clique <- c()
    for(i in 1:100){
      nodes3 <- droplevels(subset(nodes2, code_sp %in% sample(levels(nodes2$code_sp),5)))
      links3 <- subset(links, links$spi %in% levels(droplevels(nodes3$code_sp)))
      links3 <- droplevels(subset(links3, links3$spj %in% levels(droplevels(nodes3$code_sp))))
      net3  <-  graph_from_data_frame(d=links3, vertices=nodes3, directed=F)
      mean_clique[i] <- length(cliques(net3,min=3))/(2^5-(1+(5^2)/2+5/2))
    }
  }else{mean_clique <- NA}

  
  # calculate clique structure
  if(nrow(links2)<=175){
    nb_cliq_3_std_3 <- length(cliques(net,min=3))/(2^length(nodes2$code_sp)-(1+(length(nodes2$code_sp)^2)/2+length(nodes2$code_sp)/2))
  }else{nb_cliq_3_std_3 <- NA}
  
  return(c(nb_cliq_3_std_3, transitivity(net), nrow(nodes2), mean(mean_clique)))
}

deg_comm2 <- function(x){tryCatch(deg_comm(x), error=function(e) c(NA, NA, NA, NA))}

com_structure <- ddply(droplevels(dataprp),.(code_point,year),.parallel=F, .progress = "text",deg_comm2)
com_structure <- merge(com_structure, square_centroid[,-6], by="code_square")

#### Computing regional values ----

# Find best window size (best compromise between the number of windows with more than 20 sites and spatial accuracy) 

test_com_association <- com_association
test_com_association$code_square  <-  substr(test_com_association$code_point, 1, 6)
names(test_com_association)[6:7] <- c("intensity","attractiveness")

test_com_association1 <-  as.data.frame(test_com_association %>% group_by(code_square, year) %>% summarize(intensity=mean(intensity, na.rm=T),
                                                                            attractiveness=mean(attractiveness, na.rm=T)))
test_com_association2 <- as.data.frame(test_com_association %>% group_by(code_square) %>% summarize(count=n()))
test_com_association2$lon <- square_centroid$long[match(test_com_association2$code_square, square_centroid$code_square)]
test_com_association2$lat <- square_centroid$lat[match(test_com_association2$code_square, square_centroid$code_square)]
test_com_association2 <- na.omit(test_com_association2)
testsquarecoord <- test_com_association2[,c("lon","lat")]
testsquarecoord <- SpatialPoints(testsquarecoord,proj4string = CRS("+proj=longlat +datum=WGS84"))
testsquarecoordlamb <- spTransform(testsquarecoord, CRS("+init=epsg:27572"))
matdist <- data.frame(spDists(testsquarecoordlamb))
b <- na.omit(test_com_association1)
b$year <- as.numeric(as.character(b$year))

nb_site_window <- function(i, bandwidth){
  data_dist <- data.frame(code_square=test_com_association2$code_square, dist=matdist[,i])
  data_dist$weight <- gwr.bisquare(data_dist$dist^2, bandwidth)
  data_dist <- droplevels(subset(data_dist, weight!=0))
  b_20 <- droplevels(subset(b, as.factor(b$code_square) %in% levels(as.factor(data_dist$code_square))))
  b_20$weight <- data_dist$weight[match(b_20$code_square, data_dist$code_square)]
  nb_site <- length(levels(as.factor(b_20$code_square)))
  return(nb_site)
}

bandwidths <- seq(from=10000, to=200000, length.out=39)

for(j in 1:length(bandwidths)){
  colnum <- ncol(test_com_association2)
  print(j)
  for(i in 1:length(matdist)){
    #print(i)
    aa <- nb_site_window(i, bandwidths[j])
    test_com_association2[i, colnum+1] <- aa[1]}
  names(test_com_association2)[(colnum+1)] <- c(paste0("scale", j))
}

nb_site_window_result <- data.frame(scale=bandwidths, nb_site=apply(test_com_association2[,5:43], 2, function(x){sum(x>=20)}))

ggplot(nb_site_window_result, aes(x = scale/1000, y = nb_site))+
  geom_point(size=2)+
  theme_light()+
  labs(x ="Window radius (km)", y = "Number of windows with more than 20 sites")+ theme(legend.position="none")+
  geom_vline(xintercept=80, linetype="dashed", size=1)


# Intensity and attractiveness

com_association$code_square  <-  substr(com_association$code_point, 1, 6)
names(com_association)[5:6] <- c("intensity","attractiveness")

com_association1 <-  as.data.frame(com_association %>% group_by(code_square, year) %>% summarize(intensity=mean(intensity, na.rm=T),
                                                                                               attractiveness=mean(attractiveness, na.rm=T)))
com_association2 <- as.data.frame(com_association %>% group_by(code_square) %>% summarize(count=n()))
com_association2$lon <- square_centroid$long[match(com_association2$code_square, square_centroid$code_square)]
com_association2$lat <- square_centroid$lat[match(com_association2$code_square, square_centroid$code_square)]
com_association2 <- na.omit(com_association2)
testsquarecoord <- com_association2[,c("lon","lat")]
testsquarecoord <- SpatialPoints(testsquarecoord,proj4string = CRS("+proj=longlat +datum=WGS84"))
testsquarecoordlamb <- spTransform(testsquarecoord, CRS("+init=epsg:27572"))
matdist <- data.frame(spDists(testsquarecoordlamb))
b <- na.omit(com_association1)
b$year <- as.numeric(as.character(b$year))

asso_an_gwr <- function(i, bandwidth){
  data_dist <- data.frame(code_square=com_association2$code_square, dist=matdist[,i])
  data_dist$weight <- gwr.bisquare(data_dist$dist^2, bandwidth)
  data_dist <- droplevels(subset(data_dist, weight!=0))
  b_20 <- droplevels(subset(b, as.factor(b$code_square) %in% levels(as.factor(data_dist$code_square))))
  b_20$weight <- data_dist$weight[match(b_20$code_square, data_dist$code_square)]
  intensitymoy <- lm(intensity~year, weights = weight, data=b_20)
  attractivenessmoy <- lm(attractiveness~year, weights = weight, data=b_20)
  mean_intensity <- weighted.mean(b_20$intensity, b_20$weight)
  mean_attractiveness <- weighted.mean(b_20$attractiveness, b_20$weight)
  return(c(intensitymoy$coefficients[2],summary(intensitymoy)$coefficients[2,4],
           attractivenessmoy$coefficients[2],summary(attractivenessmoy)$coefficients[2,4],
           mean_intensity, mean_attractiveness))
}


bandwidths <- seq(from=10000, to=200000, length.out=20)

for(j in 1:length(bandwidths)){
  colnum <- ncol(com_association2)
  print(j)
  for(i in 1:length(matdist)){
    #print(i)
    aa <- asso_an_gwr(i, bandwidths[j])
    com_association2[i, colnum+1] <- aa[1]
    com_association2[i, colnum+2] <- aa[2]
    com_association2[i, colnum+3] <- aa[3]
    com_association2[i, colnum+4] <- aa[4]
    com_association2[i, colnum+5] <- aa[5]
    com_association2[i, colnum+6] <- aa[6]}
  names(com_association2)[(colnum+1):(colnum+6)] <- c(paste0("t_intensity", j),paste0("p_t_intensity", j),
                                                    paste0("t_attractiveness", j),paste0("p_t_attractiveness", j),
                                                    paste0("intensity", j),paste0("attractiveness", j))
}


com_association2 <- as.data.frame(com_association2)

# compute mean values along all window sizes

com_association2$t_intensity <- apply(com_association2[,5:124], 1,
                                      FUN=function(x){
                                        z <- x[seq(from=1,to=115,by=6)]
                                        z2 <- x[seq(from=2,to=116,by=6)]
                                        return(mean(z[which(z2<=0.05)],na.rm=T))
                                      })
com_association2$t_attractiveness <- apply(com_association2[,5:124], 1,
                               FUN=function(x){
                                 z <- x[seq(from=3,to=117,by=6)]
                                 z2 <- x[seq(from=4,to=118,by=6)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
com_association2$intensity <- apply(com_association2[,seq(from=9,to=123,by=6)], 1, FUN=function(x){mean(x,na.rm=T)})
com_association2$attractiveness <- apply(com_association2[,seq(from=10,to=124,by=6)], 1, FUN=function(x){mean(x,na.rm=T)})

com_association2 <- merge(com_association2, square_centroid[,c(1,4,5)], by="code_square")

# Clique structure

com_structure$code_square  <-  substr(com_structure$code_point, 1, 6)
names(com_structure)[3:6] <- c("clique","transitivity", "size", "clique_standard")
com_structure1  <-  as.data.frame(com_structure %>% group_by(code_square, year) %>% summarize(clique=mean(clique, na.rm=T),
                                                                                            trans=mean(transitivity, na.rm=T),
                                                                                            size=mean(size),
                                                                                            clique_standard=mean(clique_standard, na.rm=T)))

com_structure1$lat <- square_centroid$lat[match(com_structure1$code_square, square_centroid$code_square)]
com_structure1$lon <- square_centroid$long[match(com_structure1$code_square, square_centroid$code_square)]
com_structure1$lat2 <- square_centroid$lat2[match(com_structure1$code_square, square_centroid$code_square)]
com_structure1$lon2 <- square_centroid$lon2[match(com_structure1$code_square, square_centroid$code_square)]
com_structure2 <- as.data.frame(com_structure1 %>% group_by(code_square, lat, lon) %>% summarize(count=n()))
com_structure2 <- na.omit(com_structure2)
squarecoord <- com_structure2[,c("lon","lat")]
squarecoord <- SpatialPoints(squarecoord,proj4string = CRS("+proj=longlat +datum=WGS84"))
squarecoordlamb <- spTransform(squarecoord, CRS("+init=epsg:27572"))
matdist <- data.frame(spDists(squarecoordlamb))

bandwidths <- seq(from=10000, to=200000, length.out=20)

clique_gwr2 <- function(i, bandwidth){
  data_dist <- data.frame(code_square=com_structure2$code_square, dist=matdist[,i])
  data_dist$weight <- gwr.bisquare(data_dist$dist^2, bandwidth)
  data_dist <- droplevels(subset(data_dist, weight!=0))
  b_20 <- droplevels(subset(com_structure1, com_structure1$code_square %in% levels(as.factor(data_dist$code_square))))
  b_20$weight <- data_dist$weight[match(b_20$code_square, data_dist$code_square)]
  t_clique <- lm(sqrt(clique)~as.numeric(year), data=na.omit(b_20), weights = weight)
  clique <- weighted.mean(b_20$clique, b_20$weight, na.rm=T)
  t_trans <- lm(trans~as.numeric(year), data=na.omit(b_20), weights = weight)
  trans <- weighted.mean(b_20$trans, b_20$weight, na.rm=T)
  t_size <- lm(size~as.numeric(year), data=na.omit(b_20), weights = weight)
  size <- weighted.mean(b_20$size, b_20$weight, na.rm=T)
  t_clique_st <- lm(sqrt(clique_standard)~as.numeric(year), data=na.omit(b_20), weights = weight)
  clique_st <- weighted.mean(b_20$clique_standard, b_20$weight, na.rm=T)
  return(c(t_clique$coefficients[2],summary(t_clique)$coefficients[2,4],clique,
           t_trans$coefficients[2],summary(t_trans)$coefficients[2,4],trans,
           t_size$coefficients[2],summary(t_size)$coefficients[2,4],size,
           t_clique_st$coefficients[2],summary(t_clique_st)$coefficients[2,4],clique_st))
}

clique_gwr3 <- function(i, bandwidth){tryCatch(clique_gwr2(i, bandwidth),
                                             error=function(e) c(NA,NA,NA,
                                                                 NA,NA,NA,
                                                                 NA,NA,NA,
                                                                 NA,NA,NA))}

for(j in 1:length(bandwidths)){
  colnum <- ncol(com_structure2)
  print(j)
  for(i in 1:length(matdist)){
    aa <- clique_gwr3(i, bandwidths[j])
    com_structure2[i, colnum+1] <- aa[1]
    com_structure2[i, colnum+2] <- aa[2]
    com_structure2[i, colnum+3] <- aa[3]
    com_structure2[i, colnum+4] <- aa[4]
    com_structure2[i, colnum+5] <- aa[5]
    com_structure2[i, colnum+6] <- aa[6]
    com_structure2[i, colnum+7] <- aa[7]
    com_structure2[i, colnum+8] <- aa[8]
    com_structure2[i, colnum+9] <- aa[9]
    com_structure2[i, colnum+10] <- aa[10]
    com_structure2[i, colnum+11] <- aa[11]
    com_structure2[i, colnum+12] <- aa[12]
  }
  names(com_structure2)[(colnum+1):(colnum+12)] <- c(paste0("t_clique", j),paste0("p_t_clique", j),paste0("clique", j),
                                                   paste0("t_trans", j),paste0("p_t_trans", j),paste0("trans", j),
                                                   paste0("t_size", j),paste0("p_t_size", j),paste0("size", j),
                                                   paste0("t_clique_st", j),paste0("p_t_clique_st", j),paste0("clique_st", j))
}

com_structure2 <- as.data.frame(com_structure2)

# compute mean values along all window sizes

com_structure2$t_clique <- apply(com_structure2[,5:244], 1,
                                    FUN=function(x){
                                      z <- x[seq(from=1,to=240,by=12)]
                                      z2 <- x[seq(from=2,to=240,by=12)]
                                      return(mean(z[which(z2<=0.05)],na.rm=T))
                                    })
com_structure2$clique <- apply(com_structure2[,seq(from=7,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})
com_structure2$t_clique_ns <- apply(com_structure2[,seq(from=5,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})

com_structure2$t_trans <- apply(com_structure2[,5:244], 1,
                               FUN=function(x){
                                 z <- x[seq(from=4,to=240,by=12)]
                                 z2 <- x[seq(from=5,to=240,by=12)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
com_structure2$trans <- apply(com_structure2[,seq(from=10,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})
com_structure2$t_trans_ns <- apply(com_structure2[,seq(from=8,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})

com_structure2$t_size <- apply(com_structure2[,5:244], 1,
                               FUN=function(x){
                                 z <- x[seq(from=7,to=240,by=12)]
                                 z2 <- x[seq(from=8,to=240,by=12)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
com_structure2$size <- apply(com_structure2[,seq(from=13,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})
com_structure2$t_size_ns <- apply(com_structure2[,seq(from=11,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})

com_structure2$t_clique_st <- apply(com_structure2[,5:244], 1,
                               FUN=function(x){
                                 z <- x[seq(from=10,to=240,by=12)]
                                 z2 <- x[seq(from=11,to=240,by=12)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
com_structure2$clique_st <- apply(com_structure2[,seq(from=16,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})
com_structure2$t_clique_st_ns <- apply(com_structure2[,seq(from=14,to=244,by=12)], 1, FUN=function(x){mean(x,na.rm=T)})

com_structure2 <- merge(com_structure2, square_centroid[,c(1,4,5)], by="code_square")


#### Calculating beta diversity ----

# Beta-diversity

bdiv <- dataprp[,c("code_square","year","code_sp","abond","longitude_wgs84","latitude_wgs84","habit","zonebio")]
bdiv$zonebio <- square_centroid$zonebio[match(bdiv$code_square,square_centroid$code_square)]
habitat_maj <- bdiv %>% group_by(code_square, habit) %>% summarize(count=n()) %>% slice(which.max(count))
habitat_maj_an <- bdiv %>% group_by(code_square, year, habit) %>% summarize(count=n()) %>% slice(which.max(count))
bdiv_c <- dcast(bdiv, code_square+zonebio~code_sp, fun.aggregate = sum,value.var="abond")
bdiv_c$habit <- habitat_maj$habit
bdiv_c_an <- dcast(bdiv, code_square+year+zonebio~code_sp, fun.aggregate = sum,value.var="abond")
bdiv_c_an$habit <- habitat_maj_an$habit

data_div <- bdiv_c[,c(1,2,111)]
data_div$lon <- square_centroid$long[match(data_div$code_square, square_centroid$code_square)]
data_div$lat <- square_centroid$lat[match(data_div$code_square, square_centroid$code_square)]
data_div$lon2 <- square_centroid$lon2[match(data_div$code_square, square_centroid$code_square)]
data_div$lat2 <- square_centroid$lat2[match(data_div$code_square, square_centroid$code_square)]

testsquarecoord <- data_div[,c("lon2","lat2")]
testsquarecoord <- SpatialPoints(testsquarecoord,proj4string = CRS("+init=epsg:27572"))
matdist <- data.frame(spDists(testsquarecoord))
b <- bdiv_c_an

div_an_unif <- function(i, bandwidth){
  data_dist <- data.frame(code_square=data_div$code_square, dist=matdist[,i])
  data_dist <- droplevels(subset(data_dist, dist<bandwidth))
  b_20 <- droplevels(subset(b, as.factor(b$code_square) %in% levels(as.factor(data_dist$code_square))))
  b_c <- droplevels(subset(droplevels(bdiv_c), code_square %in% levels(as.factor(data_dist$code_square))))
  if(nrow(b_c)<10){ # choose the minimum number of square to compute the beta-diversity
    beta_div <- data.frame(rep(NA,10))
    pente <- data.frame(rep(NA,10))
    p.value <- data.frame(rep(NA,10))
  }
  if(nrow(b_c)>=10 & nrow(b_c)<20){
    beta_div <- data.frame(rep(NA,10))
    pente <- data.frame(rep(NA,10))
    p.value <- data.frame(rep(NA,10))
    sous_ech_square <- c(as.character(data_div$code_square[i]),
                      sample(levels(as.factor(data_dist$code_square))[which(levels(as.factor(data_dist$code_square))!=data_div$code_square[i])],9))
    b_20 <- droplevels(subset(b, as.factor(b$code_square) %in% sous_ech_square))
    b_c <- droplevels(subset(droplevels(bdiv_c), code_square %in% sous_ech_square))
    beta_div[1,1] <- BetaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_square,zonebio,habit)))))$Total # calculate the beta-diversity using all years
    beta_div[,1] <- rep(beta_div[1,1], 10)
    bdiv_record <- data.frame()
    for(z in 1:length(levels(droplevels(as.factor(b_20$year))))){
      if(nrow(subset(b_20, year %in% levels(droplevels(as.factor(b_20$year)))[z]))<5){
        bdiv_record[z,1] <- NA
        bdiv_record[z,2] <- levels(droplevels(as.factor(b_20$year)))[z]
      }else{
        bdiv_record[z,1] <- BetaDiversity(MetaCommunity(t(subset(b_20, year %in% levels(droplevels(as.factor(b_20$year)))[z])[,-which(names(b_20) %in% c("code_square","year","zonebio","habit"))])))$Total
        bdiv_record[z,2] <- levels(droplevels(as.factor(b_20$year)))[z]
      }
    } # calculate the beta-diversity for each year
    if(nrow(na.omit(bdiv_record))<2){
      pente[,1] <- rep(NA,10)
      p.value[,1] <- rep(NA,10)
    }else{
      beta_pente <- lm(bdiv_record[,1]~as.numeric(bdiv_record[,2])) # calculate beta-diversity trend
      pente[,1] <-  rep(beta_pente$coefficients[2],10)
      p.value[,1] <- rep(summary(beta_pente)$coefficients[2,4],10)
    }
  }
  if(nrow(b_c)>=20){ # if more than 20 square in the window, calculate several time the beta-diversity on 10 square samples
    beta_div <- data.frame()
    pente <- data.frame()
    p.value <- data.frame()
    for(k in 1:10){
      sous_ech_square <- c(as.character(data_div$code_square[i]),
                        sample(levels(as.factor(data_dist$code_square))[which(levels(as.factor(data_dist$code_square))!=data_div$code_square[i])],9))
      b_20 <- droplevels(subset(b, as.factor(b$code_square) %in% sous_ech_square))
      b_c <- droplevels(subset(droplevels(bdiv_c), code_square %in% sous_ech_square))
      beta_div[k,1] <- BetaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_square,zonebio,habit)))))$Total
      bdiv_record <- data.frame()
      for(z in 1:length(levels(droplevels(as.factor(b_20$year))))){
        if(nrow(subset(b_20, year %in% levels(droplevels(as.factor(b_20$year)))[z]))<5){
          bdiv_record[z,1] <- NA
          bdiv_record[z,2] <- levels(droplevels(as.factor(b_20$year)))[z]
        }else{
          bdiv_record[z,1] <- BetaDiversity(MetaCommunity(t(subset(b_20, year %in% levels(droplevels(as.factor(b_20$year)))[z])[,-which(names(b_20) %in% c("code_square","year","zonebio","habit"))])))$Total
          bdiv_record[z,2] <- levels(droplevels(as.factor(b_20$year)))[z]
        }
      }
      if(nrow(na.omit(bdiv_record))<2){
        pente[k,1] <- NA
        p.value[k,1] <- NA
      }else{
        beta_pente <- lm(bdiv_record[,1]~as.numeric(bdiv_record[,2]))
        pente[k,1] <- beta_pente$coefficients[2]
        p.value[k,1] <- summary(beta_pente)$coefficients[2,4]
      }
      
    }
  }
  return(data.frame(beta_div, pente, p.value))
}


bandwidths <- seq(from=10000, to=200000, length.out=20)
for(j in 1:length(bandwidths)){
  colnum <- ncol(data_div)
  print(j)
  data_div[,(colnum+1):(colnum+30)] <- NA
  for(i in 1:length(matdist)){
    print(i)
    aa <- div_an_unif(i, bandwidths[j])
    data_div[i, (colnum+1):(colnum+10)] <- aa[,1]
    data_div[i, (colnum+11):(colnum+20)] <- aa[,2]
    data_div[i, (colnum+21):(colnum+30)] <- aa[,3]}
  names(data_div)[c((colnum+1),(colnum+11),(colnum+21))] <- c(paste0("beta_div", j), paste0("beta_p", j), paste0("p_val", j))
}


data_div2 <- data_div[,1:7]
for(i in 1:60){
  data_div2[,i+7] <- rowMeans(data_div[,(8+10*(i-1)):(7+10*i)], na.rm=T)
  names(data_div2)[i+7] <- names(data_div)[(8+10*(i-1))]
}
data_div2$beta_div <- apply(data_div2[,seq(from=8,to=65,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})
data_div2$beta_p <- apply(data_div2[,8:67], 1,
                               FUN=function(x){
                                 z <- x[seq(from=2,to=59,by=3)]
                                 z2 <- x[seq(from=3,to=60,by=3)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
data_div2$beta_p_ns <- apply(data_div2[,seq(from=9,to=66,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})

#### Testing relationships ----

a <- merge(com_association2, com_structure2[,-c(2,3,4,257,258)], by.x=c("code_square"), by.y=c("code_square"))
a <- merge(a, data_div2[,-c(2:7)],by.x=c("code_square"), by.y=c("code_square"))

rela <- data.frame(cor=c(summary(gam(a[,128]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,128]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20))) # attractiveness vs intensity
for(i in 1:20){
  rela[i+1,1] <- summary(gam(a[,6*i+4]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela[i+1,2] <- summary(gam(a[,6*i+4]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela2 <- data.frame(cor=c(summary(gam(a[,127]~a[,372]+a[,379]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,127]~a[,372]+a[,379]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20))) # clique vs intensity
for(i in 1:20){
  rela2[i+1,1] <- summary(gam(a[,6*i+3]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela2[i+1,2] <- summary(gam(a[,6*i+3]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela3 <- data.frame(cor=c(summary(gam(a[,128]~a[,372]+a[,379]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,128]~a[,372]+a[,379]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20))) # attractiveness vs clique
for(i in 1:20){
  rela3[i+1,1] <- summary(gam(a[,6*i+4]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela3[i+1,2] <- summary(gam(a[,6*i+4]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_sp <- as.data.frame(cbind(scale(rela[,1],center = F), rela[,2:3], scale(rela2[,1],center = F),rela2[,2], scale(rela3[,1],center = F),rela3[,2]))
names(rela_sp) <- c("attractiveness/intensity","pval1","scale","clique/intensity","pval2","clique/attractiveness","pval3")
rela_sp1 <- melt(rela_sp, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_sp2 <- melt(rela_sp, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_sp1$pvalue <-  -sign(rela_sp2$value-0.05)
ggplot(rela_sp1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Window radius (km)", y = "Slope")+ theme(legend.position="none")

rela4 <- data.frame(cor=c(summary(gam(a[,126]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,126]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 10:20){
  aa <- a[which(a[,6*i]<=0.05 & a[,6*i+2]<=0.05),]
  rela4[i+1,1] <- summary(gam(aa[,6*i+1]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # attractiveness trend vs intensity trend
  rela4[i+1,2] <- summary(gam(aa[,6*i+1]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela5 <- data.frame(cor=c(summary(gam(a[,125]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,125]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 7:20){
  aa <- a[which(a[,6*i]<=0.05 & a[,12*i+120]<=0.05),]
  rela5[i+1,1] <- summary(gam(aa[,6*i-1]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # clique trend vs intensity trend
  rela5[i+1,2] <- summary(gam(aa[,6*i-1]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela6 <- data.frame(cor=c(summary(gam(a[,126]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,126]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 8:20){
  aa <- a[which(a[,6*i+2]<=0.05 & a[,12*i+120]<=0.05),]
  rela6[i+1,1] <- summary(gam(aa[,6*i+1]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # clique trend vs attractiveness trend
  rela6[i+1,2] <- summary(gam(aa[,6*i+1]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_trend <- as.data.frame(cbind(scale(rela4[,1],center = F), rela4[,2:3], scale(rela5[,1],center = F),rela5[,2], scale(rela6[,1],center = F),rela6[,2]))
names(rela_trend) <- c("attractiveness/intensity","pval1","scale","clique/intensity","pval2","clique/attractiveness","pval3")
rela_trend1 <- melt(rela_trend, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_trend2 <- melt(rela_trend, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_trend1$pvalue <-  -sign(rela_trend2$value-0.05)
ggplot(rela_trend1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Window radius (km)", y = "Slope")+ theme(legend.position="none")


rela_bv1 <- data.frame(cor=c(summary(gam(a[,443]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,443]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20))) # beta-diversity vs intensity
for(i in 2:20){
  rela_bv1[i+1,1] <- summary(gam(a[,3*i+380]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv1[i+1,2] <- summary(gam(a[,3*i+380]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bv2 <- data.frame(cor=c(summary(gam(a[,443]~a[,128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                     pval=c(summary(gam(a[,443]~a[,128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                     echelle=c(0,seq(from=10,to=200,length.out=20))) # beta-diversity vs attractiveness
for(i in 2:20){
  rela_bv2[i+1,1] <- summary(gam(a[,3*i+380]~a[,6*i+4]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv2[i+1,2] <- summary(gam(a[,3*i+380]~a[,6*i+4]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bv3 <- data.frame(cor=c(summary(gam(a[,443]~a[,372]+a[,378]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                     pval=c(summary(gam(a[,443]~a[,372]+a[,378]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                     echelle=c(0,seq(from=10,to=200,length.out=20))) # beta-diversity vs clique
for(i in 2:20){
  rela_bv3[i+1,1] <- summary(gam(a[,3*i+380]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv3[i+1,2] <- summary(gam(a[,3*i+380]~a[,12*i+121]+a[,12*i+127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_betadiv <- as.data.frame(cbind(scale(rela_bv1[,1],center = F), rela_bv1[,2:3], scale(rela_bv2[,1],center = F),rela_bv2[,2], scale(rela_bv3[,1],center = F),rela_bv3[,2]))
names(rela_betadiv) <- c("intensity","pval1","scale","attractiveness","pval2","clique","pval3")
rela_betadiv1 <- melt(rela_betadiv, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_betadiv2 <- melt(rela_betadiv, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_betadiv1$pvalue <-  -sign(rela_betadiv2$value-0.05)
ggplot(rela_betadiv1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Window radius (km)", y = "Slope")+ theme(legend.position="none")


rela_bvt1 <- data.frame(cor=c(summary(gam(a[,445]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,445]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa <- a[a[,6*i]<=0.05,]
  rela_bvt1[i+1,1] <- summary(gam(aa[,3*i+381]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # beta-diversity trend vs intensity trend
  rela_bvt1[i+1,2] <- summary(gam(aa[,3*i+381]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bvt2 <- data.frame(cor=c(summary(gam(a[,445]~a[,126]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,445]~a[,126]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa <- a[a[,6*i+2]<=0.05,]
  rela_bvt2[i+1,1] <- summary(gam(aa[,3*i+381]~aa[,6*i+1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # beta-diversity trend vs attractiveness trend
  rela_bvt2[i+1,2] <- summary(gam(aa[,3*i+381]~aa[,6*i+1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bvt3 <- data.frame(cor=c(summary(gam(a[,445]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,445]~a[,371]+a[,377]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 5:20){
  aa <- a[a[,3*i+129]<=0.05,]
  rela_bvt3[i+1,1] <- summary(gam(aa[,3*i+381]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] # beta-diversity trend vs clique trend
  rela_bvt3[i+1,2] <- summary(gam(aa[,3*i+381]~aa[,12*i+119]+aa[,12*i+125]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_betadiv_t <- as.data.frame(cbind(scale(rela_bvt1[,1],center = F), rela_bvt1[,2:3], scale(rela_bvt2[,1],center = F),rela_bvt2[,2], scale(rela_bvt3[,1],center = F),rela_bvt3[,2]))
names(rela_betadiv_t) <- c("intensity","pval1","scale","attractiveness","pval2","clique","pval3")
rela_betadiv_t1 <- melt(rela_betadiv_t, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_betadiv_t2 <- melt(rela_betadiv_t, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_betadiv_t1$pvalue <-  -sign(rela_betadiv_t2$value-0.05)
ggplot(rela_betadiv_t1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Window radius (km)", y = "Slope")+ theme(legend.position="none")



# Plot residuals

data_to_plot <- data.frame(res=c(residuals(gam(a[,"beta_div8"]~a[,"intensity8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="working"),
                               residuals(gam(a[,"beta_div8"]~a[,"attractiveness8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="working"),
                               residuals(gam(a[,"beta_div8"]~a[,"clique8"]+a[,"size8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="working"),
                               residuals(gam(a[a$p_t_intensity8<0.05,"beta_p8"]~a[a$p_t_intensity8<0.05,"t_intensity8"]+te(a$lon2[a$p_t_intensity8<0.05],a$lat2[a$p_t_intensity8<0.05],bs="tp",k=3)),type="working"),
                               residuals(gam(a[a$p_t_attractiveness8<0.05,"beta_p8"]~a[a$p_t_attractiveness8<0.05,"t_attractiveness8"]+te(a$lon2[a$p_t_attractiveness8<0.05],a$lat2[a$p_t_attractiveness8<0.05],bs="tp",k=3)),type="working"),
                               residuals(gam(a[a$p_t_clique8<0.05,"beta_p8"]~a[a$p_t_clique8<0.05,"t_clique8"]+a[a$p_t_clique8<0.05,"t_size8"]+te(a$lon2[a$p_t_clique8<0.05],a$lat2[a$p_t_clique8<0.05],bs="tp",k=3)),type="working")),
                         predict=c(predict(gam(a[,"beta_div8"]~a[,"intensity8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="terms")[,1],
                                   predict(gam(a[,"beta_div8"]~a[,"attractiveness8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="terms")[,1],
                                   predict(gam(a[,"beta_div8"]~a[,"clique8"]+a[,"size8"]+te(a$lon2,a$lat2,bs="tp",k=3)),type="terms")[,1],
                                   predict(gam(a[a$p_t_intensity8<0.05,"beta_p8"]~a[a$p_t_intensity8<0.05,"t_intensity8"]+te(a$lon2[a$p_t_intensity8<0.05],a$lat2[a$p_t_intensity8<0.05],bs="tp",k=3)),type="terms")[,1],
                                   predict(gam(a[a$p_t_attractiveness8<0.05,"beta_p8"]~a[a$p_t_attractiveness8<0.05,"t_attractiveness8"]+te(a$lon2[a$p_t_attractiveness8<0.05],a$lat2[a$p_t_attractiveness8<0.05],bs="tp",k=3)),type="terms")[,1],
                                   predict(gam(a[a$p_t_clique8<0.05,"beta_p8"]~a[a$p_t_clique8<0.05,"t_clique8"]+a[a$p_t_clique8<0.05,"t_size8"]+te(a$lon2[a$p_t_clique8<0.05],a$lat2[a$p_t_clique8<0.05],bs="tp",k=3)),type="terms")[,1]),
                         term=c(a[which(!is.na(a[,"beta_div8"])),"intensity8"],a[which(!is.na(a[,"beta_div8"])),"attractiveness8"],a[which(!is.na(a[,"beta_div8"])),"clique8"],
                                unlist(gam(a[a$p_t_intensity8<0.05,"beta_p8"]~a[a$p_t_intensity8<0.05,"t_intensity8"]+te(a$lon2[a$p_t_intensity8<0.05],a$lat2[a$p_t_intensity8<0.05],bs="tp",k=3))$model[2]),
                                unlist(gam(a[a$p_t_attractiveness8<0.05,"beta_p8"]~a[a$p_t_attractiveness8<0.05,"t_attractiveness8"]+te(a$lon2[a$p_t_attractiveness8<0.05],a$lat2[a$p_t_attractiveness8<0.05],bs="tp",k=3))$model[2]),
                                unlist(gam(a[a$p_t_clique8<0.05,"beta_p8"]~a[a$p_t_clique8<0.05,"t_clique8"]+a[a$p_t_clique8<0.05,"t_size8"]+te(a$lon2[a$p_t_clique8<0.05],a$lat2[a$p_t_clique8<0.05],bs="tp",k=3))$model[2])),
                         rel=c(rep("intensity",length(a[which(!is.na(a[,"beta_div8"])),"intensity8"])),
                               rep("attractiveness",length(a[which(!is.na(a[,"beta_div8"])),"attractiveness8"])),
                               rep("clique",length(a[which(!is.na(a[,"beta_div8"])),"clique8"])),
                               rep("t_intensity",length(gam(a[a$p_t_intensity8<0.05,"beta_p8"]~a[a$p_t_intensity8<0.05,"t_intensity8"]+te(a$lon2[a$p_t_intensity8<0.05],a$lat2[a$p_t_intensity8<0.05],bs="tp",k=3))$residuals)),
                               rep("t_attractiveness",length(gam(a[a$p_t_attractiveness8<0.05,"beta_p8"]~a[a$p_t_attractiveness8<0.05,"t_attractiveness8"]+te(a$lon2[a$p_t_attractiveness8<0.05],a$lat2[a$p_t_attractiveness8<0.05],bs="tp",k=3))$residuals)),
                               rep("t_clique",length(gam(a[a$p_t_clique8<0.05,"beta_p8"]~a[a$p_t_clique8<0.05,"t_clique8"]+a[a$p_t_clique8<0.05,"t_size8"]+te(a$lon2[a$p_t_clique8<0.05],a$lat2[a$p_t_clique8<0.05],bs="tp",k=3))$residuals))))

# partial residuals for a term are the whole model residuals + the corresponding estimate of the term
data_to_plot$partial_res <- data_to_plot$res+data_to_plot$predict

ggplot(data_to_plot[data_to_plot$rel %in% c("intensity","attractiveness","clique"),], aes(x = term, y = partial_res)) +
  theme_light() +
  geom_point(shape = 1, col = "blue", size = 2) +
  geom_smooth(method = "lm",formula = y~x, se = FALSE)+ facet_grid(. ~ rel, scales="free_x")

ggplot(data_to_plot[which(data_to_plot$rel %in% c("t_intensity","t_attractiveness","t_clique") & data_to_plot$partial_res<0.25),], aes(x = term, y = partial_res)) +
  theme_light() +
  geom_point(shape = 1, col = "blue", size = 2) +
  geom_smooth(method = "lm",formula = y~x, se = FALSE,
              data=data_to_plot[which(data_to_plot$rel %in% c("t_intensity","t_attractiveness") & data_to_plot$partial_res<0.25),])+ facet_grid(. ~ rel, scales="free_x")

#### Mapping indices ----

prep_krig <- function(sub_data,     # data to interpolate
                    accuracy=200){  # accuracy of the output
  names(sub_data) <- c("lon","lat","variable")
  sub_data <- na.omit(sub_data)
  coord_vars <- c("lat","lon")
  data_vars <- setdiff(colnames(sub_data), coord_vars)
  sp_points <- SpatialPoints(sub_data[,coord_vars])
  sp_df <- SpatialPointsDataFrame(sp_points, sub_data[,data_vars,drop=FALSE])
  pixels_per_side <- accuracy
  bottom.left <- apply(sp_points@coords,2,min)
  top.right <- apply(sp_points@coords,2,max)
  margin <- abs((top.right-bottom.left))/10
  bottom.left <- bottom.left-margin
  top.right <- top.right+margin
  pixel.size <- abs(top.right-bottom.left)/pixels_per_side
  g <- GridTopology(cellcentre.offset=bottom.left, cellsize=pixel.size, cells.dim=c(pixels_per_side,pixels_per_side))
  map_base_data <- map_data("france")
  map_base_data <- map_base_data[!map_base_data$region %in% c("Corse du Sud","Haute-Corse"), ]
  foo = function(x) {
    group = unique(x$group)
    Polygons(list(Polygon(x[,c("lat","long")])),ID=group)
  }
  state_pg <- SpatialPolygons(dlply(map_base_data, .(group), foo))
  grid_points <- SpatialPoints(g)
  in_points <- !is.na(over(grid_points,state_pg))
  fit_points <- SpatialPoints(as.data.frame(grid_points)[in_points,])
  krig <- autoKrige(variable~1, sp_df, new_data=fit_points, block = c(0.5,0.5))
  interp_data <- as.data.frame(krig$krige_output)
  colnames(interp_data) <- c("lat","lon","pred","var","stdev")
  return(interp_data)
}

sub_data <- com_association2[, c("lon","lat","intensity")]
interp_data1 <- prep_krig(sub_data)

sub_data <- com_association2[, c("lon","lat","attractiveness")]
interp_data2 <- prep_krig(sub_data)

sub_data <- na.omit(com_association2[, c("lon","lat","t_intensity")])
sub_data$t_intensity[abs(sub_data$t_intensity)>0.3] <- 0.3*sign(sub_data$t_intensity[abs(sub_data$t_intensity)>0.3])
interp_data3 <- prep_krig(sub_data)

sub_data <- na.omit(com_association2[, c("lon","lat","t_attractiveness")])
sub_data$t_attractiveness[abs(sub_data$t_attractiveness)>0.2] <- 0.2*sign(sub_data$t_attractiveness[abs(sub_data$t_attractiveness)>0.2])
interp_data4 <- prep_krig(sub_data)

sub_data <- com_structure2[, c("lon","lat","clique")]
sub_data$clique[sub_data$clique>0.6] <- 0.6
interp_data5 <- prep_krig(sub_data)

sub_data <- na.omit(com_structure2[, c("lon","lat","t_clique")])
sub_data$t_clique[abs(sub_data$t_clique)>0.2] <- 0.2*sign(sub_data$t_clique[abs(sub_data$t_clique)>0.2])
interp_data6 <- prep_krig(sub_data)

sub_data <- com_structure2[, c("lon","lat","clique_st")]
interp_data5b <- prep_krig(sub_data)

sub_data <- na.omit(com_structure2[, c("lon","lat","t_clique_st")])
sub_data$t_clique_st[abs(sub_data$t_clique_st)>0.05] <- 0.05*sign(sub_data$t_clique_st[abs(sub_data$t_clique_st)>0.05])
interp_data6b <- prep_krig(sub_data)

sub_data <- com_structure2[, c("lon","lat","size")]
interp_data5c <- prep_krig(sub_data)

sub_data <- na.omit(com_structure2[, c("lon","lat","t_size")])
interp_data6c <- prep_krig(sub_data)

sub_data <- data_div2[, c("lon","lat","beta_div")]
interp_data7 <- prep_krig(sub_data)

sub_data <- data_div2[, c("lon","lat","beta_p")]
interp_data8 <- prep_krig(sub_data)
interp_data8b <- prep_krig(data_div2[, c("lon","lat","beta_p_ns")])

border_fr <- readOGR("ne_10m_admin_0_map_units.shp") # https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-details/
border_fr <- border_fr[which(border_fr$NAME=="France"),]
border_fr2 <- crop(border_fr, extent(-5.14, 8.5, 42, 51.1))
border_fr3 <-  spTransform(border_fr2, CRS("+init=epsg:27572"))

ggplot(data=interp_data1, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data2, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data3, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data3$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data4, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data4$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data5, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data6, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data6$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data5b, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data6b, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data6$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data5c, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data6c, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data6$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data7, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

ggplot(data=interp_data8, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=pred),color=NA) +
  scale_fill_gradient2(low="#800080",mid="#F7F7F7",high="#FFC100", midpoint=mean(interp_data8$pred)) +
  coord_fixed(1.3)+
  geom_point(data=sub_data,color="black",size=0.3)+theme_void()+borders(border_fr2)

# Correlation maps

corr_map  <-  function(map1, map2){
  equivalent_point <- SpatialPoints(interp_data1[,c(2,1)])
  couche_1 <- rasterFromXYZ(map1[,c(2,1,3)])
  couche_2 <- rasterFromXYZ(map2[,c(2,1,3)])
  map1 <- data.frame(interp_data1[,1:2],extract(couche_1, equivalent_point))
  map2 <- data.frame(interp_data1[,1:2],extract(couche_2, equivalent_point))
  couche_1 <- rasterFromXYZ(map1[,c(2,1,3)])
  couche_2 <- rasterFromXYZ(map2[,c(2,1,3)])
  couche_3 <-  stack(couche_1, couche_2)
  names(couche_3)  <-  c("1", "2")
  couche_corr <- raster(couche_3, 1)
  values(couche_corr)  <-  1:ncell(couche_1)
  
  for(i in 1:20){
    print(i)
    weight.w <- matrix(1,2*i+1,2*i+1)
    focal_cor  <-  focal(
      x = couche_corr,
      pad=T,
      w=weight.w,
      fun = function(x, y = couche_3){
        cor(values(y)[x,1], values(y)[x,2],
            use = "na.or.complete")
      }
    )
    
    focal_cor  <-  focal_cor+couche_1-couche_1
    if(i<2){focal_cor_dat <- focal_cor}
    else{focal_cor_dat <- stack(focal_cor_dat,focal_cor)}
  }
  
  focal_1 <- focal_cor_dat
  
  focal_cor_dat <- na.omit(as.data.frame(mean(focal_1) ,xy=T))
  out_map <- ggplot(focal_cor_dat) +
    geom_tile(aes(x, y, fill = layer)) +
    scale_fill_gradient2("Corr",
                         low = "#800080",
                         mid = "#F7F7F7",
                         high = "#FFC100",
                         midpoint = 0) +
    coord_quickmap(xlim = range(focal_cor_dat$x), ylim = range(focal_cor_dat$y))+
    coord_fixed(1.3) + xlab("") + ylab("")+theme_void()+borders(border_fr2)
  return(out_map)
}

# examples:

int_att <- corr_map(interp_data1,interp_data2)
int_cli <- corr_map(interp_data1,interp_data5)
att_cli <- corr_map(interp_data2,interp_data5)
t_int_att <- corr_map(interp_data3,interp_data4)
t_int_cli <- corr_map(interp_data3,interp_data6)
t_att_cli <- corr_map(interp_data4,interp_data6)

bv_int <- corr_map(interp_data7,interp_data1)
bv_att <- corr_map(interp_data7,interp_data2)
bv_cli <- corr_map(interp_data7,interp_data5)
t_bv_int <- corr_map(interp_data8,interp_data3)
t_bv_att <- corr_map(interp_data8,interp_data4)
t_bv_cli <- corr_map(interp_data8,interp_data6)


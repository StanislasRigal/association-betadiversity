#### Data selection ----

# Raw data
data1<-read.csv2("data_2016.csv")

# Column selection
data2<-data1[, c(1,8:15,20,23:26,29,30,36:40)]

# Build a habitat column 
data3<-data.frame(data2, habit=paste(data2$p_milieu,data2$p_type,sep="_"))
data4<-data3[, -c(20,21)]

# Drop habitat classes
data4$habit[data4$habit %in% c("A_4")]<-"A_1"
data4$habit[data4$habit %in% c("B_3","B_7","B_8")]<-"B_2"
data4$habit[data4$habit %in% c("B_4","B_5","B_6")]<-"B_3"
data4$habit[data4$habit %in% c("C_5")]<-"C_1"
data4$habit[data4$habit %in% c("C_3")]<-"C_2"
data4$habit[data4$habit %in% c("C_6","C_7","C_8","C_9","C_10")]<-"C_4"
data4$habit[data4$habit %in% c("D_6")]<-"D_5"
data4$habit[data4$habit %in% c("F_1","F_2","F_3","F_4","F_5","F_6","F_7","F_8","F_9","F_10","F_11","F_12")]<-"F_1"
data4$habit[data4$habit %in% c("G_2","G_3","G_4","G_5","G_6","G_8")]<-"G_1"
data4$habit[data4$habit %in% c("A_5","B_NA","C_NA","D_7","D_NA","E_5","E_6","E_7","E_NA","P_2","X_1")]<-"NA_NA"
data4$habit<-droplevels(data4)$habit

# Use the maximal abundance value for each monitoring season (two monitorings each season)
data5<-data4 %>% group_by(code_point, code_sp, annee) %>% slice(which.max(abond))
data6<-ungroup(data5)

# use species code ("code_sp") for simplicity purpose and keep the equivalence between code and species name
species_name<-data6[,c("code_sp","scientific_name","french_name", "english_name")]
species_name<-species_name %>% arrange(code_sp) %>% group_by(code_sp,scientific_name,french_name, english_name) %>% summarize(count=n())
species_name<-species_name[,-c(5)]

# remove single square and points
tablepoint<-table(data6$code_point)
tablecarre<-table(as.factor(data6$code_carre))
data6<-subset(data6, !(code_carre %in% c(attr(which(tablecarre==1),"names"))))
data6<-subset(data6, !(code_point %in% c(attr(which(tablepoint==1),"names"))))

# remove points monitored only 1 year
data8<-dcast(data6,code_point+annee+habit~code_sp,fun.aggregate = sum,value.var="abond")
data8<-arrange(data8,code_point,annee,habit)
count1<-data8 %>% group_by(code_point, annee) %>% summarize(count = n())
count1<-count1 %>% group_by(code_point) %>% summarize(count = n())
count1<-droplevels(subset(count1, count==1))
data13<-subset(data6, !(code_point %in% count1$code_point))#data avec espece en observation

# remove square monitored only 1 year
data12<-dcast(data6,code_carre+annee+habit~code_sp,fun.aggregate = sum,value.var="abond")
count1<-data12 %>% group_by(code_carre, annee) %>% summarize(count = n())
count1<-count1 %>% group_by(code_carre) %>% summarize(count = n())
count1<-droplevels(subset(count1, count==1))
data12<-subset(data12, !(code_carre %in% count1$code_carre))
data11<-subset(data6, !(code_carre %in% count1$code_carre))

# select species representing more than 99% of the total abundance

sommeligne <-data12[ , c(4:245)]
x <- rowSums(sommeligne)
data12$somme <- x

sommecol<-data12[ , c(4:245)]
y <- colSums(sommecol)
y<-as.table(y)
y<-y[order(y, decreasing = TRUE)]
barplot(y)
total<- sum(y)
ypour<-(y/total)*100

pourcent<-0
i<-1
while(pourcent<99){
  pourcent<-pourcent+ypour[i]
  i<-i+1}
ypour99<-ypour[1:i]
y99<-y[1:i]
barplot(ypour99)

dataprp<-subset(data13, (code_sp %in% c(attr(y99,"names"))))#sans les especes peu abondantes et sans les points echantillonnes une annee
dataprc<-subset(data11, (code_sp %in% c(attr(y99,"names"))))#sans les especes peu abondantes et sans les carres echantillonnes une annee
dataprp$abond<-round(dataprp$abond)
dataprc$abond<-round(dataprc$abond) # make sure that abundances are integers

# Location data

carre<-readOGR("carres-franceshp.shp")
proj4string(carre)<-CRS("+init=epsg:27572")
carre2<-spTransform(carre,CRS("+proj=longlat +ellps=WGS84"))
carre2$NUMNAT<-str_pad(carre2$NUMNAT, 6, pad = "0")
carre2$region<-substring(carre2$NUMNAT, 1, 2)

carre_centroid<-data.frame(code_carre=carre2$NUMNAT, long=rep(NA, nrow(carre2)), lat=rep(NA, nrow(carre2)),
                           lon2=rep(NA, nrow(carre)), lat2=rep(NA, nrow(carre)))
for(i in 1:nrow(carre2)){
  carre_centroid[i,2:3]<-carre2@polygons[[i]]@labpt
  carre_centroid[i,4:5]<-carre@polygons[[i]]@labpt
}

# Biogeographic information from https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3

france<-readOGR("BiogeoRegions2016.shp")
france <- spTransform(france, CRS(proj4string(carre2)))
france2<-crop(france, extent(carre2))

zonebio<-over(gCentroid(carre2,byid=TRUE), france2 , fn = NULL)
carre_centroid$zonebio<-zonebio$short_name
carre_centroid$zonebio<-droplevels(carre_centroid)$zonebio


#### Calculating species associations ----

library(netassoc)
library(picante)
library(plyr)


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
  
  mat[,c("code_point","longitude_wgs84","latitude_wgs84","habit","zonebio")]<-list(NULL)
  
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

library(doParallel)
registerDoParallel(cores = 16)

dataperpoint<-dcast(dataprp,code_point+longitude_wgs84+latitude_wgs84+habit~code_sp,fun.aggregate = sum,value.var="abond") # use data from all years

matrixinteraction <- ddply(dataperpoint, .(habit, zonebio), .fun=g, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)


#### Computing community intensity and attractiveness ----

community_association<-function(x, nb.iter){
  hab<-names(which.max(table(droplevels(x)$habit)))
  zb<-names(which.max(table(droplevels(x)$zonebio)))
  aa<-droplevels(subset(matrixinteraction_final, habit==hab))
  aa<-droplevels(subset(aa, zonebio==zb))
  
  # calculate indices for the observed data
  b<-droplevels(x) %>% group_by(code_sp) %>% summarize(abond=sum(abond))
  
  datasso<-data.frame(spi=c(t(combn(levels(b$code_sp),2))[,1], t(combn(levels(b$code_sp),2))[,2]), spj=c(t(combn(levels(b$code_sp),2))[,2], t(combn(levels(b$code_sp),2))[,1]))
  aaa<-droplevels(subset(aa, spi %in% datasso$spi))
  aaa<-droplevels(subset(aaa, spj %in% datasso$spj))
  aaa$abond<-b$abond[match(aaa$spi,b$code_sp)]
  aaa$ses_abond<-aaa$abond*aaa$ses
  
  aaa<-data.frame(ses=sum(aaa$ses),
                  sd=sd(aaa$ses, na.rm = TRUE),
                  ab_ses=sum(abs(aaa$ses)),
                  ses_abond=sum(aaa$ses_abond),
                  ab_ses_abond=sum(abs(aaa$ses_abond)),
                  abond=sum(aaa$abond),
                  nb_p=length(which(sign(aaa$ses)==1)),
                  nb_n=length(which(sign(aaa$ses)==-1)),
                  tot=nrow(aaa))
  if((aaa$nb_p+aaa$nb_n)>0){
    result_intensity<-aaa$ab_ses_abond/aaa$abond
    result_attractiveness<-(aaa$nb_p-aaa$nb_n)/(aaa$nb_p+aaa$nb_n)
  }else{
    result_intensity<-0
    result_attractiveness<-0
  }
  
  # calculate indices for random data
  data.null<-droplevels(subset(dataprp, code_carre==levels(droplevels(x$code_carre))))

  for(i in 1: nb.iter){
    data.null2<-droplevels(data.null) %>% group_by(code_sp) %>% summarize(abond=mean(abond))
    
    data.null2<-data.null2[sample(nrow(data.null2), length(levels(droplevels(x$code_sp)))), ]
    
    b<-droplevels(data.null2) %>% group_by(code_sp) %>% summarize(abond=sum(abond))
    datasso<-data.frame(spi=c(t(combn(levels(b$code_sp),2))[,1], t(combn(levels(b$code_sp),2))[,2]), spj=c(t(combn(levels(b$code_sp),2))[,2], t(combn(levels(b$code_sp),2))[,1]))
    aaa<-droplevels(subset(aa, spi %in% datasso$spi))
    aaa<-droplevels(subset(aaa, spj %in% datasso$spj))
    aaa$abond<-b$abond[match(aaa$spi,b$code_sp)]
    aaa$ses_abond<-aaa$abond*aaa$ses
    
    aaa<-data.frame(ses=sum(aaa$ses),
                    sd=sd(aaa$ses, na.rm = TRUE),
                    ab_ses=sum(abs(aaa$ses)),
                    ses_abond=sum(aaa$ses_abond),
                    ab_ses_abond=sum(abs(aaa$ses_abond)),
                    abond=sum(aaa$abond),
                    nb_p=length(which(sign(aaa$ses)==1)),
                    nb_n=length(which(sign(aaa$ses)==-1)),
                    tot=nrow(aaa))
    
    if((aaa$nb_p+aaa$nb_n)>0){
      result_intensity[i+1]<-aaa$ab_ses_abond/aaa$abond
      result_attractiveness[i+1]<-(aaa$nb_p-aaa$nb_n)/(aaa$nb_p+aaa$nb_n)
    }else{
      result_intensity[i+1]<-0
      result_attractiveness[i+1]<-0
    }
  }
  
  # compute the SES
  
  result<-data.frame(habit=hab,
                     zonebio=zb,
                     intensity=(result_intensity[1]-mean(result_intensity[2:(nb.iter+1)]))/sd(result_intensity[2:(nb.iter+1)]),
                     attractiveness=(result_attractiveness[1]-mean(result_attractiveness[2:(nb.iter+1)]))/sd(result_attractiveness[2:(nb.iter+1)]))
  
  
  return(result)
}

community_association2<-function(x, nb.iter){tryCatch(community_association(x, nb.iter),
                                                    error=function(e) data.frame(habit="NA", zonebio="NA", intensity=NA,attractiveness=NA))}

com_association<-ddply(droplevels(dataprp), .(code_point, annee), .fun=community_association2,.parallel =T, nb.iter=200)


#### Computing community structure ----

library(igraph)
library(ggraph)

a <- matrixinteraction_final

deg_comm<-function(x, threshold){
  links<-a[,c(3,4,6,1,2)]
  links<-na.omit(links)
  links$ses<-abs(links$ses)
  
  # select association from the right habitat and biogeographic region
  links<-subset(links, links$habit==names(which.max(table(x$habit))))
  links<-subset(links, links$zonebio==tableassopoint$zonebio[match(levels(droplevels(x$code_point)), tableassopoint$code_point)])
  colnames(links)[3] <- "weight"
  rownames(links) <- NULL
  
  # specify the species
  nodes<-x[,c("code_sp","code_point","annee")]
  nodes2<-droplevels(nodes)
  nodes2<-nodes2 %>% group_by(code_sp) %>% summarize(count=n())
  links2<-subset(links, links$spi %in% levels(droplevels(nodes2$code_sp)))
  links2<-droplevels(subset(links2, links2$spj %in% levels(droplevels(nodes2$code_sp))))
  
  # make the network
  net <- graph_from_data_frame(d=links2, vertices=nodes2, directed=F)
  
  # calculate clique structure
  nb_cliq_3_std_3<-length(cliques(net,min=3))/(2^length(nodes2$code_sp)-(1+(length(nodes2$code_sp)^2)/2+length(nodes2$code_sp)/2))
  return(c(nb_cliq_3_std_3))
}

com_structure<-ddply(dataprc,.(code_point,annee),.parallel=F, .progress = "text",deg_comm)
com_structure<-merge(com_structure, carre_centroid2[,-6], by="code_carre")

#### Computing functional dissimilarity ----

library(ade4)
library(RVAideMemoire)
library(cluster)
library(geometry)

# load the bird trait database
trait<-read.csv("life_history_bird2.csv",header = TRUE)
trait<-trait[,-c(1:3)]
table.na<-apply(trait,2,function(x){sum(is.na(x))})
trait2<-na.omit(trait[,-which(table.na>12)])
trait3<-trait2[,-c(1,3,6,9,12,15,19,20)]
trait3[,c(26:66)]<-apply(trait3[,c(26:66)], 2, function(x) replace(x,x>0,1))
trait3[,c(12,26:66)]<-lapply(trait3[,c(12,26:66)], factor)
trait3[,23]<-as.numeric(trait3[,23])
mat.trait<-trait3
# the following line require a update of the scientific names used in the STOC data -> column "scientific_names2"
mat.trait$code_sp<-species_name$code_sp[match(trait2$Species,species_name$scientific_name2)]
mat.trait<-subset(mat.trait, code_sp %in% levels(droplevels(dataprp2$code_sp)))

# calculate gower distance
mat.dist<-daisy(mat.trait[,1:66], metric="gower", type=list(asymm=c(26:66)))
attr(mat.dist, "Labels")<-droplevels(mat.trait$code_sp)

# selecte clustering algorithm
library(clue)
clust1<-hclust(mat.dist, method = "ward.D")
clust2<-hclust(mat.dist, method = "ward.D2")
clust3<-hclust(mat.dist, method = "single")
clust4<-hclust(mat.dist, method = "complete")
clust5<-hclust(mat.dist, method = "average")
clust6<-hclust(mat.dist, method = "mcquitty")


list.clust<-list(clust1,clust2,clust3,clust4,clust5,clust6)

dissim<-data.frame(ward.D=rep(NA,(2^length(list.clust)-1)),ward.D2=rep(NA,(2^length(list.clust)-1)),single=rep(NA,(2^length(list.clust)-1)),
                   complete=rep(NA,(2^length(list.clust)-1)),average=rep(NA,(2^length(list.clust)-1)),mcquitty=rep(NA,(2^length(list.clust)-1)),
                   dissimilarity=rep(NA,(2^length(list.clust)-1)))

nb<-0
for(i in 1:length(list.clust)){
  combinaison<-combn(length(list.clust),i)
  for(j in 1:ncol(combinaison)){
    nb<-nb+1
    print(nb)
    ensemble<-cl_ensemble(list=list.clust[combinaison[,j]])
    if(i==1){consensus<-list.clust[combinaison[,j]]
    }else{consensus<-cl_consensus(ensemble, method = "euclidean")}
    dissim[nb,c(combinaison[,j])]<-combinaison[,j]
    dissim[nb,9]<-cl_dissimilarity(mat.dist, consensus, method = "spectral") # use 2-norm method for algorithm selection
  }
}

dissim[which.min(dissim$dissimilarity),]
dissim[order(dissim$dissimilarity),]
consensus<-cl_consensus(cl_ensemble(list=list.clust[c(5)]))
plot(consensus)

# calculate the maximal threshold  for the distosion between the cophenetic matrix and the original distance matrix
sigma <- var(mat.dist) + var(cophenetic(clust5))
seuil <- 2*sqrt(nrow(as.matrix(mat.dist))*sigma)

# extract cophenetic distances
coph.dist<-cophenetic(consensus)
coph.dist2 <- data.frame(t(combn(attr(coph.dist, "Labels"),2)), as.numeric(coph.dist))
names(coph.dist2)<-c("sp1","sp2","distance")

# compare fonctional dissimilarity and pairwise association 

a<- matrixinteraction_final
a<-data.frame(a, pair=paste(a$spi,a$spj,sep="_"))
b<-coph.dist2
b<-data.frame(b, pair=paste(b$sp1,b$sp2,sep="_"))
b<-data.frame(b, pair2=paste(b$sp2,b$sp1,sep="_"))
a$dist<-b$distance[match(a$pair, b$pair)]
a$dist2<-b$distance[match(a$pair, b$pair2)]
a$dist[which(is.na(a$dist))]<-0
a$dist2[which(is.na(a$dist2))]<-0
a$dist3<-a$dist2+a$dist

aa<-na.omit(subset(a, !(dist3==0)))
lm_asso_fd<-lm(log(abs(ses)+1)~dist3, data=aa)
summary(lm_asso_fd)

aa<-aa %>% group_by(pair) %>% summarize(dist3=mean(dist3), ses=mean(ses))
lm_asso_fd<-lm(log(abs(ses)+1)~dist3, data=aa)
summary(lm_asso_fd)


#### Computing regional values ----

# Intensity and attractiveness
require(spgwr)

com_association$code_carre <- substr(com_association$code_point, 1, 6)
names(com_association)[5:6]<-c("intensity","attractiveness")

com_association1<- as.data.frame(com_association %>% group_by(code_carre, annee) %>% summarize(intensity=mean(intensity),
                                                                            attractiveness=mean(attractiveness)))
com_association2<-com_association %>% group_by(code_carre) %>% summarize(count=n())
com_association2$lon<-carre_centroid2$long[match(com_association2$code_carre, carre_centroid2$code_carre)]
com_association2$lat<-carre_centroid2$lat[match(com_association2$code_carre, carre_centroid2$code_carre)]
com_association2<-na.omit(com_association2)
testcarrecoord<-com_association2[,c("lon","lat")]
testcarrecoord<-SpatialPoints(testcarrecoord,proj4string = CRS("+proj=longlat +datum=WGS84"))
testcarrecoordlamb<-spTransform(testcarrecoord, CRS("+init=epsg:27572"))
matdist<-data.frame(spDists(testcarrecoordlamb))
b<-na.omit(com_association1)
b$annee<-as.numeric(as.character(b$annee))

asso_an_gwr<-function(i, bandwidth){
  data_dist<-data.frame(code_carre=com_association2$code_carre, dist=matdist[,i])
  data_dist$weight<-gwr.bisquare(data_dist$dist^2, bandwidth)
  data_dist<-droplevels(subset(data_dist, weight!=0))
  b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% levels(as.factor(data_dist$code_carre))))
  b_20$weight<-data_dist$weight[match(b_20$code_carre, data_dist$code_carre)]
  intensitymoy<-lm(intensity~annee, weights = weight, data=b_20)
  attractivenessmoy<-lm(attractiveness~annee, weights = weight, data=b_20)
  mean_intensity<-weighted.mean(b_20$intensity, b_20$weight)
  mean_attractiveness<-weighted.mean(b_20$attractiveness, b_20$weight)
  return(c(intensitymoy$coefficients[2],summary(intensitymoy)$coefficients[2,4],
           attractivenessmoy$coefficients[2],summary(attractivenessmoy)$coefficients[2,4],
           mean_intensity, mean_attractiveness))
}


bandwidths<-seq(from=10000, to=200000, length.out=20)

for(j in 1:length(bandwidths)){
  colnum<-ncol(com_association2)
  print(j)
  for(i in 1:length(matdist)){
    #print(i)
    aa<-asso_an_gwr(i, bandwidths[j])
    com_association2[i, colnum+1]<-aa[1]
    com_association2[i, colnum+2]<-aa[2]
    com_association2[i, colnum+3]<-aa[3]
    com_association2[i, colnum+4]<-aa[4]
    com_association2[i, colnum+5]<-aa[5]
    com_association2[i, colnum+6]<-aa[6]}
  names(com_association2)[(colnum+1):(colnum+6)]<-c(paste0("t_intensity", j),paste0("p_t_intensity", j),
                                                    paste0("t_attractiveness", j),paste0("p_t_attractiveness", j),
                                                    paste0("intensity", j),paste0("attractiveness", j))
}


com_association2<-as.data.frame(com_association2)

com_association2$t_intensity<-apply(com_association2[,5:124], 1,
                                      FUN=function(x){
                                        z<-x[seq(from=1,to=115,by=6)]
                                        z2<-x[seq(from=2,to=116,by=6)]
                                        return(mean(z[which(z2<=0.05)],na.rm=T))
                                      })
com_association2$t_attractiveness<-apply(com_association2[,5:124], 1,
                               FUN=function(x){
                                 z<-x[seq(from=3,to=117,by=6)]
                                 z2<-x[seq(from=4,to=118,by=6)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
com_association2$intensity<-apply(com_association2[,seq(from=9,to=123,by=6)], 1, FUN=function(x){mean(x,na.rm=T)})
com_association2$attractiveness<-apply(com_association2[,seq(from=10,to=124,by=6)], 1, FUN=function(x){mean(x,na.rm=T)})

com_association2<-merge(com_association2, carre_centroid2[,c(1,4,5)], by="code_carre")

# Clique structure

com_structure$code_carre <- substr(com_structure$code_point, 1, 6)
names(com_structure)[3]<-"clique"
com_structure1 <- as.data.frame(com_structure %>% group_by(code_carre, annee) %>% summarize(clique=mean(clique, na.rm=T)))

com_structure1$lat<-carre_centroid2$lat[match(com_structure1$code_carre, carre_centroid2$code_carre)]
com_structure1$lon<-carre_centroid2$long[match(com_structure1$code_carre, carre_centroid2$code_carre)]
com_structure1$lat2<-carre_centroid2$lat2[match(com_structure1$code_carre, carre_centroid2$code_carre)]
com_structure1$lon2<-carre_centroid2$lon2[match(com_structure1$code_carre, carre_centroid2$code_carre)]
com_structure2<-as.data.frame(com_structure1 %>% group_by(code_carre, lat, lon) %>% summarize(count=n()))
com_structure2<-na.omit(com_structure2)
carrecoord<-com_structure2[,c("lon","lat")]
carrecoord<-SpatialPoints(carrecoord,proj4string = CRS("+proj=longlat +datum=WGS84"))
carrecoordlamb<-spTransform(carrecoord, CRS("+init=epsg:27572"))
matdist<-data.frame(spDists(carrecoordlamb))

bandwidths<-seq(from=10000, to=200000, length.out=20)

clique_gwr2<-function(i, bandwidth){
  data_dist<-data.frame(code_carre=com_structure2$code_carre, dist=matdist[,i])
  data_dist$weight<-gwr.bisquare(data_dist$dist^2, bandwidth)
  data_dist<-droplevels(subset(data_dist, weight!=0))
  b_20<-droplevels(subset(com_structure1, com_structure1$code_carre %in% levels(as.factor(data_dist$code_carre))))
  b_20$weight<-data_dist$weight[match(b_20$code_carre, data_dist$code_carre)]
  t_clique<-lm(log(clique+0.001)~as.numeric(annee), data=na.omit(b_20), weights = weight)
  clique<-weighted.mean(b_20$clique, b_20$weight, na.rm=T)
  return(c(t_clique$coefficients[2],summary(t_clique)$coefficients[2,4],clique))
}

clique_gwr3<-function(i, bandwidth){tryCatch(clique_gwr2(i, bandwidth),
                                             error=function(e) c(NA,NA,NA))}

for(j in 1:length(bandwidths)){
  colnum<-ncol(com_structure2)
  print(j)
  for(i in 1:length(matdist)){
    aa<-clique_gwr3(i, bandwidths[j])
    com_structure2[i, colnum+1]<-aa[1]
    com_structure2[i, colnum+2]<-aa[2]
    com_structure2[i, colnum+3]<-aa[3]
  }
  names(com_structure2)[(colnum+1):(colnum+3)]<-c(paste0("t_clique", j),paste0("p_t_clique", j),paste0("clique", j))
}

com_structure2<-as.data.frame(com_structure2)

com_structure2$t_clique<-apply(com_structure2[,4:63], 1,
                                    FUN=function(x){
                                      z<-x[seq(from=1,to=58,by=3)]
                                      z2<-x[seq(from=2,to=59,by=3)]
                                      return(mean(z[which(z2<=0.05)],na.rm=T))
                                    })
com_structure2$clique<-apply(com_structure2[,seq(from=6,to=63,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})
com_structure2$t_clique_ns<-apply(com_structure2[,seq(from=4,to=63,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})

com_structure2<-merge(com_structure2, carre_centroid2[,c(1,4,5)], by="code_carre")


#### Calculating beta (and gamma) diversity ----


# Beta-diversity

bdiv<-dataprc[,c(3,8,11,15,18:22)]
bdiv$zonebio<-carre_centroid2$zonebio[match(bdiv$code_carre,carre_centroid2$code_carre)]
habitat_maj<-bdiv %>% group_by(code_carre, habit) %>% summarize(count=n()) %>% slice(which.max(count))
habitat_maj_an<-bdiv %>% group_by(code_carre, annee, habit) %>% summarize(count=n()) %>% slice(which.max(count))
bdiv_c<-dcast(bdiv, code_carre+zonebio~code_sp, fun.aggregate = sum,value.var="abond")
bdiv_c$habit<-habitat_maj$habit
bdiv_c_an<-dcast(bdiv, code_carre+annee+zonebio~code_sp, fun.aggregate = sum,value.var="abond")
bdiv_c_an$habit<-habitat_maj_an$habit

data_div<-bdiv_c[,c(1,2,111)]
data_div$lon<-carre_centroid2$long[match(data_div$code_carre, carre_centroid2$code_carre)]
data_div$lat<-carre_centroid2$lat[match(data_div$code_carre, carre_centroid2$code_carre)]
data_div$lon2<-carre_centroid2$lon2[match(data_div$code_carre, carre_centroid2$code_carre)]
data_div$lat2<-carre_centroid2$lat2[match(data_div$code_carre, carre_centroid2$code_carre)]

testcarrecoord<-data_div[,c("lon2","lat2")]
testcarrecoord<-SpatialPoints(testcarrecoord,proj4string = CRS("+init=epsg:27572"))
matdist<-data.frame(spDists(testcarrecoord))
b<-bdiv_c_an

require(entropart)

div_an_unif<-function(i, bandwidth){
  data_dist<-data.frame(code_carre=data_div$code_carre, dist=matdist[,i])
  data_dist<-droplevels(subset(data_dist, dist<bandwidth))
  b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% levels(as.factor(data_dist$code_carre))))
  b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% levels(as.factor(data_dist$code_carre))))
  if(nrow(b_c)<10){#choose the minimum number of square to compute the beta-diversity
    beta_div<-data.frame(rep(NA,10))
    pente<-data.frame(rep(NA,10))
    p.value<-data.frame(rep(NA,10))
  }
  if(nrow(b_c)>=10 & nrow(b_c)<20){
    beta_div<-data.frame(rep(NA,10))
    pente<-data.frame(rep(NA,10))
    p.value<-data.frame(rep(NA,10))
    sous_ech_carre<-c(as.character(data_div$code_carre[i]),
                      sample(levels(as.factor(data_dist$code_carre))[which(levels(as.factor(data_dist$code_carre))!=data_div$code_carre[i])],9))
    b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% sous_ech_carre))
    b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% sous_ech_carre))
    beta_div[1,1]<-BetaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_carre,zonebio,habit)))))$Total#calculate the beta-diversity using all years
    beta_div[,1]<-rep(beta_div[1,1], 10)
    bdiv_record<-data.frame()
    for(z in 1:length(levels(droplevels(as.factor(b_20$annee))))){
      if(nrow(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z]))<5){
        bdiv_record[z,1]<-NA
        bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
      }else{
        bdiv_record[z,1]<-BetaDiversity(MetaCommunity(t(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z])[,-which(names(b_20) %in% c("code_carre","annee","zonebio","habit"))])))$Total
        bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
      }
    }#calculate the beta-diversity for each year
    if(nrow(na.omit(bdiv_record))<2){
      pente[,1]<-rep(NA,10)
      p.value[,1]<-rep(NA,10)
    }else{
      beta_pente<-lm(bdiv_record[,1]~as.numeric(bdiv_record[,2]))#calculate beta-diversity trend
      pente[,1]<- rep(beta_pente$coefficients[2],10)
      p.value[,1]<-rep(summary(beta_pente)$coefficients[2,4],10)
    }
  }
  if(nrow(b_c)>=20){# if more than 20 square in the window, calculate several time the beta-diversity on 10 square samples
    beta_div<-data.frame()
    pente<-data.frame()
    p.value<-data.frame()
    for(k in 1:10){
      sous_ech_carre<-c(as.character(data_div$code_carre[i]),
                        sample(levels(as.factor(data_dist$code_carre))[which(levels(as.factor(data_dist$code_carre))!=data_div$code_carre[i])],9))
      b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% sous_ech_carre))
      b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% sous_ech_carre))
      beta_div[k,1]<-BetaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_carre,zonebio,habit)))))$Total
      bdiv_record<-data.frame()
      for(z in 1:length(levels(droplevels(as.factor(b_20$annee))))){
        if(nrow(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z]))<5){
          bdiv_record[z,1]<-NA
          bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
        }else{
          bdiv_record[z,1]<-BetaDiversity(MetaCommunity(t(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z])[,-which(names(b_20) %in% c("code_carre","annee","zonebio","habit"))])))$Total
          bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
        }
      }
      if(nrow(na.omit(bdiv_record))<2){
        pente[k,1]<-NA
        p.value[k,1]<-NA
      }else{
        beta_pente<-lm(bdiv_record[,1]~as.numeric(bdiv_record[,2]))
        pente[k,1]<-beta_pente$coefficients[2]
        p.value[k,1]<-summary(beta_pente)$coefficients[2,4]
      }
      
    }
  }
  return(data.frame(beta_div, pente, p.value))
}


bandwidths<-seq(from=10000, to=200000, length.out=20)
for(j in 1:length(bandwidths)){
  colnum<-ncol(data_div)
  print(j)
  data_div[,(colnum+1):(colnum+30)]<-NA
  for(i in 1:length(matdist)){
    print(i)
    aa<-div_an_unif(i, bandwidths[j])
    data_div[i, (colnum+1):(colnum+10)]<-aa[,1]
    data_div[i, (colnum+11):(colnum+20)]<-aa[,2]
    data_div[i, (colnum+21):(colnum+30)]<-aa[,3]}
  names(data_div)[c((colnum+1),(colnum+11),(colnum+21))]<-c(paste0("beta_div", j), paste0("beta_p", j), paste0("p_val", j))
}


data_div2<-data_div[,1:7]
for(i in 1:60){
  data_div2[,i+7]<-rowMeans(data_div[,(8+10*(i-1)):(7+10*i)], na.rm=T)
  names(data_div2)[i+7]<-names(data_div)[(8+10*(i-1))]
}
data_div2$beta_div<-apply(data_div2[,seq(from=8,to=65,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})
data_div2$beta_p<-apply(data_div2[,8:67], 1,
                               FUN=function(x){
                                 z<-x[seq(from=2,to=59,by=3)]
                                 z2<-x[seq(from=3,to=60,by=3)]
                                 return(mean(z[which(z2<=0.05)],na.rm=T))
                               })
data_div2$beta_p_ns<-apply(data_div2[,seq(from=9,to=66,by=3)], 1, FUN=function(x){mean(x,na.rm=T)})

# Gamma diversity

gdiv_an_unif<-function(i, bandwidth){
  data_dist<-data.frame(code_carre=data_div$code_carre, dist=matdist[,i])
  data_dist<-droplevels(subset(data_dist, dist<bandwidth))
  b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% levels(as.factor(data_dist$code_carre))))
  b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% levels(as.factor(data_dist$code_carre))))
  if(nrow(b_c)<10){
    beta_div<-data.frame(rep(NA,10))
    pente<-data.frame(rep(NA,10))
    p.value<-data.frame(rep(NA,10))
  }
  if(nrow(b_c)>=10 & nrow(b_c)<20){
    beta_div<-data.frame(rep(NA,10))
    pente<-data.frame(rep(NA,10))
    p.value<-data.frame(rep(NA,10))
    sous_ech_carre<-c(as.character(data_div$code_carre[i]),
                      sample(levels(as.factor(data_dist$code_carre))[which(levels(as.factor(data_dist$code_carre))!=data_div$code_carre[i])],9))
    b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% sous_ech_carre))
    b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% sous_ech_carre))
    beta_div[1,1]<-GammaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_carre,zonebio,habit)))))
    beta_div[,1]<-rep(beta_div[1,1], 10)
    bdiv_record<-data.frame()
    for(z in 1:length(levels(droplevels(as.factor(b_20$annee))))){
      if(nrow(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z]))<5){
        bdiv_record[z,1]<-NA
        bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
      }else{
        bdiv_record[z,1]<-GammaDiversity(MetaCommunity(t(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z])[,-which(names(b_20) %in% c("code_carre","annee","zonebio","habit"))])))
        bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
      }
    }
    if(nrow(na.omit(bdiv_record))<2){
      pente[,1]<-rep(NA,10)
      p.value[,1]<-rep(NA,10)
    }else{
      beta_pente<-lm(bdiv_record[,1]~as.numeric(bdiv_record[,2]))
      pente[,1]<- rep(beta_pente$coefficients[2],10)
      p.value[,1]<-rep(summary(beta_pente)$coefficients[2,4],10)
    }
  }
  if(nrow(b_c)>=20){
    beta_div<-data.frame()
    pente<-data.frame()
    p.value<-data.frame()
    for(k in 1:10){
      sous_ech_carre<-c(as.character(data_div$code_carre[i]),
                        sample(levels(as.factor(data_dist$code_carre))[which(levels(as.factor(data_dist$code_carre))!=data_div$code_carre[i])],9))
      b_20<-droplevels(subset(b, as.factor(b$code_carre) %in% sous_ech_carre))
      b_c<-droplevels(subset(droplevels(bdiv_c), code_carre %in% sous_ech_carre))
      beta_div[k,1]<-GammaDiversity(MetaCommunity(t(subset(b_c,select=-c(code_carre,zonebio,habit)))))
      bdiv_record<-data.frame()
      for(z in 1:length(levels(droplevels(as.factor(b_20$annee))))){
        if(nrow(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z]))<5){
          bdiv_record[z,1]<-NA
          bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
        }else{
          bdiv_record[z,1]<-GammaDiversity(MetaCommunity(t(subset(b_20, annee %in% levels(droplevels(as.factor(b_20$annee)))[z])[,-which(names(b_20) %in% c("code_carre","annee","zonebio","habit"))])))
          bdiv_record[z,2]<-levels(droplevels(as.factor(b_20$annee)))[z]
        }
      }
      if(nrow(na.omit(bdiv_record))<2){
        pente[k,1]<-NA
        p.value[k,1]<-NA
      }else{
        beta_pente<-lm(bdiv_record[,1]~as.numeric(bdiv_record[,2]))
        pente[k,1]<-beta_pente$coefficients[2]
        p.value[k,1]<-summary(beta_pente)$coefficients[2,4]
      }
      
    }
  }
  return(data.frame(beta_div, pente, p.value))
}


#### Testing relationships ----

a<-merge(com_association2, com_structure2[,-c(2,3,67,68)], by.x=c("code_carre"), by.y=c("code_carre"))
a<-merge(a, data_div2[,-c(2:7)],by.x=c("code_carre"), by.y=c("code_carre"))

rela<-data.frame(cor=c(summary(gam(a[,128]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,128]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20))) #attractiveness vs intensity
for(i in 1:20){
  rela[i+1,1]<-summary(gam(a[,6*i+4]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela[i+1,2]<-summary(gam(a[,6*i+4]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela2<-data.frame(cor=c(summary(gam(a[,127]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,127]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20))) #attractiveness vs intensity
for(i in 1:20){
  rela2[i+1,1]<-summary(gam(a[,6*i+3]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela2[i+1,2]<-summary(gam(a[,6*i+3]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela3<-data.frame(cor=c(summary(gam(a[,128]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,128]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20))) #attractiveness vs intensity
for(i in 1:20){
  rela3[i+1,1]<-summary(gam(a[,6*i+4]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela3[i+1,2]<-summary(gam(a[,6*i+4]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_sp<-as.data.frame(cbind(scale(rela[,1],center = F), rela[,2:3], scale(rela2[,1],center = F),rela2[,2], scale(rela3[,1],center = F),rela3[,2]))
names(rela_sp)<-c("attractiveness/intensity","pval1","scale","clique/intensity","pval2","clique/attractiveness","pval3")
rela_sp1<-melt(rela_sp, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_sp2<-melt(rela_sp, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_sp1$pvalue<- sign(rela_sp2$value-0.05)
ggplot(rela_sp1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Scale (km)", y = "Correlation")

rela4<-data.frame(cor=c(summary(gam(a[,126]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,126]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa<-a[which(a[,6*i]<=0.05 & a[,6*i+2]<=0.05),]
  rela4[i+1,1]<-summary(gam(aa[,6*i+1]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] #attractiveness trend vs intensity trend
  rela4[i+1,2]<-summary(gam(aa[,6*i+1]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela5<-data.frame(cor=c(summary(gam(a[,125]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,125]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa<-a[which(a[,6*i]<=0.05 & a[,3*i+129]<=0.05),]
  rela5[i+1,1]<-summary(gam(aa[,6*i-1]~aa[,3*i+128]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] #clique trend vs intensity trend
  rela5[i+1,2]<-summary(gam(aa[,6*i-1]~aa[,3*i+128]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela6<-data.frame(cor=c(summary(gam(a[,126]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,126]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 1:20){
  aa<-a[which(a[,6*i+2]<=0.05 & a[,3*i+129]<=0.05),]
  rela6[i+1,1]<-summary(gam(a[,6*i+1]~a[,3*i+128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]] #clique trend vs attractiveness trend
  rela6[i+1,2]<-summary(gam(a[,6*i+1]~a[,3*i+128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_trend<-as.data.frame(cbind(scale(rela4[,1],center = F), rela4[,2:3], scale(rela5[,1],center = F),rela5[,2], scale(rela6[,1],center = F),rela6[,2]))
names(rela_trend)<-c("attractiveness/intensity","pval1","scale","clique/intensity","pval2","clique/attractiveness","pval3")
rela_trend1<-melt(rela_trend, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_trend2<-melt(rela_trend, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_trend1$pvalue<- sign(rela_trend2$value-0.05)
ggplot(rela_trend1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Scale (km)", y = "Correlation")


rela_bv1<-data.frame(cor=c(summary(gam(a[,254]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                 pval=c(summary(gam(a[,254]~a[,127]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                 echelle=c(0,seq(from=10,to=200,length.out=20))) #beta-diversity vs intensity
for(i in 2:20){
  rela_bv1[i+1,1]<-summary(gam(a[,3*i+191]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv1[i+1,2]<-summary(gam(a[,3*i+191]~a[,6*i+3]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bv2<-data.frame(cor=c(summary(gam(a[,254]~a[,128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                     pval=c(summary(gam(a[,254]~a[,128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                     echelle=c(0,seq(from=10,to=200,length.out=20))) #beta-diversity vs attractiveness
for(i in 2:20){
  rela_bv2[i+1,1]<-summary(gam(a[,3*i+191]~a[,6*i+4]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv2[i+1,2]<-summary(gam(a[,3*i+191]~a[,6*i+4]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bv3<-data.frame(cor=c(summary(gam(a[,254]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                     pval=c(summary(gam(a[,254]~a[,192]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                     echelle=c(0,seq(from=10,to=200,length.out=20))) #beta-diversity vs clique
for(i in 2:20){
  rela_bv3[i+1,1]<-summary(gam(a[,3*i+191]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]]
  rela_bv3[i+1,2]<-summary(gam(a[,3*i+191]~a[,3*i+130]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_betadiv<-as.data.frame(cbind(scale(rela_bv1[,1],center = F), rela_bv1[,2:3], scale(rela_bv2[,1],center = F),rela_bv2[,2], scale(rela_bv3[,1],center = F),rela_bv3[,2]))
names(rela_betadiv)<-c("intensity","pval1","scale","attractiveness","pval2","clique","pval3")
rela_betadiv1<-melt(rela_betadiv, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_betadiv2<-melt(rela_betadiv, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_betadiv1$pvalue<- sign(rela_betadiv2$value-0.05)
ggplot(rela_betadiv1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Scale (km)", y = "Correlation")


rela_bvt1<-data.frame(cor=c(summary(gam(a[,255]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,255]~a[,125]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa<-a[a[,6*i]<=0.05,]
  rela_bvt1[i+1,1]<-summary(gam(aa[,3*i+192]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] #beta-diversity trend vs intensity trend
  rela_bvt1[i+1,2]<-summary(gam(aa[,3*i+192]~aa[,6*i-1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bvt2<-data.frame(cor=c(summary(gam(a[,255]~a[,126]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,255]~a[,126]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 3:20){
  aa<-a[a[,6*i+2]<=0.05,]
  rela_bvt2[i+1,1]<-summary(gam(aa[,3*i+192]~aa[,6*i+1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.coef[[2]] #beta-diversity trend vs attractiveness trend
  rela_bvt2[i+1,2]<-summary(gam(aa[,3*i+192]~aa[,6*i+1]+te(aa$lon2,aa$lat2,bs="tp",k=3)))$p.pv[[2]]
}

rela_bvt3<-data.frame(cor=c(summary(gam(a[,255]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]],rep(NA,20)),
                  pval=c(summary(gam(a[,255]~a[,191]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]],rep(NA,20)),
                  echelle=c(0,seq(from=10,to=200,length.out=20)))
for(i in 2:20){
  aa<-a[a[,3*i+129]<=0.05,]
  rela_bvt3[i+1,1]<-summary(gam(a[,3*i+192]~a[,3*i+128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.coef[[2]] #beta-diversity trend vs clique trend
  rela_bvt3[i+1,2]<-summary(gam(a[,3*i+192]~a[,3*i+128]+te(a$lon2,a$lat2,bs="tp",k=3)))$p.pv[[2]]
}


rela_betadiv_t<-as.data.frame(cbind(scale(rela_bvt1[,1],center = F), rela_bvt1[,2:3], scale(rela_bvt2[,1],center = F),rela_bvt2[,2], scale(rela_bvt3[,1],center = F),rela_bvt3[,2]))
names(rela_betadiv_t)<-c("intensity","pval1","scale","attractiveness","pval2","clique","pval3")
rela_betadiv_t1<-melt(rela_betadiv_t, id.vars=c("scale"), measure.vars=c(1,4,6))
rela_betadiv_t2<-melt(rela_betadiv_t, id.vars=c("scale"), measure.vars=c(2,5,7))
rela_betadiv_t1$pvalue<- sign(rela_betadiv_t2$value-0.05)
ggplot(rela_betadiv_t1, aes(x = scale, y = value, group=variable, shape=variable))+
  geom_point(size=3, aes(col=as.character(pvalue)))+scale_color_grey()+
  theme_light(base_size = 20)+labs(x ="Scale (km)", y = "Correlation")


#### Mapping indices ----

prep_krig<-function(sub_data,     # data to interpolate
                    accuracy=200){# accuracy of the output
  names(sub_data)<-c("lon","lat","variable")
  sub_data<-na.omit(sub_data)
  coord_vars<-c("lat","lon")
  data_vars<-setdiff(colnames(sub_data), coord_vars)
  sp_points<-SpatialPoints(sub_data[,coord_vars])
  sp_df<-SpatialPointsDataFrame(sp_points, sub_data[,data_vars,drop=FALSE])
  pixels_per_side<-accuracy
  bottom.left<-apply(sp_points@coords,2,min)
  top.right<-apply(sp_points@coords,2,max)
  margin<-abs((top.right-bottom.left))/10
  bottom.left<-bottom.left-margin
  top.right<-top.right+margin
  pixel.size<-abs(top.right-bottom.left)/pixels_per_side
  g<-GridTopology(cellcentre.offset=bottom.left, cellsize=pixel.size, cells.dim=c(pixels_per_side,pixels_per_side))
  map_base_data<-map_data("france")
  map_base_data<-map_base_data[!map_base_data$region %in% c("Corse du Sud","Haute-Corse"), ]
  foo = function(x) {
    group = unique(x$group)
    Polygons(list(Polygon(x[,c("lat","long")])),ID=group)
  }
  state_pg<-SpatialPolygons(dlply(map_base_data, .(group), foo))
  grid_points<-SpatialPoints(g)
  in_points<-!is.na(over(grid_points,state_pg))
  fit_points<-SpatialPoints(as.data.frame(grid_points)[in_points,])
  krig<-autoKrige(variable~1, sp_df, new_data=fit_points, block = c(0.5,0.5))
  interp_data<-as.data.frame(krig$krige_output)
  colnames(interp_data)<-c("lat","lon","pred","var","stdev")
  return(interp_data)
}

sub_data<-com_association2[, c("lon","lat","intensity")]
interp_data1<-prep_krig(sub_data)

sub_data<-com_association2[, c("lon","lat","attractiveness")]
interp_data2<-prep_krig(sub_data)

sub_data<-na.omit(com_association2[, c("lon","lat","t_intensity")])
sub_data$t_intensity[abs(sub_data$t_intensity)>0.3]<-0.3*sign(sub_data$t_intensity[abs(sub_data$t_intensity)>0.3])
interp_data3<-prep_krig(sub_data)

sub_data<-na.omit(com_association2[, c("lon","lat","t_attractiveness")])
sub_data$t_attractiveness[abs(sub_data$t_attractiveness)>0.2]<-0.2*sign(sub_data$t_attractiveness[abs(sub_data$t_attractiveness)>0.2])
interp_data4<-prep_krig(sub_data)

sub_data<-com_structure2[, c("lon","lat","clique")]
interp_data5<-prep_krig(sub_data)

sub_data<-na.omit(com_structure2[, c("lon","lat","t_clique")])
sub_data$t_clique[abs(sub_data$t_clique)>0.2]<-0.2*sign(sub_data$t_clique[abs(sub_data$t_clique)>0.2])
interp_data6<-prep_krig(sub_data)

sub_data<-data_div2[, c("lon","lat","beta_div")]
interp_data7<-prep_krig(sub_data)

sub_data<-data_div2[, c("lon","lat","beta_p")]
interp_data8<-prep_krig(sub_data)
interp_data8b<-prep_krig(data_div2[, c("lon","lat","beta_p_ns")])

border_fr<-readOGR("ne_10m_admin_0_map_units.shp") # https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-details/
border_fr<-border_fr[which(border_fr$NAME=="France"),]
border_fr2<-crop(border_fr, extent(-5.14, 8.5, 42, 51.1))
border_fr3<- spTransform(border_fr2, CRS("+init=epsg:27572"))

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


corr_map <- function(map1, map2){
  couche_1<-rasterFromXYZ(map1[,c(2,1,3)])
  couche_2<-rasterFromXYZ(map2[,c(2,1,3)])
  couche_3<- stack(couche_1, couche_2)
  names(couche_3) <- c("1", "2")
  couche_corr<-raster(couche_3, 1)
  values(couche_corr) <- 1:ncell(couche_1)
  
  for(i in 1:20){
    print(i)
    weight.w<-matrix(1,2*i+1,2*i+1)
    focal_cor <- focal(
      x = couche_corr,
      pad=T,
      w=weight.w,
      fun = function(x, y = couche_3){
        cor(values(y)[x,1], values(y)[x,2],
            use = "na.or.complete")
      }
    )
    
    focal_cor <- focal_cor+couche_1-couche_1
    if(i<2){focal_cor_dat<-focal_cor}
    else{focal_cor_dat<-stack(focal_cor_dat,focal_cor)}
  }
  
  focal_1<-focal_cor_dat
  
  focal_cor_dat<-na.omit(as.data.frame(mean(focal_1) ,xy=T))
  out_map<-ggplot(focal_cor_dat) +
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

# example:
map_corr<-corr_map(interp_data1,interp_data2)
map_corr


#### National trend of the indices ----


model_intensity<-gamm(intensity~annee+te(lon2,lat2,bs="tp",k=3), data=com_association, random = list(code_carre= ~ 1))

model_attractiveness<-gamm(attractiveness~annee+te(lon2,lat2,bs="tp",k=3), data=com_association, random = list(code_carre= ~ 1))

model_clique<-gamm(log(clique+0.001)~annee+te(lon2,lat2,bs="tp",k=3), data=com_structure, random = list(code_carre= ~ 1), niterPQL = 50)

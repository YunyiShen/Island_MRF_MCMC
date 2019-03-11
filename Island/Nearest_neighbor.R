library(tripack)
require(RANN)
posi = read.csv(file.choose())

certain_island = posi[posi$Location=="BA",6:7]

distfull = dist(certain_island)

onlyisland = posi[posi$Location!="MD",]


require(qgraph)

distMfull = as.matrix(dist(onlyisland[,6:7]))
link_inner = (as.matrix(dist(onlyisland[,6:7]))<=1100)*1.0
diag(link_inner)=0

qgraph(link_inner, layout='spring', vsize=3)

island = unique(onlyisland$Location)


#distM_outer = distM_outer

#qgraph((distM_outer<=2000)+distM_inner, layout='spring', vsize=3)
#require(gtools)
library(tripack)

mos = voronoi.mosaic(onlyisland$X,onlyisland$Y, duplicate="remove")
plot(mos)

D_tri = mos$tri

neigh = neighbours(D_tri)

# use neighbours, now further simplify this, only related with the nearest exteral island

link_outer = 0*distMfull
for(i in 1:155){
  neighbors = neigh[[i]]
  island_on = onlyisland$Location[i]
  neighbor_ex_island = neighbors[onlyisland$Location[neighbors]!=island_on]
  if(length(neighbor_ex_island)==0){next}
  if(sum(link_inner[i])==4){next}
  distances = distMfull[i,neighbor_ex_island]
  
  nearst = neighbor_ex_island[which(distances==min(distances))]
  
  link_outer[i,nearst]=1
  link_outer[nearst,i]=1
  
  if(sum(link_inner[i,])==0){
    neighbor_in_island = neighbors[onlyisland$Location[neighbors]==island_on]
    distances = distMfull[i,neighbor_in_island]
    nearst = neighbor_in_island[which(distances==min(distances))]
    link_inner[i,nearst]=1
    link_inner[nearst,i]=1
  }
}

link_full = link_inner+exp(-0.0003*distMfull)*link_outer

qgraph(link_full, layout='spring', vsize=4,labels=onlyisland$Join1)

link_full = link_inner+link_outer
diag(link_full)=0

edges = which(link_full==1,T)

edges = edges[edges[,1]<edges[,2],]


all_edges = data.frame(edges,onlyisland[edges[,1],],onlyisland[edges[,2],])

lineWKT = apply(all_edges,1, function(X1){ paste0("LINESTRING(",paste(X1[c("X","Y")],collapse = " "),",",paste0(X1[c("X.1","Y.1")],collapse = " "),")")})

all_edges = data.frame(all_edges,lineWKT)

write.csv(all_edges,"all_link_draft.csv",row.names = F)
#https://gis.stackexchange.com/questions/225102/calculate-distance-between-points-and-nearest-polygon-in-r
# see that to takecare of mainland


install.packages("rgeos")

library(rgeos)
require(maptools)

spts = SpatialPoints(onlyisland[,6:7])

shp_Poly <- rgdal::readOGR("./GIS_yunyi/mainland.shp")

dist_to_mainland = (gDistance( shp_Poly,spts,byid=T))     

link_mainland = matrix(0,155,1)

for(i in 1:155){
  neighbors = neigh[[i]]
  island_on = onlyisland$Location[i]
  neighbor_ex_island = neighbors[onlyisland$Location[neighbors]!=island_on]
  if(length(neighbor_ex_island)==0){next}
  if(sum(link_inner[i])==4){next}
  distances = distMfull[i,neighbor_ex_island]
  
  nearst = min(distances)
  
  if(nearst>dist_to_mainland[i]) link_mainland[i]=1

  
}
write.csv(link_mainland,"link_mainland.csv",row.names = F)

# This is for graph
######
Graph_plain = read.csv("all_link_draft.csv")
posi = read.csv("CT_posi_only_island.csv")

distM_full = dist(posi[,6:7])
distM_full = as.matrix(distM_full)

row.names(distM_full) = as.character( posi$Join1)
colnames(distM_full) = row.names(distM_full)

graph_plain_inner = Graph_plain[ as.character( Graph_plain$Location)==as.character(Graph_plain$Location.1),1:2]
graph_plain_outer = Graph_plain[ as.character( Graph_plain$Location)!=as.character(Graph_plain$Location.1),1:2]

link_inner = 0*distM_full

link_inner[(graph_plain_inner$col-1)*155 + graph_plain_inner$row]=1
link_inner[(graph_plain_inner$row-1)*155 + graph_plain_inner$col]=1

link_outer = 0*distM_full


link_outer[(graph_plain_outer$col-1)*155 + graph_plain_outer$row]=1
link_outer[(graph_plain_outer$row-1)*155 + graph_plain_outer$col]=1

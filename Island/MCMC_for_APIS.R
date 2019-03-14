getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
  #eta = eta[1:nspp]
  nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
  A = list() # intralayer graphs are passed using lists
  if(int_range=="arth"){
    for(i in 1:nspp){
      A[[i]] = eta[i]*as.matrix(1/((distM)^(2+d[i])))
    }
  }
  else{
    if(int_range=="exp"){
      for(i in 1:nspp){
        A[[i]] = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
        diag(A[[i]])=0
      }
    }
    else{
      if(int_range=="nn"){
        for(i in 1:nspp){
          A[[i]] = eta[i]*as.matrix((link_map))
        }
      }
      else{
        #print("int_range must be exp or arth, will assume exp")
        for(i in 1:nspp){
          A[[i]] = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
          diag(A[[i]])=0
        }
      }
    }
  }
  return(A)
}

getfullGraph = function(A_ex,A_in,spp_mat){
  nspp = nrow(spp_mat)
  nsite = nrow(A_ex[[1]])
  A = matrix(0,nspp*nsite,nspp*nsite)
  for(i in 2:nspp-1){
    A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]] # diagonal part
    for(j in (i+1):nspp){
      
      diag(A[1:nsite + (i-1)*nsite,1:nsite + (j-1)*nsite])=spp_mat[i,j]
      diag(A[1:nsite + (j-1)*nsite,1:nsite + (i-1)*nsite])=spp_mat[j,i]
      
    }
  }
  i=nspp
  A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]]
  return(A)
}

SufStat_d_ex = function(Z,theta,link_outer,distM_full,link_mainland,distM_mainland,k){
  A = -exp(theta$d_ex[k] ) * theta$eta_ex[k] * distM_full * as.matrix(exp(-exp(theta$d_ex[k])*distM_full)) * (link_outer)
  B = -exp(theta$d_ex[k] ) * theta$eta_ex[k] * distM_mainland * as.matrix(exp(-exp(theta$d_ex[k])*distM_mainland)) * (link_mainland)
  
  return(.5*t(Z)%*%A%*%Z+sum(B*Z))
  
}

island = read.csv("CT_posi_only_island.csv")

link_inner = as.matrix( read.csv("link_inner.csv",row.names = 1))
#link_outer = as.matrix( read.csv("link_outer.csv",row.names = 1))
link_outer = as.matrix( read.csv("link_outer_full.csv",row.names = 1))
link_mainland = as.matrix( read.csv("link_mainland.csv"))
#link_outer = 0 * link_outer # this makes it a mainland-island system
#link_mainland = matrix(1,155,1)

distM_full = as.matrix( read.csv("distM_full.csv",row.names = 1))
distM_mainland = as.matrix( read.csv("dist_to_mainland.csv",row.names = 1))

intcd = min(min((distM_mainland*link_mainland)[(distM_mainland*link_mainland)>0]),
            min((link_outer*distM_full)[(link_outer*distM_full)>0]))
normd = max(max(distM_mainland*link_mainland),max(link_outer*distM_full))-intcd
  

distM_full = (distM_full-intcd)/normd # normalizing the distance
distM_mainland = (distM_mainland-intcd)/normd

spp_mat = matrix(c(0,1,1,0),2,2)

envX = matrix(1,155,1)
theta = list(beta = c(0,0),
             eta_in = c(.15,.15),
             eta_ex = c(.4,.4),
             d_ex = c(.2,.2),
             spp_mat = -0.15 * spp_mat)

A_in = getintralayerGraph(distM_full,link_inner,theta$eta_in,d,int_range = "nn",theta$spp_mat)
A_ex = getintralayerGraph(distM_full , link_outer,theta$eta_ex,theta$d_ex,int_range = "exp",theta$spp_mat)
G = getfullGraph(A_ex,A_in,theta$spp_mat)


ncov = ncol(envX)
thr = rbind(envX%*%theta$beta[1:ncov]+
              theta$eta_ex[1]*exp(-exp(theta$d_ex[1])*distM_mainland)*link_mainland,
            envX%*%theta$beta[1:ncov + ncov]+
              theta$eta_ex[2]*exp(-exp(theta$d_ex[2])*distM_mainland)*link_mainland
            )

require(IsingSampler)

set.seed(42)
Ising_sample = IsingSampler(n=500,G,thr,responses = c(-1,1),method="CFTP")
Ising_sample_mean = colMeans(((Ising_sample)+1)/2)
uniqueisland = as.character( unique(island$Location))

Z1_mean = apply(as.matrix(uniqueisland),1,function(uisland,Ising_mean,island){
  return(mean(Ising_mean[as.character(island$Location)== as.character(uisland)]))
}
                ,Ising_sample_mean[1:155]
                ,island)

Z2_mean = apply(as.matrix(uniqueisland),1,function(uisland,Ising_mean,island){
  return(mean(Ising_mean[as.character(island$Location)== as.character(uisland)]))
}
                ,Ising_sample_mean[1:155 + 155]
                ,island)

island_size = apply(as.matrix(uniqueisland),1,function(uisland,island){
  return(sum(as.character(island$Location)== as.character(uisland)))
}
,island)

island_dist = apply(as.matrix(uniqueisland),1,function(uisland,island){
  return(mean(distM_mainland[as.character(island$Location)== as.character(uisland)]))
}
,island)

plot(island_dist,(Z1_mean+Z2_mean))
plot(island_size,(Z1_mean+Z2_mean))


require(ggplot2)

num_sample_viewing = 1
tempdata = data.frame(island[,6:7],
                      Z_1 = Ising_sample_mean[1:155],
                      Z_2 = Ising_sample_mean[156:310])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1+Z_2))+
  geom_point()

I_beta_1 = (rowSums(Ising_sample[,1:155]))
I_beta_2 = (rowSums(Ising_sample[,1:155+155]))

I_eta_in_1 = (apply(Ising_sample[,1:155],1,function(Z,theta,A_in){0.5*(t(Z)%*%A_in[[1]]%*%(Z))/theta$eta_in[1]},theta,A_in))
I_eta_in_2 = (apply(Ising_sample[,1:155+155],1,function(Z,theta,A_in){0.5*(t(Z)%*%A_in[[2]]%*%(Z))/theta$eta_in[2]},theta,A_in))

I_eta_ex_1 = (apply(Ising_sample[,1:155],1,
                    function(Z,theta,A_ex,link_mainland,distM_mainland){
                      0.5*(t(Z)%*%A_ex[[1]]%*%(Z))/theta$eta_ex[1] + 
                        sum(exp(-exp(theta$d_ex[1])*distM_mainland)*link_mainland*Z)
                      },theta,A_ex,link_mainland,distM_mainland))

I_eta_ex_2 = (apply(Ising_sample[,1:155+155],1,
                    function(Z,theta,A_ex,link_mainland,distM_mainland){
                      0.5*(t(Z)%*%A_ex[[2]]%*%(Z))/theta$eta_ex[2] + 
                        sum(exp(-exp(theta$d_ex[2])*distM_mainland)*link_mainland*Z)
                    },theta,A_ex,link_mainland,distM_mainland))


I_d_ex_1 = apply(Ising_sample[,1:155],1,
                 SufStat_d_ex,
                 theta,link_outer,distM_full,link_mainland,distM_mainland,1)
I_d_ex_2 = apply(Ising_sample[,1:155+155],1,
                 SufStat_d_ex,
                 theta,link_outer,distM_full,link_mainland,distM_mainland,2)

I_eta = apply(Ising_sample,1,function(Z){sum(Z[1:155]*Z[1:155+155])})

FI_simu = data.frame(I_beta_1,I_beta_2,
                     I_eta_in_1,I_eta_in_2,
                     I_eta_ex_1,I_eta_ex_2,
                     I_d_ex_1,I_d_ex_2,
                     I_eta)

FI = cov(FI_simu)


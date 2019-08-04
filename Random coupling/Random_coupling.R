require(IsingSampler)
n_nodes = 100
sd_graph = 1
set.seed(42)
graph  = matrix(rnorm(n_nodes^2,-.3,sd_graph),n_nodes,n_nodes)
graph = 0.5 * (graph + t(graph))
diag(graph)=0
thr = matrix(0,n_nodes,1)

Zs = IsingSampler(n = 5000,graph = graph, thresholds = thr, beta = 1, responses = c(-1,1))

avg = rowMeans(Zs)
avg_spp = colMeans(Zs)
hist(avg)
hist(avg_spp)

library(trendsceek)
library(dplyr)

#Dataset 1: scRNA-seq data from study
data('scialdone')
counts = scialdone[['counts']]

#Dataset 2: scRNA-seq data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59739
"df <- read.delim('data.txt')
df[rowSums(df[, -1] > 0) != 0, ]
#format data 
temp <- df[-1:-4, ]
counts <- temp[-c(2, 865)]
temp1 <- counts$Sample.ID
counts[1] <- NULL
counts[is.na(counts)] <- 1
counts[counts==0] <- 1
counts %>% mutate_if(is.character,as.numeric)
counts <- data.frame(lapply(counts,as.numeric))
counts <- floor(counts)
row.names(counts) <- temp1"

#Tutorial by authors
min.ncells.expr = 3 #filter out genes w/ low expression
min.expr = 5
counts_filt = genefilter_exprmat(counts, min.expr, min.ncells.expr)

quantile.cutoff = 0.90 #drop lowest 10%
method = 'glm' #calculate gene variability
vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff, method = method)

n.topvar = 500 #select top 500
topvar.genes = rownames(vargenes_stats[['real.stats']])[1:n.topvar]

#plot cv2 vs avg read count
plot.ercc.points = FALSE
plot_cv2vsmean(vargenes_stats, topvar.genes, plot.ercc.points = plot.ercc.points)

#normalize counts
min.count = 1
counts_norm = deseq_norm(counts, min.count)
counts_sub = counts_norm[topvar.genes, ]

#tSNE
tsne.k = 2 
init.dims = 100 
perp.frac = 0.2 
max.iter = 300 #should converge
epoch = 50 
tsne_res = trend.tsne(counts_sub, tsne.k, init.dims, perp.frac, max.iter, epoch)

#create point pattern using tSNE positions and expression levels
pp = pos2pp(tsne_res)
log.fcn = log10
pp = set_marks(pp, counts_sub, log.fcn = log.fcn)

n.top2plot = 10 #select top 10 to run trendsceek on
topvar.genes = rownames(vargenes_stats[['real.stats']])[1:n.top2plot]
pp2plot = pp_select(pp, topvar.genes)
plot_pp_scatter(pp2plot, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)

#trendsceek
nrand = 50 #permutations, should probably set higher
ncores = 1
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

#select significant genes
alpha = 0.05
sig_list = extract_sig_genes(trendstat_list, alpha)
lapply(sig_list, nrow)

#plot significant mark correlations
sig_genes = sig_list[['markcorr']][, 'gene']
pp_sig = pp_select(pp, sig_genes)
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)
plot_pp_density(pp_sig, log_marks = FALSE)

#plot significant mark variograms
sig_genes = sig_list[['markvario']][, 'gene']
pp_sig = pp_select(pp, sig_genes)
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)
plot_pp_density(pp_sig, log_marks = FALSE)


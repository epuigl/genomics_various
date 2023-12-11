#------ HTODemux_low_high
# Context: usage of HTOdemux with large hashtag imbalance, assuming the positive signal 
# cutoffs should be different since some HTOs (also works with CMO data) have significantly
# more coverage than others. After visual inspection (e.g. signal A vs B scatterplots) and
# ONLY IF we are confident that the less covered hashtags must be kept for analyses, here
# is a modified HTOdemux (Seurat base) function so one can specify a positive.quantile_low 
#(anything below is negative) or positive.quantile_high, anything else is a singlet. 
# ! Use with caution and do not deviate too much from the default 0.99 threshold !
#------

library(cluster)
library(fitdistrplus)

MaxN <- function(x, N = 2){
len <- length(x)
if (N > len) {
warning('N greater than length(x). Setting N=length(x)')
N <- length(x)
}
sort(x, partial = len - N + 1)[len - N + 1]
}

HTODemux_low_high <- function (object, assay = "HTO", positive.quantile_low = 0.99, 
							   positive.quantile_high = 0.99, init = NULL, nstarts = 100,
							   kfunc = "clara", nsamples = 100, seed = 42, verbose = TRUE)
{
if (!is.null(x = seed)) {
set.seed(seed = seed)
}
assay <- assay %||% DefaultAssay(object = object)
data <- GetAssayData(object = object, assay = assay)
counts <- GetAssayData(object = object, assay = assay, slot = "counts")[,
colnames(x = object)]
counts <- as.matrix(x = counts)
ncenters <- init %||% (nrow(x = data) + 1)
switch(EXPR = kfunc, kmeans = {
init.clusters <- kmeans(x = t(x = GetAssayData(object = object,
assay = assay)), centers = ncenters, nstart = nstarts)
Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
}, clara = {
init.clusters <- clara(x = t(x = GetAssayData(object = object,
assay = assay)), k = ncenters, samples = nsamples)
Idents(object = object, cells = names(x = init.clusters$clustering),
drop = TRUE) <- init.clusters$clustering
}, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
average.expression <- AverageExpression(object = object,
assays = assay, verbose = FALSE)[[assay]]
if (sum(average.expression == 0) > 0) {
stop("Cells with zero counts exist as a cluster.")
}
discrete_weak <- GetAssayData(object = object, assay = assay)
discrete_weak[discrete_weak > 0] <- 0
discrete_strong <- discrete_weak

###start modified section
for (iter in rownames(x = data)) {
values <- counts[iter, colnames(object)]
values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter,
])]])]
fit <- suppressWarnings(expr = fitdist(data = values.use,
distr = "nbinom"))
cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile_low)$quantiles[1])
discrete_weak[iter, names(x = which(x = values > cutoff))] <- 1
if (verbose) {
message(paste0("Cutoff for ", iter, " : ", cutoff,
" reads"))
}
}

for (iter in rownames(x = data)) {
values <- counts[iter, colnames(object)]
values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter,
])]])]
fit <- suppressWarnings(expr = fitdist(data = values.use,
distr = "nbinom"))
cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile_high)$quantiles[1])
discrete_strong[iter, names(x = which(x = values > cutoff))] <- 1
if (verbose) {
message(paste0("Cutoff for ", iter, " : ", cutoff,
" reads"))
}
}

npositive_weak <- colSums(x = discrete_weak)
npositive_weak[npositive_weak >= 1] <- "Singlet"
npositive_weak[npositive_weak == 0] <- "Negative"

npositive_strong <- colSums(x = discrete_strong)
npositive_strong[npositive_strong >= 2] <- "Doublet"
npositive_strong[npositive_strong == 1 | npositive_strong == 0] <- "Singlet"

classification.global <- data.frame(npositive_weak, npositive_strong)
classification.global$npositive_strong[classification.global$npositive_weak == "Negative"] <- "Negative"
classification.global <- classification.global$npositive_strong
###End modififed section

donor.id = rownames(x = data)
hash.max <- apply(X = data, MARGIN = 2, FUN = max)
hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data),
FUN = function(x) {
return(which(x = data[, x] == hash.max[x])[1])
})])
hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data),
FUN = function(x) {
return(which(x = data[, x] == hash.second[x])[1])
})])
hash.margin <- hash.max - hash.second
doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])),
collapse = ""))
})
classification <- classification.global
classification[classification.global == "Negative"] <- "Negative"
classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global ==
"Singlet")]
classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global ==
"Doublet")]
classification.metadata <- data.frame(hash.maxID, hash.secondID,
hash.margin, classification, classification.global)
colnames(x = classification.metadata) <- paste(assay, c("maxID",
"secondID", "margin", "classification", "classification.global"),
sep = "")
object <- AddMetaData(object = object, metadata = classification.metadata)
Idents(object) <- paste0(assay, "_classification")
doublets <- rownames(x = object[[]])[which(object[[paste0(assay,
"_classification.global")]] == "Doublet")]
Idents(object = object, cells = doublets) <- "Doublet"
object$hash.ID <- Idents(object = object)
return(object)
}
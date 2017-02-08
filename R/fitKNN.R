#' fSVA and k-Nearest Neighbor classifier
#'
#' Couples together fSVA with Nearest Shrunken Centroid classifier.
#'
#' @param x matrix: Normalized and log-transformed expression values (samples in rows, features in columns).
#' @param y factor: Class labels.
#' @param n.sv integer: Number of surrogate variables to estimate.
#' @param k integer: Number of nearest neighbors casting votes.
#' @param ... additional arguments passed to smartsva.
#'
#' @return beefedUpKNN object.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
beefedUpKNN <- function(x, y, n.sv, k, ...){
	# Fit SVA
	svaObj <- fitSVA(x=x, y=y, n.sv=n.sv, ...)

	# Train KNN
	knnObj <- knn4(x=svaObj$corrected, y=y, k=k)

	# Return
	o <- list(sva=svaObj, knn=knnObj)
	class(o) <- "beefedUpKNN"
	o
}

#' fSVA and k-Nearest Neighbor prediction
#'
#' Predict labels or probabilities for new samples.
#'
#' @param object beefedUpKNN object.
#' @param x matrix: Normalized and log-transformed expression values (samples in rows, features in columns).
#' @param k integer: Number of nearest neighbors casting votes.
#' @param type character: Either "class" for labels or "prob" for class probabilities.
#' @param ... additional arguments, currently not used.
#'
#' @return factor if type = "class" or matrix of class probabilities if type = "posterior"
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
predict.beefedUpKNN <- function(object, x, k=NULL, type="class", ...){
	# Check k
	if(is.null(k)){
		k <- length(object$knn)
	}else{
		stopifnot(k <= length(object$knn))
	}

	# Correct with SVA
	p1 <- stats::predict(object$sva, x)

	# Predict with KNN
	p2 <- stats::predict(object$knn, p1, k=k, type=type)

	# Return
	p2
}

#' @export
print.beefedUpKNN <- function(x, ...){
	cat("# beefedUpKNN-object:\n")
	cat("Samples:", nrow(x$sva$corrected), "\n")
	cat("Features:", ncol(x$sva$corrected), "\n")
	cat("Surrogate Variables:", x$sva$n.sv, "\n")
	cat("Neighbors:", length(x$knn), "\n")
}

#' #' Plot batch-corrected centroids
#' #'
#' #' Plots the results of a fSVA and Nearest Shrunken Centroid analysis as a PCA-plot.
#' #'
#' #' @param x beefedUpSNC object.
#' #' @param threshold numeric: Amount of centroid shrinkage (controls how many features are retained in the model).
#' #' @param ... additional arguments not currently used.
#' #'
#' #' @return ggplot2 object.
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @import ggplot2
#' #' @export
#' plot.beefedUpSNC <- function(x, threshold=NULL, ...){
#' 	# Check thresholds
#' 	if(is.null(threshold)){
#' 		threshold <- x$snc$threshold
#' 	}else{
#' 		stopifnot(is.numeric(threshold))
#' 	}
#'
#' 	# Features that survive
#' 	nonzero <- pamr::pamr.predict(x$snc,
#' 																newx=t(x$sva$corrected),
#' 																threshold=x$snc$threshold,
#' 																type="nonzero")
#'
#' 	# PCA
#' 	pca <- stats::prcomp(x$sva$corrected[,nonzero], scale=TRUE)
#'
#' 	# Scores
#' 	samples <- pca$x
#' 	samples <- data.frame(samples[,1:2],
#' 												class=x$snc$y)
#' 	samples$prediction <- ifelse(samples$class == x$snc$yhat, "correct", "incorrect")
#'
#' 	# Centroids
#' 	centroids <- stats::predict(pca, t(x$snc$centroids)[,nonzero])
#' 	centroids <- data.frame(centroids[,1:2],
#' 													class=rownames(centroids),
#' 													prediction="centroid")
#'
#' 	# Assemble
#' 	P <- rbind(samples, centroids)
#' 	P$prediction <- factor(P$prediction, levels=c("centroid", "correct", "incorrect"))
#'
#' 	# Labels
#' 	vars <- summary(pca)$importance[2,1:2]
#' 	lab1 <- paste0("PC1: ", round(vars[1] * 100, digits=2), " %")
#' 	lab2 <- paste0("PC2: ", round(vars[2] * 100, digits=2), " %")
#' 	nFeatures <- paste0("Features: ", length(nonzero), " / ", ncol(x$sva$corrected))
#'
#' 	# Plot
#' 	p <- ggplot(P, aes_string(x="PC1", y="PC2",
#' 														color="class", shape="prediction")) +
#' 		geom_point(alpha=0.75) +
#' 		coord_fixed() +
#' 		scale_shape_manual(values=c(3,19,1)) +
#' 		scale_color_brewer(palette = "Set1") +
#' 		labs(title=nFeatures, x=lab1, y=lab2) +
#' 		theme_bw()
#'
#' 	# Return
#' 	p
#' }
#'
#' #' fSVA and Nearest Shrunken Centroid feature importance
#' #'
#' #' DETAILS.
#' #'
#' #' @param object beefedUpSNC object.
#' #' @param threshold numeric: Amount of centroid shrinkage (controls how many features are retained in the model).
#' #'
#' #' @return data.frame.
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @export
#' featureImportance <- function(object, threshold=NULL){
#' 	# Check thresholds
#' 	if(is.null(threshold)){
#' 		threshold <- object$snc$threshold
#' 	}else{
#' 		stopifnot(is.numeric(threshold))
#' 	}
#'
#' 	# Importance as distance to centroid
#' 	sncImp <- pamr::pamr.predict(fit=object$snc, newx=t(object$sva$corrected), type="centroid", threshold=object$snc$threshold)
#' 	sncImp <- (sncImp - object$snc$centroid.overall) / object$snc$sd
#'
#' 	# SVA statistics
#' 	svaImp <- data.frame(ProbBatch=object$sva$pprob.gam, ProbDesign=object$sva$pprob.b)
#'
#' 	# Combine
#' 	o <- data.frame(Feature=colnames(object$sva$corrected),
#' 									svaImp,
#' 									sncImp)
#'
#' 	# Return
#' 	o
#' }

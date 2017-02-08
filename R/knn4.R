#' k-Nearest Neighbor classifier with multiple k
#'
#' @inherit caret::knn3
#' @export
knn4 <- function(x, y, k){
	# All ks up to k
	o <- lapply(seq_len(k), caret::knn3, x=x, y=y)

	# Output
	class(o) <- "knn4"
	o
}

#' @export
print.knn4 <- function(x, ...){
	cat("# knn4-object:\n")
	cat("Samples:", nrow(x[[1]]$learn$X), "\n")
	cat("Features:", ncol(x[[1]]$learn$X), "\n")
	cat("Max k:", length(x), "\n")
}

#' @export
predict.knn4 <- function(object, newdata, k=NULL, type="class", ...){
	if(is.null(k)){
		k <- length(object)
	}else{
		stopifnot(k <= length(object))
	}

	stats::predict(object[[k]], newdata, type=type)
}

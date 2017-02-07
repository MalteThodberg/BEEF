caretSNC <- list(type = "Classification",
								 library = c("sva", "pamr"),
								 label="Beefed Up Shrunken Nearest Neighbor")

caretSNC$parameters <- data.frame(parameter = c("n.sv", "threshold"),
																	class = rep("integer", 2),
																	label = c("Surrogate variables", "Shrinkage threshold"))

caretSNC$grid <- function(x, y, len = NULL, search = "grid"){
	# Initial fit to obtain thresholds
	invisible(capture.output(tmp <- pamr.train(data=list(x=t(x), y=y),
																						 n.threshold=30)))

	# Expand grid
	data.frame(expand.grid(n.sv=c(0, seq_len(len)),
												 threshold=tmp$threshold))
}

caretSNC$loop <- function(grid){
	# Top levels
	loop <- subset(grid, threshold==max(threshold))

	# Submodels
	submodels <- subset(grid, threshold!=max(threshold))
	submodels <- split(submodels, submodels$n.sv)
	submodels <- lapply(submodels, subset, select="threshold")
	names(submodels) <- NULL

	# Retirn
	list(loop=loop, submodels=submodels)
}

caretSNC$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
	beefedUpSNC(x=as.matrix(x),
							y=y,
							n.sv=param$n.sv,
							threshold=param$threshold)
}

caretSNC$predict = function(modelFit, newdata, submodels = NULL) {
	out <- predict(modelFit, newdata, threshold=modelFit$tuneValue$threshold)

	if(!is.null(submodels)){
		tmp <- vector(mode = "list", length = nrow(submodels) + 1)
		tmp[[1]] <- out
		for(j in seq(along = submodels$threshold)){
			tmp[[j+1]] <- predict(modelFit, newdata, threshold=submodels$threshold[j])
		}
		out <- tmp
	}
	out
}

caretSNC$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
	out <- predict(modelFit, newdata, threshold=modelFit$tuneValue$threshold,
								 type="posterior")

	if(!is.null(submodels)){
		tmp <- vector(mode = "list", length = nrow(submodels) + 1)
		tmp[[1]] <- out
		for(j in seq(along = submodels$threshold)){
			tmp[[j+1]] <- predict(modelFit, newdata, threshold=submodels$threshold[j], type="posterior")
		}
		out <- tmp
	}
	out
}

caretSNC$sort <- function(grid){
	grid[order(grid$threshold, -grid$n.sv, decreasing=TRUE),]
}

devtools::use_data(caretSNC)

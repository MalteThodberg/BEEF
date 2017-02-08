caretKNN <- list(type = "Classification",
								 library = c("sva"),
								 label="Beefed Up k Nearest Neighbor")

caretKNN$parameters <- data.frame(parameter = c("n.sv", "k"),
																	class = rep("integer", 2),
																	label = c("Surrogate variables", "Neighbors"))

caretKNN$grid <- function(x, y, len = NULL, search = "grid"){
	data.frame(expand.grid(n.sv=c(0, seq_len(len)),
												 k=c(1,3,5,7,9,11)))
}

caretKNN$loop <- function(grid){
	# Top levels
	loop <- subset(grid, k==max(k))

	# Submodels
	submodels <- subset(grid, k!=max(k))
	submodels <- split(submodels, submodels$n.sv)
	submodels <- lapply(submodels, subset, select="k")
	names(submodels) <- NULL

	# Retirn
	list(loop=loop, submodels=submodels)
}

caretKNN$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
	beefedUpKNN(x=as.matrix(x),
							y=y,
							n.sv=param$n.sv,
							k=param$k)
}

caretKNN$predict = function(modelFit, newdata, submodels = NULL) {
	out <- predict(modelFit, newdata, k=modelFit$tuneValue$k)

	if(!is.null(submodels)){
		tmp <- vector(mode = "list", length = nrow(submodels) + 1)
		tmp[[1]] <- out
		for(j in seq(along = submodels$k)){
			tmp[[j+1]] <- predict(modelFit, newdata, k=submodels$k[j])
		}
		out <- tmp
	}
	out
}

caretKNN$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
	out <- predict(modelFit, newdata, k=modelFit$tuneValue$k,
								 type="posterior")

	if(!is.null(submodels)){
		tmp <- vector(mode = "list", length = nrow(submodels) + 1)
		tmp[[1]] <- out
		for(j in seq(along = submodels$k)){
			tmp[[j+1]] <- predict(modelFit, newdata, k=submodels$k[j], type="prob")
		}
		out <- tmp
	}
	out
}

caretKNN$sort <- function(grid){
	grid[order(grid$n.sv, grid$k),]
}

devtools::use_data(caretKNN, overwrite=TRUE)

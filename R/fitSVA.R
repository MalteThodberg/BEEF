#' Surrogate Variable Analysis (SVA) wrapper for prediction
#'
#' Wraps the smartsva and fsva functions into a standard set of fit and predict functions.
#'
#' @param x matrix: Normalized and log-transformed expression values (samples in rows, features in columns).
#' @param y factor: Class labels.
#' @param n.sv integer: Number of surrogate variables to estimate.
#' @param ... additional arguments passed to smartsva.
#'
#' @return fitSVA object (With the following components...)
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
fitSVA <- function(x, y, n.sv, ...){
	stopifnot(is.matrix(x))
	stopifnot(is.factor(y))
	stopifnot(nrow(x) == length(y))


	if(n.sv > 0){
		# Dummy design
		mod1 <- stats::model.matrix(~y)

		# Tranpose data
		xT <- t(x)

		# Fit sva
		sv <- hackedSVA(dat=xT, mod=mod1, n.sv=n.sv, ...)

		# Transform input data (From fSVA function)
		ndb <- dim(xT)[2]
		nmod <- dim(mod1)[2]
		n.sv <- sv$n.sv
		mod <- cbind(mod1, sv$sv)
		gammahat <- (xT %*% mod %*% solve(t(mod) %*% mod))[, (nmod +
																														1):(nmod + sv$n.sv)]
		db = xT - gammahat %*% t(sv$sv)
		colnames(db) <- colnames(xT)
		rownames(db) <- rownames(xT)

		# Obtain weights for transforming future samples (From fSVA function)
		wts <- (1 - sv$pprob.b) * sv$pprob.gam
		WX <- wts * xT
		svd.wx = svd(t(scale(t(WX), scale = F)))
		D <- svd.wx$d[1:n.sv]
		U <- svd.wx$u[, 1:n.sv]
		P <- t(wts * t(1/D * t(U)))

		# Sign
		sgn = rep(NA, n.sv)
		for (j in 1:sv$n.sv) {
			if (sv$n.sv > 1) {
				sgn[j] = sign(stats::cor(svd.wx$v[1:ndb, j], sv$sv[1:ndb,
																										j]))
			}
			if (sv$n.sv == 1) {
				sgn[j] = sign(stats::cor(svd.wx$v[1:ndb, j], sv$sv[1:ndb]))
			}
		}

		# Save output
		sv$P <- P
		sv$sgn <- sgn
		sv$gammahat <- gammahat
		sv$corrected <- t(db)
		sv$y <- y

	}else if(n.sv == 0){
		sv <- list(sv=NULL,
							 pprob.gam=NULL,
							 pprob.b=NULL,
							 n.sv=0,
							 P=NULL,
							 sgn=NULL,
							 gammahat=NULL,
							 corrected=x,
							 y=y)
	}else{
		stop('n.sv must be a positive integer!')
	}

	# Class and output
	class(sv) <- "fitSVA"
	sv
}

#' @export
print.fitSVA <- function(x, ...){
	message("# fitSVA-object:")
	message("Samples: ", nrow(x$corrected))
	message("Features: ", ncol(x$corrected))
	message("Surrogate Variables: ", x$n.sv)
}

#' Frozen SVA (fSVA) wrapper for prediction
#'
#' Given a SVA, batch correct new samples.
#'
#' @param object fitSVA object.
#' @param x matrix: Normalized and log-transformed expression values (samples in rows, features in columns).
#' @param ... additional arguments, not currently used.
#'
#' @return matrix with batch corrected expression for the new samples.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
predict.fitSVA <- function(object, x, ...){
	# Check dimensions match
	stopifnot(ncol(x) == ncol(object$trainData))

	if(object$n.sv > 0){
		# Transpose
		xT <- t(x)

		# Transform (from fSVA)
		newV <- object$P %*% xT
		newV <- newV * object$sgn
		newV <- t(newV)
		newV <- scale(newV) / sqrt(dim(newV)[1])
		newV <- t(newV)
		adjusted <- xT - object$gammahat %*% newV

		# Return
		o <- t(adjusted)
	}else if(object$n.sv == 0){
		o <- x
	}else{
		stop('n.sv must be a positive integer!')
	}

	# Return
	o
}

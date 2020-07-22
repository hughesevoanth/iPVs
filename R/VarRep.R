#' Variable Representation - PVA
#'
#' This function estimates how much of the total variance can be explained by each variable independently, in turn
#' @param wdata a data frame with your variable in columns and samples in rows.
#' @param cor_method the type of regression to perform
#' @keywords PVA, principal variables, variance explained
#' @export
#' @examples
#' VarRep()
VarRep <- function(wdata, cor_method = "spearman"){

	## build a correaltion matrix based on Spearman's regression
	R <- Hmisc::rcorr(as.matrix( wdata ), type = cor_method )$r
	## the sum of R-squared
	sumR2 = colSums(R*R)

	if( nrow(R) > 2) {
		
		## estimate how much of the total variance
		## can be explained by each variable - alone - in turn
		varexp = sapply(1:ncol(wdata), function(i){
			v = R[i,-i]
			remaining_R = R[-i, -i]
			remaining_R = remaining_R - (v %*% t(v))
			## estimate variance explained
			vexp = 1 - ( sum( diag(remaining_R) ) / ncol(R))
			})

		## data out  (Maximum variance explained fby each variable individually.)
		iVExp = data.frame(variable = colnames(wdata), initial_sumR2 = sumR2, VarExp_individually = varexp)

		#########################
		## Perform a proper PVA
		#########################
		## an empty object
		pva = c()
		## loop
		iterR = R
		for(i in 1:c(ncol(iterR)-1) ){
			nvar <- nrow(iterR)
			## sum of R-squared
			sumR2 = colSums(iterR*iterR)
			max_v = which(sumR2 == max(sumR2))[1] ## july 22nd
			var_name = names(sumR2)[max_v]
			## total var
			v = iterR[max_v,-max_v]
			remaining_R = iterR[-max_v, -max_v]
			remaining_R = remaining_R - (v %*% t(v))
			## estimate variance explained
			vexp = 1 - ( sum( diag(remaining_R) ) / nvar )
			## data out
			p = data.frame(var_name = var_name, vexp = vexp)
			pva = rbind( pva, p )
			## redefine iterR
			iterR = remaining_R
		}
		###  which variable remains excluded from the pva list?
		###  identify and add
		last_var = colnames(R)[!colnames(R) %in% pva$var_name]
		vexp = 1
		p = data.frame(var_name = last_var, vexp = vexp)
		pva = rbind(pva, p)
		### add the additive variance explained
		add_vexp = c( pva[1,2] , pva[2:nrow(pva),2] - pva[1:nrow(pva)-1, 2] )
		pva = data.frame(variable = pva[,1], added_vexp = add_vexp, cum_vexp = pva[,2])
		### add the Variance explained individually 
		m = match(pva$variable , iVExp$variable)
		pva = cbind(pva, iVExp[m, 2:3])
		pva = pva[, c(1,4,5,2,3)]
		} else { 
			if(nrow(R) == 2){
				vexp = R[1,2]^2
				pva = data.frame(variable = colnames(R), 
					initial_sumR2 = sumR2, VarExp_individually = vexp,
					added_vexp = c(vexp, 1-vexp), cum_vexp = c(vexp, 1) )
			}

		}

	##
	return(pva)
}




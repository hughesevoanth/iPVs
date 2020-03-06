#' A Function to extract principle variables of a data set
#'
#' This function estimates the correlation structure among two factors, using a chi-square test and the Crammer's V statsitic
#' @param responses a vector of column names to run PVA on
#' @param data a matrix of data, variables are in columns
#' @param nvarselect the number of variables to return
#' @param p.variance the proportion of variance you want you principal variable to explain
#' @param include a matrix of data, variables are in columns
#' @keywords principle variable analysis
#' @export
#' @examples
#' PVA()
PVA <- function(responses, data, nvarselect = NULL, p.variance = 1, include = NULL, 
                  ...)
  #Automatically selects variables using Principla Variables Analysis (PVA)
  #nvarselect is the number of to variables to select, counting those in include.
  #It is possible to use several criteria to control the selection:
  # 1) select all varaibles is increasing order of amount of information provided
  # 2) select nvarselect variables (this will default to the number of variables in R)
  # 3) select just enough variables, up to a maximum of nvarselect variables, to explain at least 
  #    p.variance*100 per cent of the total varaince  
{ 
  #Check response are in data
  if (!all(responses %in% names(data)))
    stop("At least one name in responses is not in data")
  
  #Get correlation matrix
  R <- rcorr(as.matrix(data[responses]))$r
  
  #Initialize
  nvar <- nrow(R)
  if (is.null(nvarselect))
    nvarselect <- nvar
  nvarselect <- min(nvarselect, nvar)
  varnotselect <- rownames(R)
  varselect <- vector(mode = "character", length = nvarselect)
  h.ordered <- vector(mode = "numeric", length = nvarselect)
  p.var <- data.frame(Variable = 1:nvarselect,
                      Cumulative.Propn = rep(0, (nvarselect))) 
  
  #Deal with variables that must be included
  ivarselect <- 0
  pinclude <- 0
  S22 <- R
  if (!(is.null(include)))
  { 
    #check include variables in responses
    if (!all(include %in% responses))
      stop("Name in include is not in responses")
    ivarselect <- length(include)
    nvarselect <- nvarselect - ivarselect
    varselect[1:ivarselect] <- include
    varnotselect <- varnotselect[!(varnotselect %in% include)]
    
    if (ivarselect == 1)
    { 
      S22 <- R[varnotselect, varnotselect] - (matrix(R[varnotselect, include], ncol=ivarselect) %*% 
                                                ginv(matrix(R[include, include], ncol=ivarselect)) %*% 
                                                matrix(R[include, varnotselect], nrow=ivarselect))
      h.ordered[1] <- sum(R[, include]*R[, include])
      p.var$Cumulative.Propn[1:ivarselect] <- 1 - (sum(diag(S22))/nvar)
    } else
    { 
      for (k in 1:ivarselect)      
      { 
        i <- match(include[k], colnames(S22))
        h.ordered[k] <-  sum((S22*S22)[,i])
        S22 <- S22update(k, i, S22, nvar)
        p.var$Cumulative.Propn[k] <- 1 - (sum(diag(S22))/nvar)
      }
    }
    pinclude <- p.var$Cumulative.Propn[ivarselect]
  }
  
  # If still have not explained enough variance or variables, select some more
  k <- 0
  if (nvarselect > 0) #still have unselected vars
  {
    if (pinclude  <= p.variance)
    { #Adjust for variables that must be included
      kvarselect <- nvarselect
      if (nvarselect+ivarselect >= nvar)
        kvarselect <- nvarselect-1
      for (k in 1:kvarselect)
      { #Select variable with maximum h
        h <- colSums(S22*S22)
        h.ordered[k+ivarselect] <- max(h)
        mh <- match(h.ordered[k+ivarselect], h)
        varselect[k+ivarselect] <- varnotselect[mh]
        varnotselect <- varnotselect[-mh]
        #Update S22
        S22 <- S22update(k+ivarselect, mh, S22, nvar)
        
        #Caclulate proportion of variance
        p.var$Cumulative.Propn[k+ivarselect] <- 1 - (sum(diag(S22))/nvar)
        if (p.var$Cumulative.Propn[k+ivarselect] >= p.variance) 
          break()
      }
      
      #Add last variable if required to exceed p.variance
      if (nvarselect+ivarselect == nvar & p.var$Cumulative.Propn[k+ivarselect] <= p.variance)
      { 
        varselect[nvar] <- varnotselect[1]
        varnotselect <- varnotselect[-1]
        h.ordered[nvar] <- h[varselect[nvar]]
        p.var$Cumulative.Propn[nvar] <- 1
      }
    }
    #crop results if required
    if (k+ivarselect < length(varselect))
    { 
      varselect <- varselect[1:(k+ivarselect)]
      h.ordered <- h.ordered[1:(k+ivarselect)]
      p.var <- p.var[1:(k+ivarselect),]
    }
  }
  
  #Finalize p.var data.frame
  nvarselect <- nrow(p.var)
  p.var <- within(p.var, 
                  { 
                    Added.Propn <- Cumulative.Propn
                    h.partial <- h.ordered
                    Selected <- factor(varselect, levels=varselect)
                  })
  if (nrow(p.var) > 1)
    p.var$Added.Propn[2:nvarselect] <- p.var$Added.Propn[2:nvarselect] - p.var$Added.Propn[1:(nvarselect-1)]
  p.var <- p.var[c("Variable","Selected", "h.partial", "Added.Propn", "Cumulative.Propn")]
  
 
  return(p.var)
}


#' Function that adds a variable and calcullates the partial correlation matrix
#'
#' This function estimates partial correlation matrix
#' @param k i do not know what this is, it comes from the PV()
#' @param i i do not know what this is, it comes from the PV()
#' @param S22 S22 to be updated
#' @param nvar the number fo variables
#' @keywords PVA
#' @export
#' @examples
#' S22update()
S22update <- function(k, i, S22, nvar)
#
#Called by PVA and PVA.manual
{ #Update S22
  si <- S22[i, i]
  if (k < (nvar - 1))
  { 
    S21 <- S22[i, -i]
    S22 <- S22[,-i] 
    S22 <- S22[-i, ]
  }
  else
  { 
    if (k == (nvar-1))
    { 
      noti <- c(1,2)[!(i == c(1,2))]
      S21 <- S22[i, noti]
      S22 <- S22[noti,noti]
    }
  }
  S22 <- S22 - (S21 %*% t(S21))/si
  return(S22)
}

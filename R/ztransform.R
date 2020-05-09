#' A Function to Z transform a vector of data
#'
#' A Function to Z transform a vector of data
#' @param data the vector of data to be transformed
#' @export
#' @examples
#' ztransform()
ztransform = function (data) {
  m = mean(data, na.rm = TRUE)
  s = sd(data, na.rm = TRUE)
  out = (data - m)/s
  return( out )
}


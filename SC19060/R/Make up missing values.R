#' @title Make up missing values 
#' @name make up missing values
#' @description Make up missing values by median
#' @param data the matrix may be has missing values
#' @return a matrix has no missing values
#' @examples
#' \dontrun{
#' x<-rnorm(20)
#' data<-rnorm(x,4,5)
#' data[2,3]<-NA
#' MakeupMissingValue(data)
#' }
#' @export
MakeupMissingValue<-function(data)
{
  for(i in seq(ncol(data)))
    if(any(a<-is.na(data[,i])))
    {
      data[a,i]<-median(data[,i],na.rm = T)
    }
  data
}
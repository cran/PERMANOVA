BinaryVectorCheck <- function (x)
{
  for (i in 1:length(x)) if (!((x[i] == 0) | (x[i] == 1)))
    return(FALSE)
  return(TRUE)
}

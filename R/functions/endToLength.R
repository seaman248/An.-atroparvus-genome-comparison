#coords(data.frame): start, end
endToLength <- function(coords){
  length <- coords[,2] - coords[,1]
  return(length)
}
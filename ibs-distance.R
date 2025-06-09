# m is a marker matrix with samples as rows and columns as marker
# that has genotype encoded as 0, 1, 2 with hom1, het, hom2
ibs.dist <- function(m, center=T){
  if(center){
    m <- m-1
  }
  m1 <- m+1
  m1[m!=0] <- 0
  cp <- tcrossprod(m)
  a1 <- tcrossprod(m1)
  a2 <- cp + a1
  1-(a2/ncol(m)+1)/2
}

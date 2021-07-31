


#==========================================================================#
# Soft.Threshold Function
#==========================================================================#
Soft.Threshold <- function(w,lam){
  # soft-thresholding function
  if(w <= -lam){
    res <- w+lam
  }else if(w >= lam){
    res <- w-lam
  }else{
    res <- 0
  }
  return(res)
} # Soft.Threshold(w,lam)








# Wrapper function to invoke C functions at the R shell

dyn.load("Gravity_new.so")
Gravity_new <- function(model_params, misfit) 
{
	misfit = 0
  	res <- .C("Gravity_new", as.double(model_params), as.double(misfit))
	return(as.numeric(res[2]))
}

dyn.load("Magneto_new.so")
Magneto_new <- function(model_params, misfit) 
{
	misfit = 0
  	res <- .C("Magneto_new", as.double(model_params), as.double(misfit))
	return(as.numeric(res[2]))
}

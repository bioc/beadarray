setGeneric("backgroundCorrect", function(object, method="subtract", offset=0) standardGeneric("backgroundCorrect"))

setMethod("backgroundCorrect" ,"BeadLevelList", function(object, method = "subtract", offset = 0)
{
     
    method <- match.arg(method, c("none", "subtract", "half", 
        "minimum",normexp))
    
    switch(method, subtract = {
        object@G <- object@G - object@Gb
        
    }, half = {
        object@G <- pmax(object@G - object@Gb, 0.5)
        
    }, minimum = {
        object@G <- as.matrix(object@G - object@Gb)
        
for (slide in 1:ncol(object@G)) {
            i <- object@G[, slide] < 1e-18
            if (any(i, na.rm = TRUE)) {
                m <- min(object@G[!i, slide], na.rm = TRUE)
                object@G[i, slide] <- m/2
            }
            
        }
    }
    ,	normexp={
	for (j in 1:ncol(object@G)) {
		x <- object@G[,j]-object@Gb[,j]
		out <- normexp.fit(x)
#		if(verbose) cat("G: bg.bias=",out$par[1]," bg.sd=",exp(out$par[2])," fg.mean=",exp(out$par[3]),"\n",sep="")
		object@G[,j] <- normexp.signal(out$par,x)

	}}       


           )

    if (offset) {
        object@G <- object@G + offset
        
    }
    object

}
)

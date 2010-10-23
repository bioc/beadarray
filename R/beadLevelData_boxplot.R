setMethod("boxplot",
    signature(x = "beadLevelData"),
    function (x, transFun=logGreenChannelTransform,...) 
    {
       tmp = list()
  	arraynms = sectionNames(x)
  	narrays = length(arraynms)


  	for(i in 1:narrays){
      	tmp[[arraynms[i]]] = transFun(x, array=i)
	}
  	boxplot(tmp,...)   
    }
)


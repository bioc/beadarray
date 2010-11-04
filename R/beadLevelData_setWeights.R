setWeights = function(BLData, wts, array, combine = FALSE, wtName = "wts"){

	if (length(array) == 1)
	{
		##one array
		if(combine)
		{
			#check for existing weights
			if(!wtName %in% colnames(BLData[[array]]))

			{
				##Try and insert the weights into the beadLevelData
	
				BLData = insertBeadData(BLData, array=array, what=wtName, data = wts)
			}
			else
			{
				wtCol = grep(wtName, colnames(BLData[[array]]))

				##Try and insert the weights into the beadLevelData

				if(length(wtCol) > 0) combinedWts = pmin(wts,BLData[[array]][,wtCol])
				else combinedWts = wts

				BLData = insertBeadData(BLData, array=array, what=wtName,data = combinedWts)
			}
		}
		else
		{
			BLData = insertBeadData(BLData, array=array, what=wtName, data = wts)
		}
	}
	else
	{
		##multiple arrays
		if(combine)
		{
			for(i in array)
			{
				#check for existing weights
				if(!wtName %in% colnames(BLData[[array]]))
				{BLData = insertBeadData(BLData, array=i, what=wtName, data = wts[[i]])}
				else{
					wtCol = grep(wtName, colnames(BLData[[array]]))

					if(length(wtCol) > 0) combinedWts = pmin(wts,BLData[[array]][,wtCol])
					else combinedWts = wts

					BLData = insertBeadData(BLData, array=i, what=wtName, data = combinedWts)
				}
			}		
		}
		else
		{
			for(i in array)
			{
				BLData = insertBeadData(BLData, array=i, what=wtName, data = wts[[i]])
			}
		}
	}

BLData

}

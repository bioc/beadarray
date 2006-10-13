DiffScore = function(BSData, QC=NULL, cond, ref){

if(is.null(QC)) QC = QCInfo(BSData)


if(length(cond)==1){
  Avg_Intensity_Ref = exprs(BSData)[,ref]
  Avg_Intensity_Cond = exprs(BSData)[,cond]

  BeadStDev_Ref = se.exprs(BSData)[,ref]
  BeadStDev_Cond = se.exprs(BSData)[,cond]

  SigmaRef = QC$StDev[ref,11]
  SigmaCond = QC$StDev[cond,11]
}


else{
Avg_Intensity_Ref = apply(exprs(BSData)[,ref],1,mean)
Avg_Intensity_Cond = apply(exprs(BSData)[,cond], 1, mean)

SigmaRef=sd(QC$Signal[ref,11])/sqrt(length(ref))
SigmaCond=sd(QC$Signal[ref,11])/sqrt(length(ref))

BeadStDev_Ref = apply(se.exprs(BSData)[,ref],1,mean)
BeadStDev_Cond = apply(se.exprs(BSData)[,cond],1,mean)


}



  lm.ref = lm(BeadStDev_Ref~Avg_Intensity_Ref)
  lm.cond = lm(BeadStDev_Cond~Avg_Intensity_Cond)

  Aref = lm.ref$coefficients[1]
  Bref = lm.ref$coefficients[2]

  Acond = lm.cond$coefficients[1]
  Bcond = lm.cond$coefficients[2]


if(length(ref)>1){

sRef = apply(cbind(BeadStDev_Ref, Aref + Bref*Avg_Intensity_Ref), 1, max)
sCond = apply(cbind(BeadStDev_Cond, Acond + Bcond*Avg_Intensity_Cond), 1, max)

Denom1 = (sRef^2 + SigmaRef^2)/length(ref)
Denom2 = (sCond^2 + SigmaCond^2)/length(cond)

}

else{  

Denom1 = (Aref +  (2.5*Bref*Avg_Intensity_Ref))^2 + SigmaRef^2

Denom2 = (Acond +  (2.5*Bcond*Avg_Intensity_Cond))^2 + SigmaCond^2

}


FullDenom = sqrt(Denom1 + Denom2)

Numerator = abs(Avg_Intensity_Cond - Avg_Intensity_Ref)

Diff_Pvalue = 2*pnorm(Numerator/FullDenom, lower.tail=FALSE)



DiffScore = -10*sign(Avg_Intensity_Cond - Avg_Intensity_Ref)*log10(Diff_Pvalue)

DiffScore

}

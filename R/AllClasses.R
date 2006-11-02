setClassUnion("AssayDataOrNull", c("AssayData", "NULL"))
setClassUnion("matrixOrNull", c("matrix", "NULL"))

setClass("ExpressionSetIllumina",
         representation(QC = "AssayDataOrNull"),
         contains="eSet"

)

setClass("BeadLevelList",

	representation(
	G="matrix",
	R="matrixOrNull",
	Rb="matrixOrNull",
	Gb="matrix",
	GrnY ="matrix",
	GrnX="matrix",
	RedX="matrixOrNull",
	RedY="matrixOrNull",
	ProbeID="matrix",
	targets ="data.frame"
)

)

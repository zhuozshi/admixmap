# Admixmap
--------------------

## Examples on UKBB

Run one chromosome for a specific trait by 

`python computeChrom.py [-o output path] [-p phenotype and covariates file path] [-r rfmix file path] 
`

Please note that phenotype and covariance should be in the same -p file. See file header for details.

--------------------

After finishing all chromosomes for one trait, compile by running

`python compileResult.py [-o outputDirectory] [-i inputDirectory] [-t traitIndex] [-p phenotypeFilePath]
`

After you finish all the trait, compile the summary table by

`python tableCompile.py [-o outputDirectory] [-n outputFileName] [-s areTableRowsSpreadInDirecotries] [-i inputDirectory] [-p phenotypeFilePath]
`




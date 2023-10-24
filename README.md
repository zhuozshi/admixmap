# Admixmap
--------------------

## Examples on UKBB

### Run all traits 
Run 22 chromosomes for one trait and generate all the data at a time

`python all_script_ukbb.py [-o outputDirectory] [-t traitIndex] [-p phenotypeFilePath] [-l lancFilePath] [-c covarianceFileDirectory] [-n numberOfAncestry]
`

After you finish all the trait, compile the summary table by

`python tableCompile_ukbb.py [-o outputDirectory] [-n outputFileName] [-s areTableRowsSpreadInDirecotries] [-i inputDirectory] [-p phenotypeFilePath]
`

Traits to summary can be given by a list in tableCompile_ukbb.py

### Run separately
Run one chromosome for a specific trait by 

`python computeChrom_ukbb.py [-o outputDirectory] [-t traitIndex] [-c chromosome] [-p phenotypeFilePath] [-l lancFilePath] [-v covarianceFileDirectory] [-n numberOfAncestry]
`

After finishing all chromosomes for one trait, compile by running

`python compileResult_ukbb.py [-o outputDirectory] [-i inputDirectory] [-t traitIndex] [-p phenotypeFilePath]
`

After you finish all the trait, compile the summary table by

`python tableCompile_ukbb.py [-o outputDirectory] [-n outputFileName] [-s areTableRowsSpreadInDirecotries] [-i inputDirectory] [-p phenotypeFilePath]
`




# Admixmap
--------------------

## Examples on UKBB

Run one chromosome for a specific trait by 

`python computeChrom.py [-o output path] [-p phenotype and covariates file path] [-r rfmix file path] 
`

Please note that phenotype and covariance should be in the same -p file. [Make sure to see file header for details.](computeChrom.py)

--------------------

After finishing all chromosomes for one trait, compile by running

`python compileResult.py [-o output path] [-i input path] [-n name displayed in figure title] 
`

[Make sure to see file header for details.](compileResult.py)

--------------------

After you finish all the trait, compile the summary table by

`python tableCompile.py [-o output path] [-i input path] [-n name of compressed csv table inside the zip file] 
`

[Make sure to see file header for details.](tableCompile.py)

--------------------

Check this [example for usage.](Example.ipynb)


--------------------

## Calculate Common Association (Harmonizing)

First generate template with chromosome locations

`python commonAssoc.py -t [assocFile 1] ... [assocFile n] [output] 
`

--------------------

Interpolate association result using template

`python commonAssoc.py -a [assocFile] [template] [output] 
`

[Make sure to see file header for details.](commonAssoc.py)


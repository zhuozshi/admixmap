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

## Calculate Common Association -- Harmonizing

First generate template with chromosome locations

`python commonAssoc.py -t [assocFile 1] ... [assocFile n] [output] 
`

--------------------

 Interpolate association result using template

`python commonAssoc.py -a [assocFile] [template] [output] 
`

[Make sure to see file header for details.](commonAssoc.py)


--------------------

## Get Peaks

Get peaks from association files directory

`python getPeaks.py [assocFiles Directory] [bandsDict] [strict p value threshold] [loose p value threshold] [output] 
`

[Make sure to see file header for details.](getPeaks.py)

--------------------

## Compare results to published data

This script compares hits from an analysis with published GWAS data. It has two steps:

1. Take in a "reference dataset" (i.e. genome-wide significant hits), and find all the variants in LD with data in the reference from LD Proxy.

`python compare_signals.py pull_ld -r [reference variants] -rc [chr column in -r] -rp [pos column in -r] -rs [delimiter in -r] -g [genome build] -anc [ancestries] -t [LD link TOKEN] -lo [LD data output path]`

* Defaults
 * `-rc`: CHR
 * `-rp`: POS
 * `-rs`: ,
 * `-g`: grch38_high_coverage
 * `-anc`: EUR
  * `anc` can take any value within `EUR`, `AFR`, `AMR`, `SAS`, `EAS`, `ALL`, `EACH`. If `EACH` is specified, data from each of the other `anc` arguments will be pulled.

2. Read in the generated results, and identify which of the generated results are in the "reference dataset" or are in LD with one of the LD Proxy outputs.

`python compare_signals.py compare -a [association file] -ac [chr column in -a] -ai [SNP id column in -a] -ap [pos column in -a] -as [delimiter in -a] -ae [end pos column in -a] -li [LD data input path] -r2 [r2 cutoff] -anc [ancestries] -o [output]

* Defaults
 * `-ac`: CHR
 * `-ap`: POS
 * `-as`: ,
 * `-anc`: EUR
 * `-r2`: 0.8
  * `anc` can take any value within `EUR`, `AFR`, `AMR`, `SAS`, `EAS`, `ALL`, `EACH`. If `EACH` is specified, data from each of the other `anc` arguments will be pulled.

* This script can take in association files in two formats. The default is for variants from an association analyses, where the chromosome and position of the variant are required. The second is for regions from an association analyses, where both the chromosome, position, and end position (`-ae`) is required.

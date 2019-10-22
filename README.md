
# Joint Likelihood Mapping (JLIM) 

JLIM is a cross-trait test of shared causal effect, which is described ([Chun et al. Nature Genetics 2017](https://www.nature.com/articles/ng.3795)). JLIM tests whether two traits – primary and secondary – are driven by shared causal effect or not. Typically, the primary trait is a large GWAS study, and the secondary trait can have a smaller sample size. For the main trait, JLIM takes only summary-level association statistics, but for secondary trait, it requires genotype-level data to generate a permutation-based null distribution. JLIM is simultaneously released at [Cotsapas lab github](https://github.com/cotsapaslab/jlim) and [Sunyaev lab website](http://genetics.bwh.harvard.edu/wiki/sunyaevlab/jlim). 

JLIM 2.0 

Has the added functionality of obtaining association statistics for the secondary trait in a cohort specific manner and then combining them before running JLIM. JLIM 2.0 is described in () 

## How to install

The core JLIM module is implemented as an R extension (**jlimR**). **jlimR** depends on **getopt** module. If it is not installed, **getopt** can be installed from CRAN by:
JLIM 2 has an added pre-porcessing step coded in pyhton 2.

```
Rscript -e 'install.packages("getopt", repos="http://cran.r-project.org")' 
```

After **getopt** has been installed, **jlimR** (included in the distribution file) can be installed by:
 
```
tar -zxvf JLIM_2.0.tar.gz 

cd JLIM_2.0

R CMD INSTALL JLIM_2.0.tar.gz
```

In case that it is preferred to install R extensions in your home directory (e.g. ~/R) instead of the default system path, please do the following instead: 
```
Rscript -e 'install.packages("getopt", "~/R", repos="http://cran.r-project.org")' 

cd JLIM_2.0

R CMD INSTALL -l ~/R JLIM_2.0.tar.gz

```
And then, add your local R library path to **R_LIBS** environment variable in .bashrc or .profile as:
```
export R_LIBS=~/R:$R_LIBS
```

The python based pre-processing step for JLIM 2 depends on numpy and scipy. Please check that your versions of the following packages are at least:

python 2.7.9
numpy 1.14.3
scipy 1.0.0


## How to run JLIM on provided example  

### Example data

In JLIM_2.0/example, we provide the following simulated dataset: 
Data from three simulated cohorts (A,B,C) in chromosome 1- Bimbam and map files for each cohort. A peaks file containing the mid point of each region in chromosome 1 that will be analyzed. A reference LD file corresponding to the two loci analyzed. An IndexSNP file (Height/Height_indexSNP.tsv) with both loci under analysis and two primary trait statistics files (Height) containing only summary statistics. We also include a phenotype file, a samples file and a covariate file for each cohort.

To run the example after JLIM is installed, run all commands in the file:

example/CommandsExample.sh

## Files needed to run JLIM

### bimbam file
-bimbam file
One genotype file in BIMBAM format is needed per cohort, per chromosome analyzed, for the secondary trait. These files should have no header and should be gzipped. The separator used in the file should be provided (see bimbam.separator). **Currently, JLIM does not allow missing genotypes.** 

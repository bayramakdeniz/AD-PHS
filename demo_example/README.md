# phs
Polygenic hazard scores for onset prediction

The main script is ``phs_generate_new.m``.
This will calculate a PHS model for a given genotype/phenotype dataset and save it as an output. 

To see how it works, run ``phs_generate_new_demo.m``.
It will use some random data saved in an example ``.mat`` file, and generate some plots.

Note that ``phs_generate_new_demo.m`` assumes an absolute incidence curve
(in fact, a population age-specific incidence for prostate cancer)
and generates some plots to demonstrate the effect of the PHS on age of prostate cancer onset.
The genotype data in the example are not real, though, so the effect looks weak.
The point is just to demonstrate how one can show the model's application to absolute risk.
The important part of this script is just to show how to call ``phs_generate_new.m``.

Note that when you run ``phs_generate_new``, you will not need the trend test portion that is meant for selecting SNPs from the larger dataset.
You can just set the ``pthresh`` input to ``1``, which will include all SNPs.
Using only ~60 SNPs, it should not take long to run. 
``pthresh2`` is an arbitrary p-value threshold to see if a given SNP adds value to the model.
This may be a useful step. Alternatively, we can just assume that all of the SNPs are worth including (i.e., set ``pthresh2`` to 1). 
An alternative to 3 would be to just use lines 132-136 of ``phs_generate_new``.
Those lines are really all you need to calculate the PHS weights for your SNPs. 

# Loading plink  genotypes into MATLAB

This can be done with ``PlinkRead_binary2.m``
To test out use a synthetic dataset (Download mostest_demo.tar.gz file from here https://1drv.ms/u/s!Ai1YZmdFa9ati40Inztrv_4erqcdWw?e=ixWDUe )

```
nsubj = 10000;
nsnps = 149454;
bfile = 'chr21';
geno_int8 = PlinkRead_binary2(nsubj, 1:nsnps, bfile);
geno = double(geno_int8);

```

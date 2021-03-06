Source code for "Cell cycle-associated expression patterns predict 
gene function in mycobacteria". Bandekar et al (2020). Current Biology.
-----------------------------------------------
copyright 2020, Thomas R. Ioerger (ioerger@cs.tamu.edu)


Curve-fitting of expression profiles with sinusoidal and quadratic 
curves to identify periodic genes.
------------------------------------------------------------------

- script: curvefit2.py
- dependencies: scipy and statsmodels (python packages)
- usage: `python curvefit2.py <normalized_expression_file> [cos|rv]` 
- input files: 
  - normalized_expression_file (![total_deseq_norm_skip0hr.txt](total_deseq_norm_skip0hr.txt)) - normalized expression values of each gene at each timepoint in each replicate of each strain)
  - corr.txt (contains correlation coefficients between cos and rv for all genes) (script automatically reads this)
  - H37Rv3.prot_table (annotation, info on ORFs) (script automatically reads this)
  - H37Rv.COG_roles.dat (functional categories for each genes) (script automatically reads this)
- output: fits curves for all genes to cos or rv cells; saves parameter estimates and goodness-of-fit in a tab-separated spreadsheet
- runtime: about 15 seconds
- Example: `python curvefit2.py total_deseq_norm_skip0hr.txt cos > curvefit2_cos.txt`
- columns in output file:
  -  orf ID
  -  gene Name
  -  COG category
  -  mean expression (over 15 timepoints, 2 reps each)
  -  standard deviation
  -  5 parameters for sin: amplitude A, frequency B, phase C, trend D, const E (formula: `A*sin(2.*PI*B*(x/55.)+C)+x*D+E)`
  -  RSS_sin - residual sum-of-squares with respect to sin fit
  -  correlation coefficient of predicted vs actual expression values
  -  period - 55/freq
  -  3 parameters for quadratic: A, B, C (`Ax^2+Bx+C`)
  -  RSS_quad - residual sum-of-squares with respect to quadratic fit
  -  RSS_sin/RSS_quad - if this ratio is >0.45, genes are considered periodic
- note: this script only generates a spreadsheet, not images, but if you plot the two fits, they look like this...

![](Rv0001_cos.png)



Clustering of Expression Profiles
---------------------------------

- script: cluster.R
- dependency: R, and gplots library
- usage: `Rscript cluster.R <inputfile> <output_spreadsheet> <output_pdf> <K>`
  - input file: spreadsheet with curvefit means of expression for each each gene at each timepoint
  - output file: tab-separated spreadsheet (.txt) with cluster numbers for each gene
  - K: number of clusters (integer)
- example: 
<BR>`Rscript cluster.R curvefit2_cos_means.txt curvefit2_cos_clust.txt curvefit2_cos_clusters.pdf 8`
<BR>output files: [curvefit2_cos_clust.txt](curvefit2_cos_clust.txt), [curvefit2_cos_clusters.pdf](curvefit2_cos_clusters.pdf)

![](cluster1.png)




Gaussian-Process Smoothing of expression profiles
-------------------------------------------------

- script: fitGP.py
- dependency: GPy python package (https://gpy.readthedocs.io/en/deploy/)
- input file: total_deseq_norm_skip0hr.txt (script automatically reads this)
- output: creates a .png file with plot of expression for a gene over 54hr, with a trend curve fit by a Gaussian Process
- usage: `python fitGP.py <orfID> <cos|rv> <output.png>`
  - orf IDs must correspond to genes in the annotation of M. tuberculosis H37Rv, NC_000962.3
  - select 'cos' if you want to see expression in cold-sensitive DnaA-mutant, or 'rv' for H37Rv control culture
- example: `python fitGP.py Rv1907c cos Rv1907c_cos.png`
- runtime: about 15 seconds

![](Rv1907c_cos.png)



Cell-cycle Modeling
-------------------

- script: sim4.py
- input file: OD.txt (contains OD and FhaA measurements at various timepoints; script automatically reads this)
- output: creates temp.png with a plot of curves for simulated OD, biomass, FhaA levels, etc.
- usage: `python sim4.py`
- runtime: about 15 seconds
- parameters: if users wish to try different parameters, they can edit the last line of the script, which calls sim()

- Example: 
```
> python sim4.py
Cycle=35(1), Shift=10, Recov=6, FHA=0.8,
ChromDupInit=20, DupWin=12+~7, GrowthRate=0.07, Delay=0.3
cc_od = 0.9195
cc_oriter = 0.9607
cc_fha = 0.9577
(see temp.png)
```

![](sim4.png)


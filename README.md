Source code for "Cell cycle-associated expression patterns predict 
gene function in mycobacteria". Bandekar et al (2020). Current Biology.
-----------------------------------------------
copyright 2020, Thomas R. Ioerger (ioerger@cs.tamu.edu)


Gaussian-Process Smoothing of expression profiles
-------------------------------------------------

* script: fitGP.py
* dependency: GPy python package (https://gpy.readthedocs.io/en/deploy/)
* input file: total_deseq_norm_skip0hr.txt (script automatically reads this)
* output: creates a .png file with plot of expression for a gene over 54hr, with a trend curve fit by a Gaussian Process
* usage: python fitGP.py <orfID> <cos|rv> <output.png>
 *  (orf IDs must correspond to genes in the annotation of M. tuberculosis H37Rv, NC_000962.2)
 *  (select 'cos' if you want to see expression in cold-sensitive DnaA-mutant, or 'rv' for H37Rv control culture)
* example: python fitGP.py Rv1907c cos Rv1907c_cos.png
* runtime: about 15 seconds



Cell-cycle Modeling
-------------------

script: sim4.py
input file: OD.txt (contains OD and FhaA measurements at various timepoints; script automatically reads this)
output: creates temp.png with a plot of curves for simulated OD, biomass, FhaA levels, etc.
usage: python sim4.py
runtime: about 15 seconds
parameters: if users wish to try different parameters, they can edit the last line of the script, which calls sim()

Example: 
> python sim4.py
Cycle=35(1), Shift=10, Recov=6, FHA=0.8,
ChromDupInit=20, DupWin=12+~7, GrowthRate=0.07, Delay=0.3
cc_od = 0.9195
cc_oriter = 0.9607
cc_fha = 0.9577
(see temp.png)



Curve-fitting of expression patterns with sinusoidal and quadratic 
curves to identify periodic genes.
------------------------------------------------------------------

script: curvefit2.py
dependenceies: scipy and statsmodels (python packages)
usage: python curvefit2.py [cos|rv] - fit curves for all genes to cos or rv cells
input file: corr.txt (contains correlation coefficients between cos and rv for all genes; script automatically reads this)
runtime: about 15 seconds
output: tab-separated spreadsheet

columns:
  orf ID
  gene Name
  COG category
  mean expression (over 15 timepoints, 2 reps each)
  standard deviation
  5 parameters for sin: amplitude A, frequency B, phase C, trend D, const E (formula: A*sin(2.*PI*B*(x/55.)+C)+x*D+E)
  RSS_sin - residual sum-of-squares with respect to sin fit
  corr - correlation coefficient of predicted vs actual expression values
  period - 55/freq
  3 parameters for quadratic: A, B, C (Ax^2+Bx+C)
  RSS_quad - residual sum-of-squares with respect to quadratic fit
  RSS_sin/RSS_quad - if this ratio is >0.45, genes are considered periodic

Example:
> python curvefit2.py cos > curvefit2_cos.txt

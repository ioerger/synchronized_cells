Source code for paper "Cell cycle-associated expression patterns
predict gene function in mycobacteria", Bandekar et al (2020). Current
Biology.

Gaussian-Process Smoothing of expression profiles
-------------------------------------------------

script: fitGP.py
dependency: GPy python package (https://gpy.readthedocs.io/en/deploy/)
input file: total_deseq_norm_skip0hr.txt (script automatically reads this)
output: creates a .png file with plot of expression for a gene over 54hr, with a trend curve fit by a Gaussian Process
usage: python fitGP.py <orfID> <cos|rv> <output.png>
  (select 'cos' if you want to see expression in cold-sensitive DnaA-mutant, or 'rv' for H37Rv control culture)
example: python fitGP.py Rv1907c cos Rv1907c_cos.png
runtime: about 15 seconds



Cell-cycle Modeling
-------------------

script: sim4.py
input file: OD.txt (script automatically reads this)
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




import sys
import numpy as np
import GPy
import scipy
from matplotlib import pyplot as plt

Rv = sys.argv[1]

hours = [3,6.5,9,12,18.5,21,27,31,33,36,39.5,42,45.5,52,55]
times = hours+hours+hours+hours

for line in open("total_deseq_norm_skip0hr.txt"):
  if Rv in line:
    w = line.rstrip().split('\t')
    vals = [float(x) for x in w[1:]]

######################

titlesize = 36
labelsize = 28
ticsize = labelsize

params = {
     'xtick.labelsize':ticsize,
     'ytick.labelsize':ticsize}
plt.rcParams.update(params)

def plot_gp(X, m, C, title,training_points=None):
    """ Plotting utility to plot a GP fit with 95% confidence interval """
    # Plot 95% confidence interval 
    plt.fill_between(X[:,0],
                     m[:,0] - 1.96*np.sqrt(np.diag(C)),
                     m[:,0] + 1.96*np.sqrt(np.diag(C)),
                     alpha=0.5)
    # Plot GP mean and initial training points
    plt.plot(X, m, "-")
    #plt.legend(labels=["GP fit"])
    #
    plt.xlabel("Time (hrs)",fontsize=labelsize)
    plt.ylabel(title,fontsize=labelsize)
    #
    # Plot training points if included
    if training_points is not None:
        X_, Y_ = training_points
        plt.plot(X_, Y_, "kx", mew=2)
        #plt.legend(labels=["GP fit", "sample points"])

######################

Xnew = np.linspace(0 , 55 , 100)[:, None]
k = GPy.kern.RBF(1, variance=1., lengthscale=1., name="rbf")

def eval_ll(m,X,Y):
  means,vars = m.predict(X)
  stdevs = np.sqrt(vars)
  ll = scipy.stats.norm.logpdf(Y,means,stdevs)
  vals = [[X[i],Y[i],means[i],stdevs[i],ll[i]] for i in range(len(X))]
  for x in vals: print '\t'.join([str(y) for y in x])
  print("sum_ll=%s" % sum(ll))
  return sum(ll)[0]

######################

if "cos" in sys.argv: X,Y = np.array([[x] for x in times[:30]]),np.array([[x] for x in vals[:30]])
elif "rv" in sys.argv: X,Y = np.array([[x] for x in times[30:]]),np.array([[x] for x in vals[30:]])
else: print "error: give 'cos' or 'rv' as second arg"; sys.exit(0)

title = Rv
yaxis = "Expression (RPKM)"
fname = sys.argv[3]

for x,y in zip(X,Y): print x,y

m = GPy.models.GPRegression(X, Y, k)
m.optimize()
print(m)
ll1 = eval_ll(m,X,Y)

# Use GPy model to calculate the mean and covariance of the fit at Xnew
mean, Cov = m.predict_noiseless(Xnew, full_cov=True)

# Plot the GP fit mean and covariance
plt.figure(figsize=(14, 8))
plot_gp(Xnew, mean, Cov, yaxis, training_points=(X,Y))
plt.title(title,fontsize=titlesize)
plt.savefig(fname)
#plt.show()


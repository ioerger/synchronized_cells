import numpy,scipy.optimize,scipy.stats
import sys,math,random
import numpy as np
import statsmodels.stats.multitest

def stats(vals):
  sum,ss = 0,0
  for x in vals: sum += x; ss += x*x
  N = float(len(vals))
  mean = sum/N
  var = ss/N-mean*mean
  if var<0: var = 0
  stdev = math.sqrt(var)
  return mean,stdev

correlations = {}
for line in open("corr.txt"):
  w = line.rstrip().split('\t')
  correlations[w[0]] = float(w[1])

genenames = {}
for line in open("H37Rv3.prot_table"):
  w = line.rstrip().split("\t")
  genenames[w[8]] = w[7]

COG = {}
for line in open("H37Rv.COG_roles.dat"):
  w = line.rstrip().split('\t')
  COG[w[0]] = w[3]

PI = 3.1415927
TP = [3., 6.5, 9., 12., 18.5, 21., 27., 31., 33., 36., 39.5, 42., 45.5, 52., 55.]

def F(x,ampl,freq,phas,trend,const):
  return ampl*numpy.sin(2.*PI*freq*(x/55.)+phas)+x*trend+const

def Q(x,a,b,c): return a*x*x+b*x+c

###############################

if len(sys.argv)<2 or sys.argv[1] not in ["cos","rv","permute"]:
  print("usage: python curvefit2.py [cos|rv|permute]")
  sys.exit(0)
series = sys.argv[1]

data = []
skip = 1
for line in open("total_deseq_norm_skip0hr.txt"):
  if skip>0: skip -= 1; continue
  w = line.split()
  data.append(w)
N = len(data)

cos1 = [[float(w[i]) for i in range(1,16)] for w in data]
cos2 = [[float(w[i]) for i in range(16,31)] for w in data]
rv1 = [[float(w[i]) for i in range(31,46)] for w in data]
rv2 = [[float(w[i]) for i in range(46,61)] for w in data]

# standard-normalize expr levels
cos1norm,cos2norm,rv1norm,rv2norm = [],[],[],[]
for i in range(len(data)):
  Y1,Y2 = cos1[i],cos2[i]
  m,s = stats(Y1+Y2)
  cos1norm.append([(y-m)/s for y in Y1])
  cos2norm.append([(y-m)/s for y in Y2])
  Y1,Y2 = rv1[i],rv2[i]
  m,s = stats(Y1+Y2)
  rv1norm.append([(y-m)/s for y in Y1])
  rv2norm.append([(y-m)/s for y in Y2])

if series=="permute":
  Ng,Nt = len(data),15
  for i in range(Ng):
    for j in range(Nt):
      # swap both replicates of cos
      temp1,temp2 = cos1norm[i][j],cos2norm[i][j]
      p,q = random.randint(0,Ng-1),random.randint(0,Nt-1)
      cos1norm[i][j],cos2norm[i][j] = cos1norm[p][q],cos2norm[p][q]
      cos1norm[p][q],cos2norm[p][q] = temp1,temp2

########################

results = []
for i,w in enumerate(data):
  orf = w[0]
  cog = COG.get(orf,'?')
  gene = genenames.get(orf,"?")
  #sys.stderr.write(orf+"\n")

  if correlations[orf]>0.9: 
    sys.stderr.write("skipping %s because corr(cos,rv)>0.9\n" % orf)
    continue

  X = numpy.array(TP)

  if series=="cos" or series=="permute": 
    Y1,Y2 = cos1norm[i],cos2norm[i]
    m,s = stats(cos1[i]+cos2[i])
  else: 
    Y1,Y2 = rv1norm[i],rv2norm[i]
    m,s = stats(rv1[i]+rv2[i])

  # if series=="cos":
  #  Y1 = [float(w[i]) for i in range(1,16)]
  #  Y2 = [float(w[i]) for i in range(16,31)]
  #else:
  #  Y1 = [float(w[i]) for i in range(31,46)]
  #  Y2 = [float(w[i]) for i in range(46,61)]
  #
  #m,s = stats(Y1+Y2)
  #Y1 = [(y-m)/s for y in Y1]
  #Y2 = [(y-m)/s for y in Y2]

  Y = [(y1+y2)/2. for y1,y2 in zip(Y1,Y2)]
  delta = [0.25*(y1-y2)**2 for y1,y2 in zip(Y1,Y2)]
  #delta = [abs (y1 - y2) for y1, y2 in zip(Y1, Y2)]
  ss = sum(delta)-max(delta)
  #absdelta = [abs (y1 - y2) for y1, y2 in zip(Y1, Y2)]
  kurt = [(y1 - y) ** 4 for y1,y in zip(Y1,Y)]
  krt = sum(kurt)-max(kurt)


  combine_x = np.array(list(X)+list(X))
  combine_y = np.array(list(Y1)+list(Y2))

  ################################
  # fit data with sinusoidal and quadratic curves; do statistical analysis of goodness-of-fit (calculate residuals sum of squares, RSS)

  try: params,covar = scipy.optimize.curve_fit(F,combine_x,combine_y,bounds=([0.,1.0,-PI,-0.2,-2.],[5.,2.0,PI,0.2,2.]))
  except: sys.stderr.write("sinusoidal fitting failed for %s\n" % orf); print orf,"fitting-failed"; continue

  stderrs = numpy.sqrt(numpy.diag(covar))

  Yhat = F(combine_x,*params)
  RSS_sin = sum([(y-yhat)**2 for y,yhat in zip(combine_y,Yhat)])

  corr = np.corrcoef(Yhat,combine_y)[0,1]

  freq = params[1]
  per = 55.0/freq
  if freq<1.00001 or freq>1.99999: continue
  ampl = params[0]
  if ampl<0.7: continue

  try: params2,covar2 = scipy.optimize.curve_fit(Q,combine_x,combine_y)
  except: sys.stderr.write("quadratic fitting failed for %s\n" % orf); print orf,"fitting-failed"; continue

  Yhat2 = Q(combine_x,*params2)
  RSS_quad = sum([(y-yhat)**2 for y,yhat in zip(combine_y,Yhat2)])

  ###############################
  # save values in results

  vals = [orf,gene,cog]
  vals += ["%0.1f" % m,"%0.1f" % s]

  for j in range(len(params)): vals.append("%0.3f" % params[j])
  vals.append("%0.2f" % RSS_sin)
  vals.append("%0.3f" % corr)
  vals.append("%0.1f" % per)

  vals += ["%0.6f" % x for x in params2]
  vals.append("%0.1f" % RSS_quad)
  vals.append("%0.3f" % (RSS_sin/RSS_quad))

  results.append(vals)

##########################################
# print out results

for res in results: print '\t'.join([str(x) for x in res])




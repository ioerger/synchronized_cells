import sys,numpy,scipy.stats
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

OD,OriTer,FhaA,skip = {},{},{},1
for line in open("OD.txt"):
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  tp,od,oriter,fha = float(w[0]),float(w[1]),float(w[2]),float(w[3])
  OD[tp] = od
  OriTer[tp] = oriter
  FhaA[tp] = fha

N = 10000
T = 100
#Shift = 20
#ChromDupInit = 20 # hours into cell cycle
#FHAstart = 0.8 # frac of cell cycle

def make_plot(chromosomes,FHA,biomass,popsize,fname,title):
  plt.plot(range(T),chromosomes,label="Ori/Ter")
  plt.plot(range(T),FHA,label="FHA %")
  plt.plot(range(T),biomass,label="biomass")
  plt.plot(range(T),popsize,label="population size")
  plt.ylim((0,4))
  plt.legend()
  plt.title(title)
  
  plt.xlim((0,55))
  plt.savefig(fname)

def sim(Shift=20,cycleMn=45,cycleSd=1,ChromDupInit=20,dupWinA=10,dupWinB=5,FHAstart=0.8,recovTime=5,growthRate=0.08,growthDelay=0.2,fname="temp.png"):
  chromosomes = [0]*T
  biomass = [0]*T
  popsize = [0]*T
  FHA = [0]*T

  recovTimes = numpy.random.normal(recovTime,recovTime,N) # should these be positive? mean=sd
  #dupWin = dupWinA+dupWinB*numpy.random.normal(0,1,N) # variable length duplication times
  dupWin = numpy.random.normal(dupWinA,dupWinB,N) # variable length duplication times
  FHAlevels = [1.0]*N
  cellSizes = [1.0]*N
  cellCounts = [1]*N
  cycleTimes = numpy.random.normal(cycleMn,cycleSd,N) 
  growthRates = growthRate*numpy.random.normal(1,0.3,N) 

  for t in range(T):
    for i in range(N):
      cellTime = t+Shift-recovTimes[i]
      frac = (cellTime%cycleTimes[i])/cycleTimes[i]
  
      dup = 1
      if t<recovTimes[i]: dup = 1
      elif t>recovTimes[i] and t<recovTimes[i]+dupWin[i]: dup = 2
      elif cellTime>cycleTimes[i]:
        a = cellTime%cycleTimes[i]
        b = ChromDupInit
        c = b+dupWin[i]
        if a>b and a<c: dup = 2
      chromosomes[t] += dup
  
      if frac>FHAstart: FHAlevels[i] = 1
      else: FHAlevels[i] *= 0.87 # decay
      FHA[t] += FHAlevels[i]
  
      if cellSizes[i]<2 and (t<20 or frac>growthDelay): 
      #if t>recovTime and cellSizes[i]<2 and frac>growthDelay: 
        cellSizes[i] = min(2,cellSizes[i]+growthRates[i])
      biomass[t] += cellCounts[i]*cellSizes[i]
  
      if t>recovTime and cellTime%cycleTimes[i]>(cellTime+1)%cycleTimes[i]:
        cellSizes[i] = 1
        cellCounts[i] *= 2
        #FHAlevels[i] = 0 # instead of turning them off,let them decay
      popsize[t] += cellCounts[i]
  
  biomass = [x/float(N) for x in biomass]
  chromosomes = [x/float(N) for x in chromosomes]
  FHA = [x/float(N) for x in FHA]
  popsize = [x/float(N) for x in popsize]

  title="Cycle=%s(%s), Shift=%s, Recov=%d, FHA=%s," % (cycleMn,cycleSd,Shift,recovTime,FHAstart)
  title += "\nChromDupInit=%s, DupWin=%s+~%s, GrowthRate=%s, Delay=%s" % (ChromDupInit,dupWinA,dupWinB,growthRate,growthDelay)
  #chromosomes = [x-0.25 for x in chromosomes]
  if fname!=None: make_plot(chromosomes,FHA,biomass,popsize,fname,title)  

  print title
  ODact,ODmod,OriTerAct,OriTerMod,FhaAct,FhaMod = [],[],[],[],[],[]
  for key,od in OD.items():
    tp = int(key)
    bm = biomass[tp] # model
    oriter = OriTer[key] # actual
    fhaa = FhaA[key] # actual
    #print key,tp,od,round(bm,2),oriter,chromosomes[tp],fhaa,round(FHA[tp],3)
    ODact.append(od); ODmod.append(bm)
    OriTerAct.append(oriter); OriTerMod.append(chromosomes[tp])
    FhaAct.append(fhaa); FhaMod.append(FHA[tp])
    
  cc_od = scipy.stats.pearsonr(ODact,ODmod)[0]
  cc_oriter = scipy.stats.pearsonr(OriTerAct,OriTerMod)[0]
  cc_fha = scipy.stats.pearsonr(FhaAct,FhaMod)[0]
  print "cc_od =",round(cc_od,4)
  print "cc_oriter =",round(cc_oriter,4)
  print "cc_fha =",round(cc_fha,4)
  return (cc_od,cc_oriter,cc_fha)

##################################################
# Phase1
#Optimize the params for Shift, cycleTime, growthRate, growthDelay, and FHAstart first, based on OD and FhaA

def optimize_cell_cycle():
 best_od,best_oriter,best_fha,best_prod1= 0,0,0,0
 # 5.5.3.3.3=675
 for cycleMn in numpy.linspace(30,50,5):
  for Shift in numpy.linspace(10,30,5):
    for growthRate in [0.07,0.08,0.09]:
      for growthDelay in [0.1,0.2,0.3]:
        for FHAstart in [0.7,0.8,0.9]:
          print
          cc_od,cc_oriter,cc_fha = sim(cycleMn=cycleMn,Shift=Shift,growthRate=growthRate,growthDelay=growthDelay,FHAstart=FHAstart,fname=None)
          if cc_od>best_od: print "*cc_od"; best_od = cc_od
          if cc_oriter>best_oriter: print "*cc_oriter"; best_oriter = cc_oriter
          if cc_fha>best_fha: print "*cc_fha"; best_fha = cc_fha
          prod1 = cc_od*cc_fha
          if prod1>best_prod1: print "*prod1"; best_prod1 = prod1

if len(sys.argv)>1 and sys.argv[1]=="phase1":
  optimize_cell_cycle()
  sys.exit(0)

##################################################
# Phase2
# Then optimize recovTime, ChromDupINit, and dupWin based on OriTer

def optimize_chrom_dup():
  best_od,best_oriter,best_fha,best_prod1= 0,0,0,0
  for recovTime in [2,3,4,5,6]:
    for ChromDupInit in numpy.linspace(12,24,4):
      for dupWinA in numpy.linspace(6,12,2):
        for dupWinB in [3,4,5,6,7]:
          print
          #cc_od,cc_oriter,cc_fha = sim(cycleMn=35,Shift=10,growthRate=0.08,growthDelay=0.3,FHAstart=0.8,recovTime=recovTime,ChromDupInit=ChromDupInit,dupWinA=dupWinA,dupWinB=dupWinB,fname=None)
          #cc_od,cc_oriter,cc_fha = sim(cycleMn=50,Shift=20,growthRate=0.08,growthDelay=0.3,FHAstart=0.8,recovTime=recovTime,ChromDupInit=ChromDupInit,dupWinA=dupWinA,dupWinB=dupWinB,fname=None)
          cc_od,cc_oriter,cc_fha = sim(cycleMn=35,Shift=10,growthRate=0.07,growthDelay=0.3,FHAstart=0.8,recovTime=recovTime,ChromDupInit=ChromDupInit,dupWinA=dupWinA,dupWinB=dupWinB,fname=None)
          if cc_od>best_od: print "*cc_od"; best_od = cc_od
          if cc_oriter>best_oriter: print "*cc_oriter"; best_oriter = cc_oriter
          if cc_fha>best_fha: print "*cc_fha"; best_fha = cc_fha
          prod1 = cc_od*cc_oriter
          if prod1>best_prod1: print "*prod1"; best_prod1 = prod1

if len(sys.argv)>1 and sys.argv[1]=="phase2":
  optimize_chrom_dup()
  sys.exit(0)

###################################

# users can try varying these parameters...

sim(cycleMn=35,Shift=10,growthRate=0.07,growthDelay=0.3,FHAstart=0.8,recovTime=6,ChromDupInit=20,dupWinA=12,dupWinB=7,fname="temp.png")

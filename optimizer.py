# Created by Karoline Kallis 2021-22
# Class to run optimized DTs to mimick a reference dose

import numpy as np
import random
import copy
from scipy import stats
from itertools import repeat
from scipy import optimize as opt
from Verification import*
import time
from line_profiler import LineProfiler
from skimage.metrics import structural_similarity as ssim
import gurobipy as gp
from gurobipy import GRB
#from doseKernel import*

class Optimizer():
    """class to optimize dose distribution to fit to prediction"""
    #dosekernel = []
    #doseDistribution = []
    loss = []
    DTdiff = []
    #optimizer = 'basin-hopping'
    optimizer = ''
    #optimizer = 'Powell'
    #optimizer = 'BFGS'
    #optimizer = 'Annealing'
    #optimizer = 'COBYLA'
    #costfunction = 'RRMSE'
    #costfunction = 'RelMSE'
    costfunction = ''
    prescribedDose = 0
    downsampleFactor = 4 # dose is downsmapled to every 4th voxel for the optimization # This paramter probably should be added to the .json input file
    niter=0
    verbose = False
    name=''

    def __init__(self, name, optimizer, costfunction, prescribedDose, verbose):
        self.name = name
        self.optimizer = optimizer
        self.costfunction = costfunction
        self.prescribedDose = prescribedDose
        self.verbose = verbose
        self.niter =0 

    # Mehtod to find the best input paramters for basin-hopping
    def gridSearch(self, doseKernel, prediction, doseDistribution, airKermaStrength, temperature,prescribedDose, step):

        print(doseDistribution[:,1])
        b =(0,60)
        bnds = tuple(repeat(b,len(doseDistribution[:,1])))
        print(bnds)
        start = time.time()
        idx = self.getPoints4Optimization(prediction, prescribedDose)
        
        
        kwargs = {'method':'COBYLA','bounds': bnds, 'tol' :1e-05, 'options' : {'rhobeg': 5.00, 'maxiter': 500, 'disp': False},'args':(prediction, doseDistribution, doseKernel, airKermaStrength,idx)}

        res = opt.basinhopping(self.func, doseDistribution[:,1], 
                 niter =20, T=temperature, stepsize =step, interval = 4, niter_success = 7, callback=None,
                 minimizer_kwargs=kwargs)
        print("global minimum: x0 = %.4f, x1 = %.4f, x2 = %.4f, f(x0) = %.4f" % (res.x[0],res.x[1],res.x[2], res.fun))

        end = time.time()
        print("Time ", end-start)
        doseDistribution[:,1] = res.x
        return res, end-start

    # Optimization function
    def minimize(self, doseKernel, prediction, doseDistribution, airKermaStrength, DT):

        loss=[]

        b =(0,90)
        bnds = tuple(repeat(b,len(doseDistribution[:,1])))

        start = time.time()
        idx = self.getPoints4Optimization(prediction, self.prescribedDose)
        
        #Constraints
        # DT>=0
        # DT <100
        cons = []
        for k in range(0,len(doseDistribution[:,1])):
            con1 = {'type':'ineq','fun':lambda x, a=-0.01,i=k : x[i] - a}
            con2 = {'type':'ineq','fun':lambda x, b=100, i=k : b - x[i]}
            cons.append(con1)
            cons.append(con2)

        # Different optimization functions! Cobyla mostly used
        if(self.optimizer =='basin-hopping'):
        #kwargs = {'method':'COBYLA','bounds': bnds, 'tol' :1e-05, 'options' : {'rhobeg': 5.00, 'maxiter': 500, 'disp': False},'args':(prediction, doseDistribution, doseKernel, airKermaStrength,idx, DP, dcmParam)}
            kwargs = {'method':'COBYLA','constraints': cons, 'tol' :1e-02, 'options' : {'rhobeg': 5.00, 'maxiter': 500, 'disp': False},'args':(prediction, doseDistribution, doseKernel, airKermaStrength,idx, DT, self.downsampleFactor, False)}

            res = opt.basinhopping(self.func, doseDistribution[:,1], 
                     niter =20, T=271, stepsize =1.0, interval = 4, niter_success = 5, callback=None,
                     minimizer_kwargs=kwargs)

        if(self.optimizer =='COBYLA'):
            res = opt.minimize(self.func, doseDistribution[:,1], args = (prediction, doseDistribution, doseKernel, airKermaStrength,idx, DT, self.downsampleFactor, False),
                               method = 'COBYLA',constraints=cons, tol =1e-3, callback = None, options = {'rhobeg': 5, 'maxiter': 1000, 'disp': False})
        if(self.optimizer =='Powell'):
            res = opt.minimize(self.func, doseDistribution[:,1], args = (prediction, doseDistribution, doseKernel, airKermaStrength,idx, DT, self.downsampleFactor, False), 
                       method = 'Powell',constraints=(), bounds = bnds, tol =1e-02, options = {'ftol': 1e-08,'xtol': 1e-08,'maxfev':500, 'maxiter': 500, 'disp': True})
        if(self.optimizer =='BFGS'):
            res = opt.minimize(self.func, doseDistribution[:,1], args = (prediction, doseDistribution, doseKernel, airKermaStrength, idx, DT, self.downsampleFactor, False), 
                    method = 'L-BFGS-B',constraints=(), bounds = bnds, tol =1e-03, options = {'ftol': 1e-08,'gtol': 1e-08, 'maxcor':20,'maxfun':1000, 'maxiter': 1000, 'disp': True, 'eps':2.5e-2, 'iprint':-1})
        if(self.optimizer =='Annealing'):
            res = opt.dual_annealing(self.func, bounds=bnds, args=(prediction, doseDistribution, doseKernel, airKermaStrength, idx, DT, self.downsampleFactor, False),
                             maxiter=1000, local_search_options={}, initial_temp=15.0, restart_temp_ratio=2e-05, visit=2.62, accept=-5.0, maxfun=10000000.0)
       
        end = time.time()
        print("Time ", end-start)



        if(self.optimizer != 'Gurobi'):
            res.x[res.x<0.1]=0
            doseDistribution[:,1] = res.x
            
        if(self.verbose):
           self.func(res.x, prediction, doseDistribution, doseKernel, airKermaStrength, idx, DT, self.downsampleFactor, True)
        return doseDistribution
    
    def cback(self, x):
        self.loss.append(self.func(x))


    #simulated annealing according to Lessart et al. ! 
    #Slow but runs! 
    def optimizeSimulatedAnnealing(self, doseKernel, prediction, doseDistribution, airKermaStrength, mask):
        run = True
        deltaE = 10000.0           
        doseDistributionPotential = copy.deepcopy(doseDistribution)
        k = 0
        alpha = 0.6 # empirical value based on Lessard et al. 
        T0 = 15 # Choose a useful starting parameter

        #while k!=500:
        while run:
            k+=1
            Tk = T0/k**alpha
            print("Tk ", Tk)
            print(k)
            print("Iteration ", k, "Time ", doseDistribution[:,1])
            timePotential = self.updateDoseDistribution(doseDistribution[:,1])
            doseDistributionPotential[:,1] = timePotential[:]
            sumDose = doseKernel.creatSumDose(doseDistribution, airKermaStrength, True)
            print("Iteration ", k, "TimePot ", doseDistributionPotential[:,1])
            sumDosePotential = doseKernel.creatSumDose(doseDistributionPotential, airKermaStrength, True)

            Ek = self.costFunction(prediction*mask, sumDose*mask) # think of a smart optimization function!! Let's start with absoulte differences between dose
            print("EK ", Ek)
            Ek1 = self.costFunction(prediction*mask, sumDosePotential*mask) # other option would be E/m finite dose points! 
            print("Ek1 ", Ek1)

            #compareDoseFiles(sumDose*mask, sumDosePotential*mask, 6.0)
            deltaE = Ek1 - Ek
            print("Delta ", deltaE)



            if(deltaE < 0):
               doseDistribution = copy.deepcopy(doseDistributionPotential)
            else:
                p = math.exp(-deltaE/Tk)
                print("p ", p)
                if(p>0.7):
                    doseDistribution = copy.deepcopy(doseDistributionPotential)

            if (deltaE==0 and k<1000):
                run= True
            elif(deltaE < 0.001 or k==5000):
                run=False

        print("k", k," Delta ", deltaE, " Time ", doseDistribution[:,1])
        return doseDistribution
    
    # To create a cirbular mask around the DPs and only use this reagion for optimization
    # not used
    def create_circular_mask(self,shape, spacing, center, radius1, radius2):
        Y, X, Z = np.ogrid[:shape[0], :shape[1], :shape[2]]
        dist_from_center = spacing*np.sqrt((X - center[0])**2 + (Y-center[1])**2+ (Z-center[2])**2)
        mask = (dist_from_center <= radius1) & (dist_from_center >= radius2)
        return mask

    # Defines the region for optimzation
    # toDO maket the bounds not hard coded
    def getPoints4Optimization(self, prediction, prescribedDose):
        idx = np.where(((prediction>=0.8*prescribedDose)&(prediction<=1.2*prescribedDose)),1,0)
        return idx


    def costFunction(self, prediction, calcDose, idx):

       ## Realtive MSE
        if(self.costfunction =='RelMSE'):
            Eges = 100*(np.square(np.subtract(calcDose[idx], prediction[idx])/prediction[idx])).mean()
        ## Realtive RMSE
        if(self.costfunction =='RRMSE'):
            Eges = 100*np.sqrt(((np.square(np.subtract(calcDose[idx], prediction[idx])/prediction[idx]))).mean())

        if(self.costfunction =='logMSE'):
            Eges = math.log((np.square(np.subtract(calcDose[idx], prediction[idx]))).mean())
        if(self.costfunction =='MSE'):
            Eges = (np.square(np.subtract(prediction[idx==1],calcDose[idx==1]))).mean()
        ## RMSE
        if(self.costfunction =='RMSE'):
            Eges = np.sqrt(np.square(np.subtract(calcDose[idx], prediction[idx])).mean())
        if(self.costfunction == 'absE'):
            Eges = (np.abs(np.subtract(calcDose[idx], prediction[idx]))).mean()

        if(self.costfunction == 'SSD'):
             Eges = np.square(np.subtract(calcDose[idx], prediction[idx])).sum()
        ## Structural-Similarity-Index
        if(self.costfunction =='SSIM'):
            Eges = - ssim(prediction[idx], calcDose[idx])

        ### Kullbeck-Leibler Divergence

        #p = np.random.normal(prediction[idx].mean(), prediction[idx].std(), len(prediction[idx]))
        #q = np.random.normal(calcDose[idx].mean(), calcDose[idx].std(), len(calcDose[idx]))

        #p = stats.gaussian_kde(prediction[idx])
        #q = stats.gaussian_kde(calcDose[idx])

        #valp = p.pdf(prediction[idx])
        #valq = p.pdf(calcDose[idx])
        ##print(valp * np.log2(valp / valq))
        ##plotDoseHistogram(prediction[idx], calcDose[idx], 6.0)
        ##plotDoseHistogram(valp, valq, 6.0)
        
        #Eges = np.abs(np.sum(np.where(valp != 0, valp * np.log2(valp /valq), 0)))

        ## Hellinger Distance
        
        #plotDoseHistogram(prediction[idx], calcDose[idx], 6.0)
        #print("Mean prediction ",prediction[idx].mean() )
        #print("Mean prediction ",prediction[idx].std() )

        #print("Mean calcDose ",calcDose[idx].mean())
        #print("Mean calcDose ",calcDose[idx].std())

        #p = np.random.normal(prediction[idx].mean() , prediction[idx].std() , len(prediction[idx]))
        #p_pdf=stats.multivariate_normal.pdf(p, mean= p.mean(), cov=np.cov(p))

        #q = np.random.normal(calcDose[idx].mean(), calcDose[idx].std(), len(calcDose[idx]))
        #q_pdf=stats.multivariate_normal.pdf(q, mean= q.mean(), cov=np.cov(q))

        #Eges = (1.0 / np.sqrt(2.0)) * np.sqrt(np.sum(np.sqrt(p_pdf) - np.sqrt(q_pdf))**2)

        # Bhattacharyya distance
        #Eges = 0.25 *np.log(0.25*((prediction[idx].std()**2/calcDose[idx].std()**2)+(calcDose[idx].std()**2/prediction[idx].std()**2)+2))+0.25*((prediction[idx].mean()-calcDose[idx].mean())**2/(prediction[idx].std()**2+calcDose[idx].std()**2))

        if(self.verbose):
            print('Eges',np.round(Eges,4))        

        #Eges = Eoar + Ehrctv
        return np.round(Eges,3)

    # function for simulated annealing
    # not used at the moment
    def updateDoseDistribution(self, time):
        #Forbit to have large jumps between DPs
        #Difference Ovoids...Tandem
        dwellUpdate = random.choices([-0.1, 0, 0.1], k=(time.shape[0]))
        print("Update ", dwellUpdate)

        dwellTimes = np.zeros(time.shape, dtype= time.dtype)

        for i in range(0, time.shape[0],1):
            dwellTimes[i] = time[i] + dwellUpdate[i]
            if(dwellTimes[i]<=0.1):
                dwellTimes[i]=0
        return dwellTimes

    # Function for optimzation
    # time = array of dwell times
    # prediction = 3D dose array
    # doseDistribution = list of dose arrays for each DP
    # idx = indicis of points used for the optimization
    # DT = Dicom DT -- as of now a DP.dcm for optimzation is used
    # factor = downsample factor which is set to 4
    # final = boolean to decide if it is the last iteration
    
    def func(self, time, prediction, doseDistribution, doseKernel, airKermaStrength, idx, DT, factor, final):
        self.niter+=1
        if(self.verbose):
            print("Niter: ", self.niter)
            print(time)

        doseDistribution[:,1] = time
        # create one dose file out of doseDistributions
        sumDose = doseKernel.createSumDose(doseDistribution, airKermaStrength, True, factor)
        # only used if intermediate results should be saved
        if(self.verbose):
            iterList=[50,100,125, 150, 200, 300, 400, 500]
            biter = self.niter in iterList
        # MESe between reference and current itreation are compared
        loss = self.costFunction(prediction[::factor,::factor,::factor], sumDose, idx[::factor, ::factor,::factor]) 
        self.loss.append(loss)
        DTDiff = []
        DTDiff = np.sqrt(((DT.reshape(1,DT.shape[0])-np.round(time,2))**2).mean())
        self.DTdiff.append(np.round(DTDiff,2))

        if(self.verbose):
           if(biter|final):
              compareDoseDistribution(prediction[::factor,::factor,::factor], sumDose, DT, time, loss, DTDiff, self.niter, self.name,  self.prescribedDose, factor)
           print ('Diff DT ', DTDiff, ' [s]')
        return loss

    def getDoseDistribution(self):
        return self.doseDistribution

    def setDoseDistribution(self, doseDistribution):
        self.doseDistribution = doseDistribution

    def getDoseKernel(self):
       return self.doseKernel


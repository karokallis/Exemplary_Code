# Created by Karoline Kallis 2021-22 
# Function to compare dose or image files!
import numpy as np

# Mean squared error
def MSE(origDose, calcDose,prescribedDose):
    diff = np.square(np.subtract(origDose[origDose < 1.5*prescribedDose], calcDose[origDose < 1.5*prescribedDose])).mean()
    return diff

# absolute differences
def absoluteDifference(origDose, calcDose,prescribedDose):
    diff = np.mean(np.absolute(origDose-calcDose))
    return diff

# realtive root mean squared error / realtive to origin Dose
def RRMSE(origDose, calcDose):
    Eges = 100*np.sqrt(((np.square(np.subtract(calcDose, origDose)/origDose))).mean())
    return Eges

# root mean squared error in a certain bound 
# lb = lower bound --1 =100%
# ub = upper bound
def isoRMSE(origDose, calcDose, lb,ub,  prescribedDose):
    idx = np.where((origDose>=lb*prescribedDose)&(origDose<=ub*prescribedDose))
    Eges = np.sqrt(((np.square(np.subtract(calcDose[idx], origDose[idx])))).mean())
    return Eges

# mean absolute error in a certain iso dose range
# lb = lower bound --1 =100%
# ub = upper bound
def isoMAE(origDose, calcDose, lb,ub,  prescribedDose):
    idx = np.where((origDose>=lb*prescribedDose)&(origDose<=ub*prescribedDose))
    Eges = np.mean(np.absolute(origDose[idx]-calcDose[idx]))
    return Eges

# root mean squared error
def RMSE(origDose, calcDose):
    Eges = np.sqrt(((np.square(np.subtract(calcDose, origDose)))).mean())
    return Eges

# mean error and standard deviation
def ME(origDose, calcDose):
    Diff = np.subtract(calcDose, origDose)
    Eges = Diff.mean()
    std = Diff.std() 
    return Eges, std

# mean error and standard deviation in certain isodose range
# lb = lower bound --1 =100%
# ub = upper bound
def isoME(origDose, calcDose, lb,ub,  prescribedDose):
    idx = np.where((origDose>=lb*prescribedDose)&(origDose<=ub*prescribedDose))
    Diff = np.subtract(calcDose[idx], origDose[idx])
    Eges = Diff.mean()
    std = Diff.std()
    return Eges, std

# Dice similarity coefficient
def dice(maskPred, maskQA):
    maskQA = np.asarray(maskQA).astype(np.bool)
    maskPred = np.asarray(maskPred).astype(np.bool)

    if maskPred.shape != maskQA.shape:
       raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")

    # Compute Dice coefficient
    intersection = np.logical_and(maskPred,maskQA)
    dice = 2.*intersection.sum()/(maskQA.sum() + maskPred.sum())
    return np.round(dice,2)

# false negative dice
def FND(maskPred, maskQA): # potential near misses
    maskPred = np.asarray(maskPred).astype(np.bool)
    maskQA = np.asarray(maskQA).astype(np.bool)

    if maskPred.shape != maskQA.shape:
       raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")

    intersection = np.logical_and(np.invert(maskPred), maskQA)
    FND = 2.*intersection.sum()/(maskPred.sum() + maskQA.sum())
    
    return np.round(FND,2)

# false positive dice
def FPD(maskPred, maskQA): # overtreatement 
    maskPred = np.asarray(maskPred).astype(np.bool)
    maskQA = np.asarray(maskQA).astype(np.bool)

    if maskPred.shape != maskQA.shape:
       raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")

    intersection = np.logical_and(maskPred, np.invert(maskQA))
    FPD = 2.*intersection.sum()/(maskPred.sum() + maskQA.sum())
    
    return np.round(FPD,2)

# valume differences
def volumetricDifference(volPre, volQA):
    VD = 100*(volPre - volQA)/volQA
    #print("VD" , VD)
    return np.round(VD,2)

# conformity index as defined in baltas et al
def COIN(PTVref, PTV, Vref):
    # c1 = PTVref / PTV --- relation of PTV covered by Dref and  PTV -- undercovered -- ideal c1 =1
    # c2 = PTVref / Vref --- relation of PTV covered by Dref and total volume of isodose -- overcovered -- ideal c2 =1 
    # COIN = c1*c2
    # Vref = volume of reference dose
    # PTVref = volume of the PTV covered by Dref

    c1 = PTVref/PTV
    c2 = PTVref/Vref
    COIN = c1*c2
    return np.round(COIN,2), np.round(c1,2), np.round(c2,2)



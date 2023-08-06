#Created by Karoline Kallis 2021-22

import numpy as np
#import dicom_numpy
import os
import pydicom as dcm
import math 
import copy
import time

import scipy.interpolate as interp
import skimage.transform
from scipy import ndimage
from plotResults import*
from Applicator import*
from Contour import*
from numpy.linalg import inv
from Verification import*

# Class to define the dose rate kernel
class doseKernel:
    """"""
    doseKernel = []
    doseKernelTime = 0 # dwell time of the kernel
    doseKernelArray = [] # dcm input saved as npy array
    doseKernelPosition = [] # the DPs defining the kernel 
    doseKernelAirKermaStrength = []
    kernelSpacing = [] # dose spacing the kernel is defined on
    kernelOrigin = [] # dose grid origin 
    kernelSize = []  # dose grid size of the kernel
    #pixelDimension = []
    prescribedDose = 0
    verbose = False
    def __init__(self, doseFile, positionFile, prescribedDose, verbose):

        self.prescribedDose = prescribedDose # has to be set before reading in the dose in case there should be a cut-off
        self.doseKernelArray, kernelParam = self.readDCM(doseFile, True, 0.25, False)
        self.kernelSpacing = kernelParam[1,:] 
        if(self.verbose):
            print("Spacing ", kernelParam[1,:])
        self.kernelOrigin = kernelParam[0,:]
        self.kernelSize = kernelParam[2,:]
        self.doseKernelPosition = self.getSourcePosition(positionFile) # 0 = Position / 1 = Orientation
        self.doseKernelPositionPixel = self.getMaxPixel(self.doseKernelArray)
        self.verbose = verbose

    def setkernelOrigin(self, origin):
        self.kernelOrigin= origin

    def setkernelSpacing(self, gridsize):
        self.kernelSpacing = (gridsize, gridsize, gridsize)

    # Key function of the program
    # Creates out of one dose rate kernel the reconstruction for the whole applicator
    # Dp = array of all dwell poitions in one applicator
    # DP Orientation = each of the DP is defined wiht a orientation since it is a approximated as line source
    # DwellTime = initial dwell time for each dwell position
    # dcmParam = array defining size, origin and spacing
    def createDoseDistribution(self, kernelDoseArray, DP, DPOrientation, DwellTimes, dcmParam):
        doseDistribution = np.empty([DP.shape[0],2], dtype = np.ndarray) 
        # n*2 array where 
        # n= amount of DPs
        # 1. element is the dose kernel moved to the position of the dwell position
        # 2. element is the dwell time of this DP
        for i in range(0, DP.shape[0],1):
            tic = time.time()
            #function to rotate and translate the dose rate kernel to specific DPs
            dose = self.moveDoseKernel(kernelDoseArray, DP[i,:], DPOrientation[i,:])
            # maps the dose array with a size of [400,400,400] to the size of the input
            # toDO: make function better to deal with large shifts when the patient position of the dicom differs a lot from the dose kernel
            doseDistribution[i,0] = self.mapDose2Actual(dose, dcmParam)
            doseDistribution[i,1] = DwellTimes[i]
            toc = time.time()
            if(self.verbose):
                print("i: " + str(i))
                print("Dwell Position: ", DP[i,:])
                print("Dwell Orientation: ", DPOrientation[i,:])
                print("Time Move/Map ", toc-tic)
        return doseDistribution
    

    # same function as createDoseDistribution but for multiprocessign to make it quciker
    def createDoseDistributionParallel(self, kernelDoseArray, DP, DPOrientation, DwellTime, dcmParam):
    
        dose = self.moveDoseKernel(kernelDoseArray, DP, DPOrientation)
        return np.array([self.mapDose2Actual(dose, dcmParam), DwellTime], dtype=object)


    # Function to update the dwell time in the dose distribution
    def updateDoseDistribution(self, doseDistribution, DwellTimes):
        updatedDose = copy.deepcopy(doseDistribution) # to make a copy and not just a pointer
        for i in range(0, doseDistribution.shape[0],1):
            updatedDose[i,1] = DwellTimes[i]

        return updatedDose
   
   # Function to add all single dose arrays to one sum dose to compare it tot he input dicom
    def createSumDose(self, doseDistribution, airKermaStrength,doseCutOff, factor):
        sumDose = np.zeros((int(np.ceil(doseDistribution[0,0].shape[0]/factor)), int(np.ceil(doseDistribution[0,0].shape[1]/factor)),int(np.ceil(doseDistribution[0,0].shape[2]/factor))), dtype = np.float) 
        for i in range(0, doseDistribution.shape[0],1):
            # scale dose kernel according to dwell time
            doseArray = self.scaleDoseKernel(doseDistribution[i,0], doseDistribution[i,1], airKermaStrength)
            sumDose = np.add(sumDose, doseArray[::factor, ::factor, ::factor])

        #toDo make dose cut of flexible and not only 200% of prescribed dose
        if(doseCutOff==True):
            sumDose[sumDose > 2.0*self.prescribedDose] =2.0*self.prescribedDose

        return sumDose

    # function to resample dose to a different grid and size
    # sumDose = dose array
    # dcmParam = parameters definingin size, sapcing and origin
    def mapDose2Actual(self, sumDose, dcmParam):
        #Get ratio of new and old dose Resolution
        OriginNew = self.world2pixelNearestNeighbor(dcmParam[0,:])[0]
        ratio = self.kernelSpacing/dcmParam[1]
        
        #Resample Dose Distribution
        outputArray = np.zeros((math.ceil(dcmParam[2,0]/ratio[0]),math.ceil(dcmParam[2,1]/ratio[1]),math.ceil(dcmParam[2,2]/ratio[2])), dtype = sumDose.dtype)
        # toDO: make the mapping smarter
        if((sumDose.shape[0]>=dcmParam[2,0]/ratio[0])&(sumDose.shape[1]>=dcmParam[2,1]/ratio[1])&(sumDose.shape[2]>=dcmParam[2,2]/ratio[2])):
            
            
            if(OriginNew[0]>=0):
                x = OriginNew[0]
                xend = outputArray.shape[0]+ OriginNew[0]
                xendOut = outputArray.shape[0]

                if(xend>sumDose.shape[0]):
                    xend = sumDose.shape[0]
                    xendOut = outputArray.shape[0] - (outputArray.shape[0]+ OriginNew[0] - sumDose.shape[0])

                xNew = 0
            else:
                x = 0
                xNew = (-1)*OriginNew[0]
                xend = outputArray.shape[0]+ OriginNew[0]
                xendOut = outputArray.shape[0]

                if(xend>sumDose.shape[0]):
                    xend = sumDose.shape[0]
                    xendOut = outputArray.shape[0] - (outputArray.shape[0]+ OriginNew[0] - sumDose.shape[0])


            if(OriginNew[1]>=0):
                y = OriginNew[1]
                yNew = 0
                yendOut = outputArray.shape[1]

                yend = outputArray.shape[1]+ OriginNew[1]
                if(yend>sumDose.shape[1]):
                   yend = sumDose.shape[1]                   
                   yendOut = outputArray.shape[1] - (outputArray.shape[1]+ OriginNew[1] - sumDose.shape[1])


            else:
                y =0
                yNew = (-1)*OriginNew[1]
                yend = outputArray.shape[1]+ OriginNew[1]
                yendOut = outputArray.shape[1]
                if(yend>sumDose.shape[1]):
                   yend = sumDose.shape[1]
                   yendOut = outputArray.shape[1] - (outputArray.shape[1]+ OriginNew[1] - sumDose.shape[1])


            if(OriginNew[2]>=0):
                z = OriginNew[2]
                zNew = 0
                zend = outputArray.shape[2]+ OriginNew[2]
                zendOut = outputArray.shape[2]

                if(zend>sumDose.shape[2]):
                   zend = sumDose.shape[2]
                   zendOut = outputArray.shape[2] - (outputArray.shape[2]+ OriginNew[2] - sumDose.shape[2])

            else:
                z =0
                zNew = (-1)*OriginNew[2]
                zend = outputArray.shape[2]+ OriginNew[2]
                zendOut = outputArray.shape[2]
                if(zend>sumDose.shape[2]):
                   zend = sumDose.shape[2]
                   zendOut = outputArray.shape[2] - (outputArray.shape[2]+ OriginNew[2] - sumDose.shape[2])

            outputArray[xNew:xendOut,yNew:yendOut,zNew:zendOut] = sumDose[x:xend,y:yend,z:zend]

        if((ratio[0]==0) & (ratio[1]==0) & (ratio[2]==0)):
            return outputArray
        else:
            resampledDose =  self.resample(outputArray, ratio)
            return resampledDose

    # Function to read in Dicom Dose File
    # doseFile = RD.dcm
    # kernel = boolean to decide if a kernel is read in or the reference dose -- no dose cut off if kernel
    # gridsize = desired girdsize -- will be igonred if resampel == false
    # resample = boolean to define gridsize to be used; if false the original grid size is used
    def readDCM(self, doseFile, kernel, gridsize, resample): 
        arrayParam = np.zeros([3,3], dtype = np.float)
        data = dcm.read_file(doseFile)
        doseArray = data.pixel_array * data.DoseGridScaling
        doseArray = np.array(doseArray)
        outputArray = np.transpose(np.moveaxis(doseArray,0,-1),(1,0,2))
        arrayParam[0,:] = data.ImagePositionPatient

        if(resample ==True): 
            spacing = np.array((data.PixelSpacing[0], data.PixelSpacing[1],data.GridFrameOffsetVector[1]))
            ratio = spacing/gridsize # newSpacing / oldSpacing
            outputArray = self.resample(outputArray, ratio)
            arrayParam[1,:]= gridsize
            arrayParam[2,:] = outputArray.shape
        else:
            arrayParam[1,:] = [data.PixelSpacing[0], data.PixelSpacing[1],data.GridFrameOffsetVector[1]] # Assuming that the dose grid is regtangular
            arrayParam[2,:] = outputArray.shape

        # High dose regions are not considered in optimization
        # Optimization in 0.8-1.5 isodose range -- dose in higher regions neglected
        # probably unnecssary since we only consider a certain isodose range for the optimisation
        if(kernel==False):
            outputArray[outputArray > 2.0*self.prescribedDose] = 2.0*self.prescribedDose
            #outputArray[outputArray > 1.5*self.prescribedDose] = -1
        return outputArray, arrayParam


    # Function to get position, orientation and time of dose kernel as stated in the dcm
    def getSourcePosition(self, positionFile):
        applicator = Applicator(positionFile)
        # applicator._init_(positionFile)
        DP = applicator.getApplicator() # source center
        DOrientation = applicator.getOrientation() # directional vector of the source
        time = applicator.getCumDwellTimes() # in seconds
        self.setdoseKernelAirKermaStrength(applicator.getAirKermaStrength())
        self.setdoseKernelTime(time[0])

        position = np.vstack((DP[0],DOrientation[0]))
        return position

    def getdoseKernelTime(self):
        if(self.doseKernelTime == 0):
            raise Exception("KK:: Dwell Time Kernel equal to zero and obivously not set!") 
        return self.doseKernelTime

    def setdoseKernelTime(self, time):
        self.doseKernelTime = time

    def setdoseKernelAirKermaStrength(self, airKermaStrength):
        self.doseKernelAirKermaStrength = airKermaStrength

    def getDoseKernelArray(self):
        if(self.doseKernelArray.size == 0):
            raise Exception("KK:: Dose kernel array not set yet!!!") 
        return self.doseKernelArray

    #Function not used but would theoretically move 2D arrays
    def moveDoseKernel2D(self,DP, DPOrientation,slice):
        transMatrix = self.calculateTranformationMatrix2D(DP, DPOrientation)
        warpedImage = ndimage.affine_transform(self.doseKernelArray[:,:,slice],transMatrix,mode='nearest', cval= 0.0, prefilter =True)
        return warpedImage

    #Method to calcualted the Transformation matrix moving the Kernel to a specific DP in 2D
    def calculateTranformationMatrix2D(self, DP, DPOrientation): # 2D rotation around axes =0
        # get Translation
        t = (int)((DP - self.doseKernelPosition[0])[0] / self.kernelSpacing)
       
        T = np.matrix(((1,0,t[0]),(0,1,t[1]),(0,0,1)))
        print(t)
        #getRotation
        beta = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation[0], 0)
        R = np.matrix(((math.cos(beta),-math.sin(beta),0),(math.sin(beta),math.cos(beta),0),(0,0,1)))

        #Rotated around center
        T_center = np.matrix(((1,0,self.doseKernelPositionPixel[0]),(0,1,self.doseKernelPositionPixel[1]),(0,0,1)))
        #print(T_center)
        T_inv_center = np.linalg.inv(T_center)
        tranMatrix = np.linalg.inv(T*T_center*R*T_inv_center)

        return tranMatrix
    #Function to move Kernel position to DPs
    # Transformation done in homohgenous coordinates
    def moveDP(self, DP, DPOrientation):
        transMatrix = self.calculateTranformationMatrixWorld(DP, DPOrientation) 
        DP_homo = np.reshape([self.doseKernelPosition[0][0],self.doseKernelPosition[0][1],self.doseKernelPosition[0][2],1],(4,1))
        DP_moved = np.dot(transMatrix,DP_homo)
        DP_result = np.reshape([DP_moved[0],DP_moved[1],DP_moved[2]],(1,3))
        return DP_result

    #
    def movePoint(self, DP, DPOrientation):
        transMatrix = self.calculateTranformationMatrixWorld(DP, DPOrientation)
        P_homo = np.reshape([DP[0],DP[1],DP[2],1],(4,1))
        P_moved = np.dot(transMatrix,P_homo)
        P_result = np.reshape([P_moved[0],P_moved[1],P_moved[2]],(1,3))
        return P_result

    # function to move the dose kernel to a particular DP
    # dosearray = 3D np.array of dose
    def moveDoseKernel(self,doseArray, DP, DPOrientation):
        transMatrix = self.calculateTranformationMatrix(DP, DPOrientation)
        doseFile = np.zeros(doseArray.shape, dtype = doseArray.dtype) 
        # backwardswarp!
        warpedImage = ndimage.affine_transform(doseArray,transMatrix, output_shape = doseFile.shape, output = doseFile, mode='constant', order = 2, prefilter =True)
        return warpedImage
     

    # def getMaxExtend(self, DP,center):
    #     x,y,z = 0,0,0
    #     for i in range(0, len(DP), 1):
    #         t = (DP[i] - center) / self.kernelSpacing
    #         if (abs(t[0])>x):
    #             x = t[0]
    #         if (abs(t[1])>y):
    #             y = t[1]
    #         if (abs(t[2])>z):
    #             z = t[2]

    #     return [math.ceil(x),math.ceil(y),math.ceil(z)]

    # Retunrs a Transformation matrix for a given translation in homogenous coordiantes
    def calculateTranformationMatrixTranslation(self, tranlsation): # 3D rotation in homogenous coordinates
        T = np.matrix(((1,0,0,tranlsation[0]),(0,1,0,tranlsation[1]),(0,0,1,tranlsation[2]),(0,0,0,1)))
        tranMatrix = np.linalg.inv(T)

        return tranMatrix

    # Returns a Transformation matrix for a given translation adn rotation in homogenous coordiantes
    def calculateTranformationMatrixTranslationRotation(self, center, tranlsation, rotation): # 3D rotation in homogenous coordinates
        T = np.matrix(((1,0,0,tranlsation[0]),(0,1,0,tranlsation[1]),(0,0,1,tranlsation[2]),(0,0,0,1)))
        
        Rx = np.matrix(((1,0,0,0),(0, math.cos(rotation[0]), -math.sin(rotation[0]),0),(0, math.sin(rotation[0]), math.cos(rotation[0]),0),(0,0,0,1)))
        Ry = np.matrix(((math.cos(rotation[1]),0,math.sin(rotation[1]),0),(0, 1, 0,0),(-math.sin(rotation[1]), 0, math.cos(rotation[1]),0),(0,0,0,1)))
        Rz = np.matrix(((math.cos(rotation[2]),-math.sin(rotation[2]),0,0),(math.sin(rotation[2]), math.cos(rotation[2]), 0,0),(0, 0, 1,0),(0,0,0,1)))
        R = Rz * Ry * Rx
        center = self.world2pixelSpace(center)[0]
        T_center = np.matrix(((1,0,0,center[0]),(0,1,0,center[1]),(0,0,1,center[2]),(0,0,0,1)))
        T_inv_center = np.linalg.inv(T_center)
        tranMatrix = np.linalg.inv(T*T_center*R*T_inv_center)
        # tranformation matrix has to be inverted because the pyhton function is backwards-warping
        return tranMatrix

    # calucalted the tranformation matrix to move the kernel to a certain DP in world coordiantes instead of pixel
    def calculateTranformationMatrixWorld(self, DP, DPOrientation): # 3D rotation in homogenous coordinates
        # get Translation
        t =  (DP-self.doseKernelPosition[0])
        T = np.matrix(((1,0,0,t[0]),(0,1,0,t[1]),(0,0,1,t[2]),(0,0,0,1)))

        alpha = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 0)
        beta = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 1)
        gamma = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 2)

        Rx = np.matrix(((1,0,0,0),(0, math.cos(alpha), -math.sin(alpha),0),(0, math.sin(alpha), math.cos(alpha),0),(0,0,0,1)))
        Ry = np.matrix(((math.cos(beta),0,math.sin(beta),0),(0, 1, 0,0),(-math.sin(beta), 0, math.cos(beta),0),(0,0,0,1)))
        Rz = np.matrix(((math.cos(gamma),-math.sin(gamma),0,0),(math.sin(gamma), math.cos(gamma), 0,0),(0, 0, 1,0),(0,0,0,1)))
        R = Rx * Ry * Rz
        
        #Rotated around center
        t_center = self.doseKernelPosition[0]
        T_center = np.matrix(((1,0,0,t_center[0]),(0,1,0,t_center[1]),(0,0,1,t_center[2]),(0,0,0,1)))
        T_inv_center = np.linalg.inv(T_center)
        tranMatrix =np.linalg.inv(T*T_center*R*T_inv_center)
        return tranMatrix

    def caluculateTranslation(self, DP):
        t = (DP - self.doseKernelPosition[0])/ self.kernelSpacing
        print("Translation: ", t)
        return t

    # This function is actually used for the optimization
    # The method calucalted transformation matrix to move the dose kernel to a specific DP
    def calculateTranformationMatrix(self,DP, DPOrientation): # 3D rotation in homogenous coordinates
        t = self.world2pixelDicomSpaceNearestNeighbor(DP, self.kernelOrigin, self.kernelSpacing)[0] - self.world2pixelNearestNeighbor(self.doseKernelPosition[0])[0]
        if(self.verbose):
            print("Translation", t , " Length ",np.linalg.norm(t) )

        T = np.matrix(((1,0,0,t[0]),(0,1,0,t[1]),(0,0,1,t[2]),(0,0,0,1)))
        #getRotation #TODO double check if rotation direction is correct!!!!
        alpha = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 0)
        beta = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 1)
        gamma = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 2)

        Rx = np.matrix(((1,0,0,0),(0, math.cos(alpha), -math.sin(alpha),0),(0, math.sin(alpha), math.cos(alpha),0),(0,0,0,1)))
        Ry = np.matrix(((math.cos(beta),0,math.sin(beta),0),(0, 1, 0,0),(-math.sin(beta), 0, math.cos(beta),0),(0,0,0,1)))
        Rz = np.matrix(((math.cos(gamma),-math.sin(gamma),0,0),(math.sin(gamma), math.cos(gamma), 0,0),(0, 0, 1,0),(0,0,0,1)))
        R = Rz * Ry * Rx
        
        #Rotated around center
        t_center = self.world2pixelNearestNeighbor(self.doseKernelPosition[0])[0]

        T_center = np.matrix(((1,0,0,t_center[0]),(0,1,0,t_center[1]),(0,0,1,t_center[2]),(0,0,0,1)))
        T_inv_center = np.linalg.inv(T_center)
        tranMatrix = np.linalg.inv(T*T_center*R *T_inv_center)
        return tranMatrix


    def calculateTranformationMatrixNonHomogenous(self,DP, DPOrientation): # 3D rotation in homogenous coordinates
        # get Translation
        print("Hello World")
        t = (DP - self.doseKernelPosition[0])/self.kernelSpacing
        #T = np.matrix(((1,0,0,t[0]),(0,1,0,t[1]),(0,0,1,t[2]),(0,0,0,1)))

        #getRotation #TODO double check if rotation direction is correct!!!!
        alpha = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 0)
        beta = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 1)
        gamma = self.calculateRotationAngle(self.doseKernelPosition[1], DPOrientation, 2)

        print("Alpha: ", np.degrees(alpha))
        print("Beta: ", np.degrees(beta))
        print("Gamma: ", np.degrees(gamma))
        #alpha = math.pi/4
        #beta = 0
        gamma = 0
        Rx = np.matrix(((1.0,0.0,0.0),(0.0, math.cos(alpha), -math.sin(alpha)),(0.0, math.sin(alpha), math.cos(alpha))))
        Ry = np.matrix(((math.cos(beta),0.0,math.sin(beta)),(0.0, 1.0, 0.0),(-math.sin(beta), 0.0, math.cos(beta))))
        Rz = np.matrix(((math.cos(gamma),-math.sin(gamma),0.0),(math.sin(gamma), math.cos(gamma), 0.0),(0.0, 0.0, 1.0)))
        R = Rz*Ry*Rx
        
        #Rotated around center
        t_center = self.world2pixelNearestNeighbor(self.doseKernelPosition[0])[0]
        #T_center = np.matrix(((1,0,0,t_center[0]),(0,1,0,t_center[1]),(0,0,1,t_center[2]),(0,0,0,1)))

        return t, R, t_center
     
    #Function to scale the dose array to particular dwell time
    # dose Array n*dimensional array 
    # time = new dwell time
    # ariKermaStrength of plan to be optimized
    def scaleDoseKernel(self, doseArray,  time, airKermaStrength):

        scalingFactor = (np.round(time,4)/np.round(self.doseKernelTime,4))*(airKermaStrength/self.doseKernelAirKermaStrength)
        scaledDoseKernelArray = doseArray*scalingFactor
        #if(self.verbose):
        #   print("ScalingFactor ", scalingFactor)
        return scaledDoseKernelArray

    def calculateAngle(self, v1, v2):
        uv1 = v1 / np.linalg.norm(v1)
        uv2 = v2 / np.linalg.norm(v2)

        dot_product = np.dot(uv1,np.transpose(uv2))
        angle = np.arccos(dot_product) * 180 / math.pi
        return angle
    
    def getValue(self, doseArray, x,y,z):
        return doseArray[x,y,z]

    def resample(self, doseArray, ratio):
        zoomed = ndimage.zoom(doseArray,ratio,output = doseArray.dtype, order=3)
        return zoomed


    def calculateRotationAngle(self, v1, v2, axes):
       if axes == 0:
           v1 = np.hstack((round(v1[1],2),round(v1[2],2)))
                #v1 =v1.reshape(1,2)[0]
           v2 = np.hstack((round(v2[1],2),round(v2[2],2)))
       elif axes == 1:
            v1 = np.hstack((round(v1[0],2),round(v1[2],2)))
            v2 = np.hstack((round(v2[0],2),round(v2[2],2)))
       elif axes == 2:
            v1 = np.hstack((round(v1[0],2),round(v1[1],2)))
            v2 = np.hstack((round(v2[0],2),round(v2[1],2)))
       else:
            raise Exception("KK:: Axes not defined!") 


       #if(np.linalg.norm(v1 - v2) != 0):
       if(np.linalg.norm(v1) != 0):
           v1_unit = v1 / np.linalg.norm(v1)
       else:
           v1_unit = v1

       if(np.linalg.norm(v2) != 0):
           v2_unit = v2 / np.linalg.norm(v2)
       else:
           v2_unit = v2

       dot_product = np.dot(v1_unit,v2_unit)
       angle = np.math.atan2(np.linalg.det([v1,v2]),np.dot(v1,v2))
       return angle

    # not used
    # function to move dose kernel into center of dcm
    # might actually solve the dose function mapping problem! 
    def center(self, kernelDoseArray, dcmParam):

        centerKernel = self.getCenter(self.kernelOrigin, self.kernelSize, self.kernelSpacing)
        centerDCM = self.getCenter(dcmParam[0,:], dcmParam[2,:],dcmParam[1,:])
        translation = (centerDCM-centerKernel)/ self.kernelSpacing
        T = np.matrix(((1,0,0,np.round(translation[0],2)),(0,1,0,np.round(translation[1],2)),(0,0,1,np.round(translation[2])),(0,0,0,1)))

        self.doseKernelPosition[0] = self.doseKernelPosition[0] + translation
        self.setkernelOrigin(self.kernelOrigin + translation )
        translatedDoseKernel = ndimage.affine_transform(kernelDoseArray,np.linalg.inv(T), mode='nearest', order = 2, prefilter =True)

        return translatedDoseKernel 

    def getCenter(self, origin, size, spacing):
        center = origin + spacing*size/2
        return center


    def pixel2worldSpace(self,pixel):
        world = np.zeros((3,1))
        world[0] = self.kernelSpacing[0] * pixel[0] + self.kernelOrigin[0]
        world[1] = self.kernelSpacing[1] * pixel[1] + self.kernelOrigin[1]
        world[2] = self.kernelSpacing[2] * pixel[2] + self.kernelOrigin[2]

        return np.transpose(world)

    def world2pixelDicomSpace(self,world, origin, spacing):
        pixel = np.zeros((3,1),dtype = int)
        pixel[0] = (int)(world[0] - origin[0]) / spacing[0]
        pixel[1] = (int)(world[1] - origin[1]) / spacing[1]
        pixel[2] = (int)(world[2] - origin[2]) / spacing[2]
        return np.transpose(pixel)[0]

    def world2pixelSpace(self,world):
        pixel = np.zeros((3,1),dtype = int)
        pixel[0] = (int)(world[0] - self.kernelOrigin[0]) / self.kernelSpacing[0]
        pixel[1] = (int)(world[1] - self.kernelOrigin[1]) / self.kernelSpacing[1]
        pixel[2] = (int)(world[2] - self.kernelOrigin[2]) / self.kernelSpacing[2]
        return np.transpose(pixel)

    # function to map world coordinate point into space defined by origin/spacing
    def world2pixelDicomSpaceNearestNeighbor(self,world, origin, spacing):
        pixel = np.zeros((3,1),dtype = int)
        pixel[0] = np.round((world[0] - origin[0]) / spacing[0],0)
        pixel[1] = np.round((world[1] - origin[1]) / spacing[1],0)
        pixel[2] = np.round((world[2] - origin[2]) / spacing[2],0)
        return np.transpose(pixel)

    # function to map world coordinate point into kernel space
    def world2pixelNearestNeighbor(self,world):
        pixel = np.zeros((3,1),dtype = int)
        pixel[0] = np.round((world[0] - self.kernelOrigin[0]) / self.kernelSpacing[0],0)
        pixel[1] = np.round((world[1] - self.kernelOrigin[1]) / self.kernelSpacing[1],0)
        pixel[2] = np.round((world[2] - self.kernelOrigin[2]) / self.kernelSpacing[2],0)
        return np.transpose(pixel)

    # def get_grid(self,x, y, z, homogenous=False):
    #     coords = np.indices((x, y, z),dtype = int).reshape(3, -1)
    #     return np.vstack((coords, np.ones(coords.shape[1]))) if homogenous else coords

    # def get_grid2D(self,x, y, homogenous=False):
    #     coords = np.indices((x, y),dtype = int).reshape(2, -1)
    #     return np.vstack((coords, np.ones(coords.shape[1]))) if homogenous else coords

    def getMaxPixel(self, array):
        index = np.unravel_index(np.argmax(array, axis=None), array.shape)
        # print(array[index])
        return index


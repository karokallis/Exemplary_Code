from tkinter import image_names
import numpy as np
import SimpleITK as sitk
import json
from scipy import ndimage
import tifffile as tf
import nibabel as nib

import os

from adaptTiff import*
# Function to tranform 3D world coordinates to pixel space
# world = input array in pixel space
def world2pixelArray(world, pixelOrigin, pixelSpacing):
    # [Px Py Py 1] = [Xx*Spacing Yx*Spacing 0 Originx *[i j 0 1]
    #                Xy*Spacing Yy*Spacing 0 Originy
    #                Xz*Spacing Yz*Spacing 0 Originz]

    output = np.zeros((len(world[:, 0]), 3), dtype=world.dtype)

    origin = np.full((len(world[:, 0]), 3), (pixelOrigin[0],
                     pixelOrigin[1], pixelOrigin[2]), dtype=float)

    output[:, 0] = np.round(
        (world[:, 0] - origin[:, 0]) / pixelSpacing[0], 0).astype(int)
    output[:, 1] = np.round(
        (world[:, 1] - origin[:, 1]) / pixelSpacing[1], 0).astype(int)
    output[:, 2] = np.round(
        (world[:, 2] - origin[:, 2]) / 1.0, 0).astype(int)

    return output.astype(int)

def pixel2Array3D(pixel, pixelOrigin, pixelSpacing):
    # [Px Py Py 1] = [Xx*Spacing Yx*Spacing 0 Originx *[i j 0 1]
    #                Xy*Spacing Yy*Spacing 0 Originy
    #                Xz*Spacing Yz*Spacing 0 Originz]

    output = np.zeros((len(pixel[:, 0]), 3), dtype=np.float16)
    origin = np.full((len(pixel[:, 0]), 3), (pixelOrigin[0],
                     pixelOrigin[1], pixelOrigin[2]), dtype=float)

    output[:, 0]  = pixel[:, 0]*pixelSpacing[0]+origin[:, 0]
    output[:, 1]  = pixel[:, 1]*pixelSpacing[1]+origin[:, 1]
    output[:, 2]  = pixel[:, 2]*pixelSpacing[2]+origin[:, 2]

    return output

def pixel2Array2D(pixel, slice, pixelOrigin, pixelSpacing):
    # [Px Py Py 1] = [Xx*Spacing Yx*Spacing 0 Originx *[i j 0 1]
    #                Xy*Spacing Yy*Spacing 0 Originy
    #                Xz*Spacing Yz*Spacing 0 Originz]

    output = np.zeros((len(pixel), 3), dtype=np.float16)
    origin = np.full((len(pixel[:, 0]), 3), (pixelOrigin[0],
                     pixelOrigin[1],pixelOrigin[2]), dtype=float)
    
    # repSlice = np.full((len(pixel[:, 0]), 1), slice, dtype=float)

    output[:, 0] = pixel[:,0]*pixelSpacing[0]+origin[:, 0]
    output[:, 1] = pixel[:,1]*pixelSpacing[1]+origin[:, 1]
    output[:, 2] = slice*pixelSpacing[2]+origin[:, 2]
    #output[:, 2] = repSlice[1,]

    return output

def rgb2gray(rgb):
    #return rgb[:,:,0]
    #return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140]).astype(np.uint8)
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

def resampleSITKImage(image, ratio, outputpath ,name, write):
    outputSpacing = [ratio, ratio, 1.0]* np.asarray(image.GetSpacing())

    inputSize = image.GetSize()
    inputSpacing = image.GetSpacing()

    outputSize = []
    outputSize.append(int(inputSize[0] * inputSpacing[0] / outputSpacing[0]+.5))
    outputSize.append(int(inputSize[1] * inputSpacing[1] / outputSpacing[1]+.5))
    outputSize.append(int(inputSize[2] * inputSpacing[2] / outputSpacing[2]+.5))
    # resampler = sitk.ResampleImageFilter()
    # resampler.SetSize(outputSize)
    # resampler.SetOutputSpacing(outputSpacing)
    # resampler.SetOutputOrigin(image.GetOrigin())
    # resampler.SetOutputDirection(image.GetDirection())
    # resampler.SetInterpolator(sitk.sitkLinear)
    # resampler.SetDefaultPixelValue(0)
    # image = resampler.Execute(image)

    
    image = sitk.Resample(image, outputSize, sitk.ScaleTransform(3, [1/ratio ,1/ratio, 1]),
                sitk.sitkNearestNeighbor, image.GetOrigin(), outputSpacing, image.GetDirection(),0.0, image.GetPixelID())

    if(write):
        #sitk.WriteImage(image, 'testl.prostate.nii.gz')
        sitk.WriteImage(image, os.path.join(outputpath, name))
    return image


    # Adated by Karoline Kallis 7/5/2022
    # Function to find slice correspondence of MRI to Histopatology
    # Input MRI prostate imC, reference Histopatholgy slice
    # should return an array of correspondences
    # similarity maybe established with dice and hellinger or only hellinger?
def findCorrespondence(self, moving2D, fixed3D, fixedMask, ref, refMask, refIdx):

    dice = 0.0
    idx = -1
    overlapMeasures = sitk.LabelOverlapMeasuresImageFilter()

    for ifix in range(fixed3D.GetSize()[2]-(moving2D.refSize[2]-1)):
        moving2D.registerToConstrait(
            fixed3D[:, :, ifix], self.refWoContraints, self.mskRefWoContraints, ref, refMask, ifix, True)
        # mask = moving2D.loadMask(refIdx)
        mask = moving2D.setTransformedMask(refMask, 0, 0)[:, :, refIdx]

        mask = sitk.Cast(sitk.Resample(mask, fixedMask[:, :, ifix], sitk.Transform(),
                                       sitk.sitkNearestNeighbor, 0.0, fixedMask[:, :, ifix].GetPixelID()) > 0, sitk.sitkUInt16)

        overlapMeasures.Execute(fixedMask[:, :, ifix], mask)

        sim = overlapMeasures.GetDiceCoefficient()

        if (sim > dice):
            print("--------")
            print("Idx optimized " + str(ifix))
            print("Dice optimized " + str(sim))
            dice = sim
            idx = ifix
    return idx

# Function to find correspondance between MRI an dhistopatholgy
# return a restricted image of histpathology
def restrictConstraints(image,imgMask, histo):
    refidx = 0 # Use first histo slice for now

    moving_image = histo[refidx].loadRgbImage()
    moving_image = histo[refidx].getGrayFromRGB(moving_image)
    mask = histo[refidx].loadMask(refidx)

    reg = sitk.ImageRegistrationMethod()

    # initial_transform = sitk.CenteredTransformInitializer(sitk.Cast(image[:,:,3],sitk.sitkFloat64), 
    #                                                   moving_image, 
    #                                                   sitk.Euler2DTransform(), 
    #                                                   sitk.CenteredTransformInitializerFilter.GEOMETRY) # the registration problem is purely a shift
    # Similarity metric settings.
    #reg.SetMetricAsMeanSquares()
    #reg.SetMetricAsCorrelation()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50) # MSE as similarity metric
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(1.00)
    reg.SetInterpolator(sitk.sitkLinear)

    reg.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=500, convergenceMinimumValue=1e-5, convergenceWindowSize=10)

    dice = 0.0
    idx = 0
    overlapMeasures = sitk.LabelOverlapMeasuresImageFilter()

    for fidx in range(image.GetSize()[2]-len(histo)):
        #final_transform = reg.Execute(imgMask[:,:,fidx], mask[:,:,refidx])
        initial_transform = sitk.CenteredTransformInitializer(sitk.Cast(image[:,:,fidx],sitk.sitkFloat64), 
                                                      moving_image,
                                                      sitk.Euler2DTransform(),
                                                      sitk.CenteredTransformInitializerFilter.GEOMETRY)
        # initial_transform = sitk.CenteredTransformInitializer(sitk.Cast(image[:,:,fidx],sitk.sitkFloat64), 
        #                                               moving_image,
        #                                               sitk.AffineTransform(2), 
        #                                               sitk.CenteredTransformInitializerFilter.GEOMETRY)
    
        reg.SetInitialTransform(initial_transform)

        # final_transform = reg.Execute(sitk.Cast(imgMask[:,:,fidx],sitk.sitkFloat32), sitk.Cast(mask,sitk.sitkFloat32))
        final_transform = reg.Execute(sitk.Cast(image[:,:,fidx],sitk.sitkFloat64), moving_image)
        movedHistoMask = sitk.Resample(mask, imgMask[:,:,fidx], final_transform, sitk.sitkNearestNeighbor, 0.0, mask.GetPixelID())
        movedHisto = sitk.Resample(moving_image, image[:,:,fidx], final_transform, sitk.sitkLinear, 0.0, moving_image.GetPixelID())
        
        # sitk.WriteImage(image[:,:,fidx], 'fixed.nii.gz')
        # sitk.WriteImage(imgMask[:,:,fidx], 'fixedMask.nii.gz')
        # sitk.WriteImage(movedHistoMask, 'movedHistoMask.nii.gz')
        # sitk.WriteImage(movedHisto, 'movedHisto.nii.gz')

        overlapMeasures.Execute(sitk.Cast(imgMask[:,:,fidx],sitk.sitkUInt8), movedHistoMask)
        sim = overlapMeasures.GetDiceCoefficient()

        if (sim > dice):
            print("--------")
            print("Idx optimized " + str(fidx))
            print("Dice optimized " + str(sim))
            dice = sim
            idx = fidx
            sitk.WriteImage(movedHistoMask, 'movedHistoMask.nii.gz')
            sitk.WriteImage(movedHisto, 'movedHisto.nii.gz')
    
    sitk.WriteImage(image[:,:,idx], 'fixed.nii.gz')
    sitk.WriteImage(imgMask[:,:,idx], 'fixedMask.nii.gz')

    return image[:,:,idx:idx+len(histo)]


def pickCorrespondence(img, corres):
    imgA = sitk.GetArrayFromImage(img)
    #img = np.rollaxis(img, 2, 0)

    imstack = np.zeros([len(corres), imgA.shape[1], imgA.shape[2]], dtype=imgA.dtype)

    for i, idx in enumerate(corres):
        imstack[i,:,:]=imgA[idx,:,:]
    imstack = sitk.GetImageFromArray(imstack)
    imstack.SetSpacing(img.GetSpacing())
    imstack.SetOrigin(img.GetOrigin())
    imstack.SetDirection(img.GetDirection())

    return imstack



def centerHistopathology(patID, path, resolution, forceFlag):
    with open(path+'Histopathology/'+patID+'.json') as jF:
        data = json.load(jF)
        for i in data.keys():
            dict = data[i]
            img = readTiff(dict['filename'])
            prostate = readTiff(dict['regions']['region00']['filename'])
            tumor = readTiff(dict['regions']['region01']['filename'])

            center = ndimage.center_of_mass(prostate)
            translation = [int(prostate.shape[0]/2-center[0]), int(prostate.shape[1]/2-center[1])]
            print("Translation: ", np.linalg.norm(translation))

            if(np.linalg.norm(translation)>10):
                translationMatrix = np.linalg.inv(np.matrix(((1,0,translation[0]),(0,1,translation[1]),(0,0,1))))

                prostateWarped = ndimage.affine_transform(prostate,translationMatrix, output_shape = prostate.shape, order = 0 ,mode='constant')
                tumorWarped = ndimage.affine_transform(tumor,translationMatrix, output_shape = prostate.shape, mode='constant')
                imageWarped = np.zeros(img.shape, dtype = img.dtype)
                for j in range(0, img.shape[2]):
                    imageWarped[:,:,j] = ndimage.affine_transform(img[:,:,j],translationMatrix, output_shape = prostate.shape, mode='constant')


                writeTiff(dict['regions']['region00']['filename'], '', prostateWarped, resolution)
                writeTiff(dict['regions']['region01']['filename'], '', tumorWarped, resolution)
                if(not forceFlag):
                    writeTiff(dict['filename'], '', imageWarped, resolution)
            else:
                writeTiff(dict['regions']['region00']['filename'], '', prostate, resolution)
                writeTiff(dict['regions']['region01']['filename'], '', tumor, resolution)
                if(not forceFlag):
                    writeTiff(dict['filename'], '', img, resolution)
                print(patID+' :Images already centered!!!')


def preProcessingRSI(path, dataset):
    image = nib.load(f"{dataset}") 
    data = image.get_fdata()


    if('rsi' in dataset.lower()):
        #img = sitk.GetImageFromArray(100*np.transpose(data[:,:,:,0]))
        img = sitk.GetImageFromArray(np.rollaxis(100*data[:,:,:,0].transpose(1,0,2), 2, 0))
        img.SetSpacing(image.header['delta'].astype(np.double))
        img.SetDirection([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,image.header['Mdc'].astype(np.double).flatten()[-1]]) #
        #img.SetDirection(image.header['Mdc'].astype(np.double).flatten())
        img.SetOrigin((-1.0, -1.0,1.0)*image.affine[0:3,3].astype(np.double))# Not sure if that is correct

        sitk.WriteImage(img, os.path.join(path,"rsi.nii.gz"))
        #nib.save(100*img[:,:,:,0], "rsi.nii.gz") # Save C1 map! # I believe the randomly multiply by 100 to make the values larger
    elif('contour' in dataset.lower()):
        img = sitk.GetImageFromArray(np.rollaxis(data.transpose(1,0,2), 2, 0))
        img.SetSpacing(image.header['delta'].astype(np.double))
        img.SetDirection([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,image.header['Mdc'].astype(np.double).flatten()[-1]]) #
        #img.SetDirection([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]) # Yeah! sitk!!!
        #img.SetDirection(image.header['Mdc'].astype(np.double).flatten())
        img.SetOrigin((-1.0, -1.0,1.0)*image.affine[0:3,3].astype(np.double))# Not sure if that is correct
        #resampleSITKImage(img, ratio, path ,"rsi.nii.gz", True)
        sitk.WriteImage(img, os.path.join(path,"dwi_prostate_contour.nii.gz"))

        #nib.save(img, "contour.nii.gz") #toDO make the naming smarter or at least more accurate!

# Function write registration results to dcm to be imported to MIM
# Created by Karoline Kallis 01/19/2023
import os
from datetime import date
import pydicom
from pydicom.dataset import Dataset, FileDataset
import numpy as np
#import math
from scipy.spatial import ConvexHull, Delaunay
from ypstruct import struct
import SimpleITK as sitk
from helperFunction import*
from skimage.measure import approximate_polygon
import alphashape 

# from rt_utils import RTStructBuilder

def writeDicomSeries(img, contour_list, contour_name_list, patID, name, outdcmpath,dicomMRIFolder, resolution, origin, grayScale):
    # You should change 'test' to your preferred folder.
    MYDIR = (outdcmpath)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)
    else:
        print(MYDIR, "folder already exists.")

    duid_list = []
    vol = struct()
        # TODO: Find smarter solution to import values out of .json
    vol.dimr = img.GetSize()[1]
    vol.dimc = img.GetSize()[0]
    vol.vx = resolution[0]  # in mm
    vol.vy = resolution[1]
    vol.vz = resolution[2]
    vol.origin = img.GetOrigin()
    vol.Mvx2lph = np.multiply(np.eye(4), np.array(
            [[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1]]))

    ds =setMetaData(patID, img.GetSize()[2], vol, name, grayScale)
    img_array = sitk.GetArrayFromImage(img)
    # if(img_array.ndim>3 and grayScale):
    #     img_array = np.moveaxis(img_array, [0, 1], [-2, -3])

    for i in range(0,img.GetSize()[2]): 
        if(img_array.ndim>3 and grayScale):
            tmp_img = rgb2gray(img_array[i,:,:,:])
            #tmp_img = np.transpose(tmp_img)
        elif(img_array.ndim>3 and not grayScale):
            tmp_img = img_array[i,:,:,:].reshape((img_array.shape[3],img_array.shape[2],img_array.shape[1]))
        else:
            tmp_img = img_array[i,:,:]
        duid_list = writeImage(tmp_img,outdcmpath, name, ds, i, vol, duid_list) 
    # writeRTStructReference(dcmFile, contour_list, outdcmpath, ds)
    writeDicomREG(outdcmpath, patID, dicomMRIFolder)
    writeRTSTruct(outdcmpath, patID, contour_list, contour_name_list, ds, duid_list, resolution, origin)
    # writeRTSTructReference(outdcmpath, patID, contour_list, contour_name_list, ds, duid_list)
    return 1


#def writeRTStructReference(dcmFile, contour_list, outdcmpath, ds):
    data = pydicom.dcmread(dcmFile)
    # toDO: Update UIDs to new dicom otherwise that will not be working
    # Decide if you want to just add ROI or delete and create new?
    data.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID = ds.FrameOfReferenceUID
    for i, contour in enumerate(data.StructureSetROISequence):
        contourName = contour.ROIName
        if 'prostate' in contourName.lower():
            data.StructureSetROISequence[i].ROIName = 'HP_Prostate_deformed'
            prostate_img = sitk.GetArrayFromImage(contour_list[0])
            prostate_array = np.moveaxis(prostate_img, [0], [2])

            for slice in range(0, prostate_array.shape[2]):
                print("Slice: "  + slice)
                contourData2D = getListBoundary(prostate_array[:,:,slice])
                worldcontourData = pixel2Array2D(contourData2D, slice, contour_list[0].GetOrigin(), contour_list[0].GetSpacing())
                #worldcontourData_Sorted = sortContourPoints(worldcontourData)              
#  # contourData needs to be reformated!
                data.ROIContourSequence[i].ContourSequence[slice].ContourData = np.ravel(worldcontourData).tolist()
            # put code for tumor here
        else:
            data.StructureSetROISequence[i].ROIName = ''
            data.ROIContourSequence[i].ContourSequence[slice].ContourData = ''
    
    print('RS_hist_moved successfully saved!')

    pydicom.filewriter.dcmwrite(
        outdcmpath+'RS_hist_moved.dcm', data, write_like_original=False)

def writeRTSTruct(outdcmpath, patID, contour_list, contour_name_list, dsImage, duid_list,spacing, origin,):
    #Adapt dicom tags to write a structure set
    #ds.TransferSyntaxUID = '1.2.840.10008.1.2';
    # data = pydicom.dcmread(dcmFile)
    duid = pydicom.uid.generate_uid()

    file_meta = Dataset()
    #file_meta.MediaStorageSOPClassUID = '1.2.246.352.70.2.1.160.3'
    file_meta.FileMetaInformationGroupLength = 202
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    file_meta.MediaStorageSOPInstanceUID = duid
    file_meta.ImplementationClassUID = '2.16.840.1.114362.1'

    #file_meta.MediaStorageSOPInstanceUID ???
    # file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
    file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds = FileDataset(patID, {}, file_meta=file_meta, preamble=b"\0"*128)
    #ds.ImageType = "DERIVATED\SECONDARY"
    ds.SpecificCharacterSet = 'ISO_IR 192'
    ds.InstanceCreationDate = date.today()
    ds.InstanceCreationTime = '000000'
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3' # RTStructureSetStorage
    #ds.SOPInstaneUID = '2.16.840.1.114362.1.12156802.25748924467.629960731.346.26'
    duid = pydicom.uid.generate_uid()
    ds.SOPInstanceUID = duid
    ds.StudyDate = date.today()
    ds.StudyTime = '000000'
    ds.SeriesDate = date.today()
    ds.SeriesTime = '000000'

    ds.AccessionNumber = ''
    ds.Modality = 'RTSTRUCT'
    #ds.Manufacture = 'MIM' # probably doens't matter
    #ds.SecondaryCaptureDeviceManufacturer = 'Python 3.x'
    ds.StudyDescription = 'PATH VOLUMES Registered'
    ds.SeriesDescription = 'Moved Contour'

    ds.PatientID = patID
    ds.PatientName = 'ProsPath'
    ds.PatientBirthDate = '20010101'
    ds.PatientBirthTime = '000000'
    ds.PatientSex = 'M'
    ds.PatientAge = ''
    
    ds.StudyInstanceUID = dsImage.StudyInstanceUID
    ds.SeriesInstanceUID = duid
    
    ds.StudyID = dsImage.StudyID
    ds.SeriesNumber = 1
    ds.StructureSetLabel = 'RTstruct_deformed'
    ds.StructureSetName = ''
    ds.StructureSetDate = date.today()
    ds.StructureSetTime = '000000'

    referenceSequence = Dataset()
    referenceSequence.FrameOfReferenceUID = dsImage.FrameOfReferenceUID

    RTReferenced = Dataset()
    RTReferenced.ReferencedSOPClassUID = dsImage.SOPClassUID 
    RTReferenced.ReferencedSOPInstanceUID=dsImage.StudyInstanceUID
    
    RTSequence = Dataset()
    RTSequence.SeriesInstanceUID = dsImage.SeriesInstanceUID

    RTSequence.ContourImageSequence = pydicom.sequence.Sequence() ### Fill something in here
    for item in duid_list:   
        contourImageSequence = Dataset()
        contourImageSequence.ReferencedSOPClassUID = dsImage.SOPClassUID 
        contourImageSequence.ReferencedSOPInstanceUID = item
        RTSequence.ContourImageSequence.append(contourImageSequence)


    RTReferenced.RTReferencedSeriesSequence = pydicom.sequence.Sequence([RTSequence])
    referenceSequence.RTReferencedStudySequence = pydicom.sequence.Sequence([RTReferenced])
    ds.ReferencedFrameOfReferenceSequence = pydicom.sequence.Sequence([referenceSequence])

    ds.StructureSetROISequence = pydicom.sequence.Sequence()

    for i in range(len(contour_list)):
        #print(i)
        ROI = Dataset()
        ROI.ROINumber = i+1
        ROI.ReferencedFrameOfReferenceUID = dsImage.FrameOfReferenceUID
        ROI.ROIName = contour_name_list[i]
        ROI.ROIDescription = 'moved'
        ROI.ROIGenerationAlgorithm = 'MANUAL'
        ds.StructureSetROISequence.append(ROI) 
    
    ds.ROIContourSequence = pydicom.sequence.Sequence()

    for i, item in enumerate(contour_list):
        ROIContourSequence = Dataset()
        ROIContourSequence.ROIDisplayColor = np.random.randint(1, [255,255,255]).tolist()
        ROIContourSequence.ReferencedROINumber = i+1
        ROIContourSequence.ContourSequence =pydicom.sequence.Sequence()

        # img = sitk.GetArrayFromImage(item)
        # img_array = np.moveaxis(img, [0], [2])
        
        for slice in range(item.shape[2]):
            ContourSequence = Dataset()
            ContourSequence.ContourGeometricType = 'CLOSED_PLANAR'

            if np.sum(item[:,:,slice]!=0):
                contourData2D = getListBoundary(item[:,:,slice])
                worldcontourData = pixel2Array2D(contourData2D, slice, origin, spacing)

                ContourSequence.NumberOfContourPoints = len(worldcontourData) ### might not be working!
                ContourSequence.ContourData = np.ravel(worldcontourData).tolist()
                #ContourSequence.NumberOfContourPoints = data.ROIContourSequence[i].ContourSequence[slice].NumberOfContourPoints
                #ContourSequence.ContourData = data.ROIContourSequence[i].ContourSequence[slice].ContourData # Not correct at the moment slice is not 0-length
                ContourImageSequence = Dataset()
                ContourImageSequence.ReferencedSOPClassUID = dsImage.SOPClassUID
                ContourImageSequence.ReferencedSOPInstanceUID = duid_list[slice]
                ContourSequence.ContourImageSequence = pydicom.sequence.Sequence([ContourImageSequence])
                ROIContourSequence.ContourSequence.append(ContourSequence)
            else:
                continue

        ds.ROIContourSequence.append(ROIContourSequence)

    # ds.RTROIObservationsSequence = pydicom.sequence.Sequence()
    # for i in range(len(contour_list)):
    #     ObservationSequence = Dataset()
    #     ObservationSequence.ObservationNumber = i
    #     ObservationSequence.ReferencedROINumber = i
    #     ObservationSequence.ROIObservationDescription = "Type:Soft,Range:*/*,Fill:0,Opacity:0.5,Thickness:1,LineThickness:2,read-only:false"
    #     ObservationSequence.ROIInterpreter = ''
    #     ObservationSequence.RTROIInterpretedType = ''
    #     ds.ROIObservationsSequence.append(ObservationSequence)

    ds.ApprovalStatus = 'UNAPPROVED'

    pydicom.filewriter.dcmwrite(
        os.path.join(outdcmpath,'RS_hist_moved.dcm'), ds, write_like_original=False)
    
    print('RS_hist_moved successfully saved!')



def writeImage(img, outdcmpath, name, ds, instance, vol, duid_list):
        # setmetadata information -- needs to be done smarter! Maybe using MRI header?
        # The values are very randoms
        # # Metadata
        # if(grayScale):
        #     img = rgb2gray(img)
    if img.dtype != np.uint16:
        print("Datatype: " + str(img.dtype))
        img = img.astype(np.uint16)
    else:
        if img.dtype != np.uint8:
            print("Datatype: " + str(img.dtype))
            img = img.astype(np.uint8)
        #img = img.convert('RGB')
    duid = pydicom.uid.generate_uid()
    ds.SOPInstanceUID = duid
    duid_list.append(duid)
    ds.file_meta.MediaStorageSOPInstanceUID = duid
    ds.InstanceNumber = instance
    ds.SliceLocation = instance+1
    #ds.Modality = 'CT'
        #ds.ImagePositionPatient = [0, 0, instance]
    ds.ImagePositionPatient = [vol.origin[0],vol.origin[1],vol.origin[2]+instance*vol.vz]
    pydicom.dataset.validate_file_meta(ds.file_meta, enforce_standard=True)

    ds.PixelData = img.tobytes()
    saveName = name+'_'+str(instance)+'.dcm'
    print(saveName)
    pydicom.filewriter.dcmwrite(
        os.path.join(outdcmpath,saveName), ds, write_like_original=False)
    return duid_list

def setMetaData(patID, nimages, vol, name, grayScale):
        # setmetadata information -- needs to be done smarter! Maybe using MRI header?
        # The values are very random
        # Metadata

    file_meta = Dataset()
    if ~grayScale:
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.6.3'
    else:
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
    file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    ds = FileDataset(patID, {}, file_meta=file_meta, preamble=b"\0"*128)
    ds.ImageType = "DERIVATED\SECONDARY"
    ds.SecondaryCaptureDeviceManufacturer = 'Python 3.x'
    ds.PatientID = patID
    ds.PatientName = 'ProsPath'
    ds.PatientBirthDate = '20010101'
    ds.PatientBirthTime = '000000'
    ds.PatientSex = 'M'
    ds.StudyDate = date.today()
    ds.StudyDescription = 'PATH VOLUMES Registered'
    ds.StudyInstanceUID = pydicom.uid.generate_uid()
    ds.StudyTime = '000000'
    ds.FrameOfReferenceUID = pydicom.uid.generate_uid()
    ds.ImagesInAcquisition = nimages
    ds.PixelSpacing = [vol.vx, vol.vy]
    ds.Columns = vol.dimc
    ds.Rows = vol.dimr
    ds.SliceThickness = str(vol.vz)
    ds.AcquisitionMatrix = [vol.dimc, vol.dimr, 0, 0]
    #ds.ImagePositionPatient = [vol.origin[0] ,vol.origin[1],vol.origin[2]]

    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    if(grayScale):
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
        ds.Modality = 'CT'
        ds.StudyID = '7001'
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.PixelRepresentation = 0
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
    else:
        ds.Modality = 'US'
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.6.3' 
        ds.StudyID = '7002'
        ds.NumberOfFrames = '1'
        ds.SamplesPerPixel = 3
        ds.PhotometricInterpretation = "RGB"
        ds.PixelRepresentation = 0
        ds.PlanarConfiguration = 0
        ds.BitsAllocated = 8
        ds.BitsStored = 8
        ds.HighBit = 7
        #ds.PurposeOfReferenceCodeSequence ='121322' #"Source of Image Processing Operation 0040,A170
        #ds.DerivationCodeSequence='113091' #Spatially-related frames extracted from the volume

    ds.SeriesDescription = name
    ds.RescaleIntercept = '0'
    ds.RescaleSlope = '1'
    ds.SeriesNumber = 7000
    ds.SeriesInstanceUID = pydicom.uid.generate_uid()
    ds.SeriesDate = date.today()
    ds.AcquisitionDate = date.today()
    ds.ContentDate = date.today()

    ds.is_little_endian = True
    ds.is_implicit_VR = False

    return ds

def getListBoundary(mask):
    # Create Retangle to see if the code would work in theory
    mask = np.transpose(mask)
    # prediction = np.zeros(mask.shape)
    # prediction[200:500,150:250]=1
    grad = np.gradient(mask)
    grad = np.abs(grad[0])+ np.abs(grad[1]) # Probably could be expressed nicer or more in programming style! But who cares really!
    # bound = np.argwhere(grad>0)
    bound = np.argwhere(grad>0)
    # boundSorted = sortContourPoints(bound)# In order to write the contour correctly we need to order it like we would probably contour it based on the distance to neigbhoring points

    # try using a convex hull instead potentially the order of the points are better
    #hull = convex_hull(bound,np.int(len(bound)/3))
    #print(len(bound))
    alpha_shape = alphashape.alphashape(bound, 0.1)
    vertices = np.array(list(alpha_shape.boundary.coords))
    # hull = ConvexHull(bound, incremental= True)
    # vertices = find_concave_hull(bound)
    return vertices[::2,:]

def sortContourPoints(bound):
    # Points are in Pixel space as of now!!! Is that a problem?
    sortedBound = np.zeros(bound.shape)
    sortedBound[0,:] = bound[0,:] # random sart point! Could be anywhere
    bound_new = np.delete(bound,0,0) # delete the first row
    for i in range(0, len(sortedBound)):
        #print(i)
        if(len(bound_new)>1):
            rep =  np.full(bound_new.shape, sortedBound[i,:], dtype=float)
            idx = (np.linalg.norm(rep-bound_new, axis = 1)).argmin()
            sortedBound[i,:] = bound_new[idx,:]
            bound_new = np.delete(bound_new,idx,0)
        else:
            sortedBound[i,:] = bound_new
    return sortedBound

def writeDicomREG(outdcmpath, patID,  dicomMRIFolder):
    duid = pydicom.uid.generate_uid()
    for filename in os.listdir(dicomMRIFolder):
        f = os.path.join(dicomMRIFolder, filename)
        if os.path.isfile(f):
            #print(f)
            MRI_slice_meta = pydicom.dcmread(f)
            break

    file_meta = Dataset()
    #file_meta.MediaStorageSOPClassUID = '1.2.246.352.70.2.1.160.3'
    file_meta.FileMetaInformationGroupLength = 198
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.66.1'
    file_meta.MediaStorageSOPInstanceUID = duid
    file_meta.ImplementationClassUID = '2.16.840.1.114362.1'
    file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian # in MIM implicit-- i dont think it makes a difference
    ds = FileDataset(patID, {}, file_meta=file_meta, preamble=b"\0"*128)

    ds.SpecificCharacterSet = 'ISO_IR 192'
    ds.InstanceCreationDate = date.today()
    ds.InstanceCreationTime = '000000'
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.66.1' # REGFile
    #ds.SOPInstaneUID = '2.16.840.1.114362.1.12156802.25748924467.629960731.346.26'
    duid = pydicom.uid.generate_uid()
    ds.SOPInstanceUID = duid
    ds.StudyDate = date.today()
    ds.StudyTime = '000000'
    ds.SeriesDate = date.today()
    ds.SeriesTime = '000000'
    ds.ContentDate = date.today()
    ds.ContentTime = '000000'
    ds.AccessionNumber = ''
    ds.Modality = 'REG'
    #ds.Manufacture = 'MIM' # probably doens't matter
    #ds.SecondaryCaptureDeviceManufacturer = 'Python 3.x'
    ds.StudyDescription = 'PATH VOLUMES Registered'
    ds.SeriesDescription = 'REG dummy'
    ds.PatientID = patID
    ds.PatientName = 'ProsPath'
    ds.PatientBirthDate = '20010101'
    ds.PatientBirthTime = '000000'
    ds.PatientSex = 'M'
    ds.PatientAge = ''
    
    ds.StudyInstanceUID = MRI_slice_meta.StudyInstanceUID
    ds.SeriesInstanceUID = duid
    
    ds.StudyID = MRI_slice_meta.StudyID
    ds.SeriesNumber = 1
    ds.InstanceNumber = 1
    ds.FrameOfReferenceUID = MRI_slice_meta.FrameOfReferenceUID
    ds.PositionReferenceIndicator = ''
    ds.ContentLabel = 'Registration'
    ds.ContentDescription = 'Fusion Registration'

    ds.ReferencedSeriesSequence = pydicom.sequence.Sequence()
    RSS_MRI = Dataset()
    RSS_MRI.SeriesInstanceUID = MRI_slice_meta.SeriesInstanceUID
    RSS_MRI.ReferencedInstanceSequence = pydicom.sequence.Sequence()

    for filename in os.listdir(dicomMRIFolder):
        RIS = Dataset()
        f = os.path.join(dicomMRIFolder, filename)
        if os.path.isfile(f):
            #print(f)
            MRI_slice_meta = pydicom.dcmread(f)
            RIS.ReferencedSOPClassUID = MRI_slice_meta.SOPClassUID
            RIS.ReferencedSOPInstanceUID = MRI_slice_meta.SOPInstanceUID
            RSS_MRI.ReferencedInstanceSequence.append(RIS)
    ds.ReferencedSeriesSequence.append(RSS_MRI)


    ds.StudiesContainingOtherReferencedInstancesSequence = pydicom.sequence.Sequence()
    SCORIS = Dataset()
    SCORIS.ReferencedSeriesSequence = pydicom.sequence.Sequence()
    RSS_WMHP = Dataset()
    RSS_WMHP.ReferencedInstanceSequence = pydicom.sequence.Sequence()
    for filename in os.listdir(outdcmpath):
        f = os.path.join(outdcmpath, filename)
    # checking if it is a file
        if os.path.isfile(f) and ("RS_hist" not in f) and ("REG" not in f):
            #print(f)
            WMHP_slice_meta = pydicom.dcmread(f)
            RIS = Dataset()
            RIS.ReferencedSOPClassUID = WMHP_slice_meta.SOPClassUID
            RIS.ReferencedSOPInstanceUID = WMHP_slice_meta.SOPInstanceUID
            RSS_WMHP.ReferencedInstanceSequence.append(RIS)

    try:
        RSS_WMHP.SeriesInstanceUID = WMHP_slice_meta.SeriesInstanceUID
        SCORIS.ReferencedSeriesSequence.append(RSS_WMHP)
        SCORIS.StudyInstanceUID = WMHP_slice_meta.StudyInstanceUID
        ds.StudiesContainingOtherReferencedInstancesSequence.append(SCORIS)
    except:
        print("No moved Histopath found!")
        
    WMHP_Block = Dataset()
    WMHP_Block.FrameOfReferenceUID = WMHP_slice_meta.FrameOfReferenceUID
    WMHP_Block.ReferencedImageSequence = pydicom.sequence.Sequence()
    for filename in os.listdir(outdcmpath):
        f = os.path.join(outdcmpath, filename)
    # checking if it is a file
        if os.path.isfile(f) and ("RS_hist" not in f) and ("REG" not in f):
            #print(f)
            WMHP_slice_meta = pydicom.dcmread(f)
            RIS = Dataset()
            RIS.ReferencedSOPClassUID = WMHP_slice_meta.SOPClassUID
            RIS.ReferencedSOPInstanceUID = WMHP_slice_meta.SOPInstanceUID
            WMHP_Block.ReferencedImageSequence.append(RIS)

    WMHP_Block.MatrixRegistrationSequence = pydicom.sequence.Sequence()

    MRS = Dataset()
    MRS.MatrixSequence = pydicom.sequence.Sequence()
    MRS.RegistrationTypeCodeSequence = pydicom.sequence.Sequence()

    MS = Dataset()
    MS.FrameOfReferenceTransformationMatrixType = 'RIGID'
    MS.FrameOfReferenceTransformationMatrix = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    MRS.MatrixSequence.append(MS)

    RTCS = Dataset()
    RTCS.CodeValue = '125025'
    RTCS.CodingSchemeDesignator ='DCM'
    RTCS.CodeMeaning ='Alignment'
    MRS.RegistrationTypeCodeSequence.append(RTCS)

    WMHP_Block.MatrixRegistrationSequence = pydicom.sequence.Sequence([MRS])
    
    MRI_Block = Dataset()
    MRI_Block.FrameOfReferenceUID = MRI_slice_meta.FrameOfReferenceUID
    MRI_Block.ReferencedImageSequence = pydicom.sequence.Sequence()
    for filename in os.listdir(dicomMRIFolder):
        f = os.path.join(dicomMRIFolder, filename)
    # checking if it is a file
        if os.path.isfile(f):
            #print(f)
            MRI_slice_meta = pydicom.dcmread(f)
            RIS = Dataset()
            RIS.ReferencedSOPClassUID = MRI_slice_meta.SOPClassUID
            RIS.ReferencedSOPInstanceUID = MRI_slice_meta.SOPInstanceUID
            MRI_Block.ReferencedImageSequence.append(RIS)
    
    MRI_Block.MatrixRegistrationSequence = pydicom.sequence.Sequence([MRS])

    
    ds.RegistrationSequence = pydicom.sequence.Sequence([MRI_Block,WMHP_Block])
    pydicom.filewriter.dcmwrite(
        os.path.join(outdcmpath,'REG.dcm'), ds, write_like_original=False)
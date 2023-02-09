#!/usr/local/bin/python3
'''
 Example Usage:
 python3 parquet_from_thermo_raw.py ~/path/MA4358_FFPE_HPVpos_01_071522.raw ~/path_to_dlls/ ['ITMS','full ms'];
 
'''

############
#Parse Input
############
import argparse
import sys
import os 
from os import listdir
from os.path import isfile, join

#dir_path = os.path.dirname(os.path.realpath(__file__))
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    
    #A path to a thermo .raw file OR a path to a directory containing thermo raw files
    #If a path to a directory, all files ending in '.raw' in that directory are converted
    parser.add_argument("raw_dir", help="transition_list.", type=str)
    
    import math
    # Optional arguments
    parser.add_argument("-d", "--thermo_dlls",
                        help="path containing thermo dll's", 
                        type=str, 
                        default=os.path.dirname(os.path.realpath(__file__)))
    
    #Maximum angle between light/heavy peak area pairs. If transitions y3+, y4+, and b4+ have peak areas [10, 20, 30] for the light precursor
    #and [100, 200, 300] for the heavy precursor, then the angle is zero radians. 
    parser.add_argument("-sf", "--scan_filter_regex_list", 
                        help="list of regex that match scans to omit", 
                        nargs = "*",
                        type=str, 
                        default=['ITMS'])
    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()

    return args

#Want to report the analysis time
import time
time0 = time.time()


args = parseArguments()
out_path = './'
#Print arguments 
print("Raw Files Directory:", args.raw_dir)
print("Thermo dlls Directory: ", args.thermo_dlls)
print("Scan Filter Regex List: ", args.scan_filter_regex_list)
#print("apex percent: ", args.apex_percent)

############
#Import Dependencies
############

#From pythonnet package. Needed to use ThermoFisher.CommonCore dll's
#import pythonnet
#from pythonnet import load
#print(pythonnet.get_runtime_info())
#load("coreclr")
#print(pythonnet.get_runtime_info())
import clr
#Paths to thermo dlls
thermo_data_path = args.thermo_dlls + '/ThermoFisher.CommonCore.Data.dll'
thermo_rawfilereader_path = args.thermo_dlls + '/ThermoFisher.CommonCore.RawFileReader.dll'
#Throw an exception if the dll's cannot be found and print the path
#Otherwise, import the dll
for path in [thermo_data_path, thermo_rawfilereader_path]:
    if not os.path.isfile(path):
        raise Exception("The file: " + path + " does not exist!")
    else:
        clr.AddReference(path)

#Import required functions from the Thermo Dlls
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data.Business import Device, Scan
from ThermoFisher.CommonCore.Data.Interfaces import IScanFilter
from ThermoFisher.CommonCore.Data.FilterEnums import MassAnalyzerType, MSOrderType

############
#Get and Check Raw File Paths
############

#Make sure the path to the raw file 
#if not os.path.isfile(args.raw_file):
#    raise Exception("The raw file: " + args.raw_file + " could not be found!")
raw_file_paths = []
#The path is to a single .raw file
if args.raw_dir.endswith(".raw"):
    raw_file_paths += [args.raw_dir]
#The path is to a folder containing many raw files
else:
    raw_file_paths = [join(args.raw_dir, f) for f in listdir(args.raw_dir) if isfile(join(args.raw_dir, f)) and f.endswith('.raw')]
    if len(raw_file_paths) == 0:
        raise Exception("There were no .raw files in the path: " + args.raw_dir)
    

############
#Format scan filters
############
import re
scan_filters = [re.compile(scanFilter) for scanFilter in args.scan_filter_regex_list]
############
#Convert Raw Files
############
import time
import struct
import pyarrow as pa
import pyarrow.parquet as pq
import gc
import numpy as np

#Scan table holds all scans for a .raw file. 
#For example the key 'masses' has a value, that is a 
#list of lists. That is, one mass list for each scan in the .raw file
#Give a copy of scan table to each instance of "readRawFile"
class Scans(object):
    def __init__(self):
        self.__Dict__ = {
                                'fileName' : [],
                                'basePeakMass' : [], 
                                'scanType' : [], 
                                'basePeakIntensity' : [],
                                'packetType' : [],
                                'scanNumber' : [], 
                                'retentionTime' : [],
                                'masses' : [],
                                'intensities' : [],
                                'lowMass' : [],
                                'highMass' : [],
                                'TIC' : [],
                                'FileID' : [],
                                'precursorMZ' : [],
                                'precursorCharge' : [], 
                                'msOrder': [],
                            }
        self.__PaSchema__ = pa.schema([
                                        pa.field('fileName', pa.string()),
                                        pa.field('basePeakMass', pa.float32()),
                                        pa.field('scanType', pa.string()),
                                        pa.field('basePeakIntensity', pa.float32()),
                                        pa.field('packetType', pa.int32()),
                                        pa.field('scanNumber', pa.int32()),
                                        pa.field('retentionTime', pa.float32()),
                                        pa.field('masses', pa.large_list(pa.float32())),
                                        pa.field('intensities', pa.large_list(pa.float32())),
                                        pa.field('lowMass', pa.float32()),
                                        pa.field('highMass', pa.float32()),
                                        pa.field('TIC', pa.float32()),
                                        pa.field('FileID', pa.string()),
                                        pa.field('precursorMZ', pa.float32()),
                                        pa.field('precursorCharge', pa.int32()),
                                        pa.field('msOrder', pa.int32())
                                    ])      
        self.__PaTable__ = None

    def addScan(self, scan_stats, centroid_stream, file_name, ms_order, scan_number, precursor_MZ = None, precursor_Charge = None):
        try:
            self.__Dict__['fileName']  += [file_name]
            self.__Dict__['basePeakMass']  += [scan_stats.BasePeakMass]
            self.__Dict__['scanType']  += [scan_stats.ScanType]
            self.__Dict__['basePeakIntensity']  += [scan_stats.BasePeakIntensity]
            self.__Dict__['packetType']  += [scan_stats.PacketType]
            self.__Dict__['scanNumber']  += [ scan_stats.StartTime]
            self.__Dict__['retentionTime']  += [scan_stats.StartTime]
            self.__Dict__['masses']  += [centroid_stream.Masses]
            self.__Dict__['intensities']  += [centroid_stream.Intensities]
            self.__Dict__['lowMass']  += [scan_stats.LowMass]
            self.__Dict__['highMass']  += [scan_stats.HighMass]
            self.__Dict__['TIC']  += [scan_stats.TIC]
            self.__Dict__['FileID']  += [file_name]
            self.__Dict__['precursorMZ']  += [precursor_MZ]
            self.__Dict__['precursorCharge']  += [precursor_Charge]
            self.__Dict__['msOrder']  += [ms_order]
        except:
            print("scan number ", str(scan_number))
        return self
    
    def __toPaTable__(self, f_out):
        self.__PaTable__ = pa.Table.from_pydict(self.__Dict__, schema = self.__PaSchema__)
        return self 

    def writeParquet(self, f_out):
        self.__toPaTable__(f_out)
        pq.write_table(self.__PaTable__, 
                       where = f_out, #File name out
                       compression = 'SNAPPY' 
                       )
        return

def filterScan(scan_event_string, scan_filters):
    '''
    Test if the scan filter string matches any of the 
    "scan_filters". These are regular expressions
    If there is a match, the scan will be skipped. 
    Could alternatively use "ThermoFisher.CommonCore.Data.FilterEnums"
    '''
    for scan_filter in scan_filters:
        if scan_filter.search(scan_event_string):
            return True
    return False

def GetCentroidedScan(rawFile, scan_filter, scan_number):

    centroid_stream = rawFile.GetCentroidStream(scan_number, False)

    if scan_filter.MassAnalyzer ==MassAnalyzerType.MassAnalyzerITMS:
        scan_data = rawFile.GetSimplifiedScan(scan_number)
        centroid_stream.Masses = np.array(scan_data.Masses)
        centroid_stream.Intensities = np.array(scan_data.Intensities)

    return centroid_stream


from tqdm import tqdm

def convertRawFile(raw_file_path, scan_filters, out_path):
    f_out = raw_file_path.split('/')[-1].split('.')[0]
    print("Processing raw file: ", f_out)
    print("Reading Scans...")
    #Instance of Scans class
    scans = Scans()

    #Import the raw file and select instrument
    rawFile = RawFileReaderAdapter.FileFactory(raw_file_path)
    rawFile.SelectInstrument(Device.MS, 1)
    #Get the firt and last scan indices
    first_scan_number = rawFile.RunHeaderEx.FirstSpectrum
    last_scan_number = rawFile.RunHeaderEx.LastSpectrum

    #for scan_number in range(45000, last_scan_number): no error
    for scan_number in tqdm(range(first_scan_number, last_scan_number)): 
        #print(scan_number)
    #for scan_number in range(37500, 45000): no error
    #for scan_number in range(30000, 37500): error
        scan_filter = rawFile.GetFilterForScanNumber(scan_number)
        scan_stats = rawFile.GetScanStatsForScanNumber(scan_number)
                #If scan matches one of the filters, skip it
        if filterScan(rawFile.GetScanEventStringForScanNumber(scan_number), scan_filters):
            continue
        

        centroid_stream = GetCentroidedScan(rawFile, scan_filter, scan_number)

        msOrder = rawFile.GetScanEventForScanNumber(scan_number).MSOrder

        trailers = rawFile.GetTrailerExtraValues(scan_number)
        trailer_labels = rawFile.GetTrailerExtraInformation(scan_number)
        charge_state = 0

        for i in range(trailer_labels.Labels.Length):
            if trailer_labels.Labels[i] == "Charge State:":
                charge_state = rawFile.GetTrailerExtraValue(scan_number, i)
                break

        if msOrder == MSOrderType.Ms:
            scans.addScan(scan_stats, centroid_stream, raw_file_path, msOrder, scan_number,
                          None, charge_state)
        else:
            scans.addScan(scan_stats, centroid_stream, raw_file_path, msOrder, scan_number,
                          rawFile.GetScanEventForScanNumber(scan_number).GetReaction(0).PrecursorMass,
                          charge_state)
    print(join(out_path, f_out)+'.parquet')
    print("Writing scans to parquet for ", f_out, " ...")
    t0 = time.time()
    scans.writeParquet(join(out_path, f_out)+'.parquet')
    print("Writing raw file ", f_out, " took ", str(time.time()-t0) , " seconds")
    rawFile.Dispose()
    del scans #still error
    gc.collect()
    return 

from itertools import repeat
#with mp.Pool(4) as pool: 
#    args = zip(raw_file_paths, repeat(scan_filters), repeat(out_path))
#    pool.starmap(convertRawFile, args)
    #pool.map(lambda p: (p, convertRawFile(p, scan_filters, out_path)), raw_file_paths)

initial = time.time()
for raw_file_path in raw_file_paths:
    convertRawFile(raw_file_path, scan_filters, out_path)
print("Converted " + str(len(raw_file_paths)) + " raw files in " + str((time.time() - initial)/60) + " minutes")
# pandas
#import psutil
#############
#WARNING! DO NOT REMOVE THE FOLLOWING LINE OF CODE!
#ON MAC OS, REMOVING "import psutil" CAUSES A SEGV ERROR IN 
#THE MONO RUNTIME. PART OF THE ERROR MESSAGE IS PASTED BELOW. 
#############
print("Memory in use (MB): ")
import psutil; print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
#PART OF THE ERROR RECIEVED IF YOU REMOVE "import psutil"
#####
#=================================================================
#        Native Crash Reporting
#=================================================================
#Got a segv while executing native code. This usually indicates
#a fatal error in the mono runtime or one of the native libraries 
#used by your application.
#####
sys.exit()
#import multiprocessing as mp
#with mp.Pool(4) as pool: 
#    dict(pool.map(lambda p: (p, convertRawFile(p, scan_filters, out_path)), raw_file_paths))
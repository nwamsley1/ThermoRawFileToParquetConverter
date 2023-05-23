#!/usr/local/bin/python3
'''
 Example Usage:
 python3 parquet_from_thermo_raw.py ~/path/MA4358_FFPE_HPVpos_01_071522.raw ~/path_to_dlls/ -sf ITMS cid 
 python3 raw_to_parquet.py args.json

 Uses Thermo RawFileReader .NET Assemblies (NetStandard20 .dlls) found in https://github.com/thermofisherlsms/RawFileReader
 to convert Thermo '.raw' files to Apache '.parquet files'. Features parallel processing. 

At present there are at least 2 major limitations:

    1) Will not centroid profile scans. Reading, convertine, and writing, profile mode scans (ITMS scans)
    is slow and not efficient in terms of memory or disk space. 

    2) Does not put info about pressure traces into the .parquet file. Only retrieves and writes
    information from the first 'MS' device. 



To convert Thermo Raw files to. See "parseArguments()" for details on implementation
'''

############
#Parse Input
############
import argparse
import sys
import os 
from os import listdir
from os.path import isfile, join
import time

#dir_path = os.path.dirname(os.path.realpath(__file__))
def parseArguments():
    #In the future, need to make optional to provide reference to a single .json file with all the arguments 
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
    
    #Omit any scans where the scan filter matches at least one of these regex expressions
    #Example ['ITMS'] means that the scan 
    #ITMS + p NSI t Full ms [300.0000-1100.0000]
    #Would be skipped 
    parser.add_argument("-sf", "--scan_filter_regex_list", 
                        help="list of regex that match scans to omit", 
                        nargs = "*",
                        type=str, 
                        default=['ITMS'])

    #The number of raw files to process in parallel. Uses 'multiprocessing' python module
    parser.add_argument("-o", "--parquet_out",
                        help="path to folder in which to write parquet files", 
                        type=str, 
                        default=os.path.dirname(os.path.realpath(__file__)))

    parser.add_argument("-n", "--num_workers",
                        help="number of workers to use", 
                        type=int, 
                        default=4)


    parser.add_argument

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()

    return args




#Two opsions
import json 
args = parseArguments()
if args.raw_dir.split('.')[-1] == "json":

    json_args_f = open(args.raw_dir)
    json_args = json.load(json_args_f)

    try:
        args.raw_dir = json_args['raw_dir']
        args.thermo_dlls = json_args['thermo_dlls']
        args.scan_filter_regex_list = json_args['scan_filter_regex_list']
        args.num_workers = json_args['num_workers']
        args.parquet_out= json_args['parquet_out']
    except:
        print("Could not convert json_args to properly formated arguments")
else:
    #If json was not provided as input, make a json file of the arguments. 
    #Could simply provide a path to this file in the future as the sole argument
    #To replicate the conversion
    with open('args.json', 'w') as outfile:
        json_args = {
                 'raw_dir': args.raw_dir,
                 'thermo_dlls': args.thermo_dlls,
                 'scan_filter_regex_list': args.scan_filter_regex_list,
                 'num_workers': args.num_workers,
                 'parquet_out': args.parquet_out            
                }
        json.dump(json_args, outfile)

#Print arguments 
def printArguments():
    print("Raw Files Directory:", args.raw_dir)
    print("Thermo dlls Directory: ", args.thermo_dlls)
    print("Scan Filter Regex List: ", args.scan_filter_regex_list)
    print("Number of Workers: ", args.num_workers)
    print("Parquet File Output Folder: ", args.parquet_out)
    return 


############
#Import Thermo Dependencies
############
import clr #From pythonnet package
#Paths to thermo dlls
import os
from os.path import abspath
thermo_data_path = abspath(os.path.join(args.thermo_dlls, 'ThermoFisher.CommonCore.Data.dll'))
thermo_rawfilereader_path = abspath(os.path.join(args.thermo_dlls,  'ThermoFisher.CommonCore.RawFileReader.dll'))

#Throw an exception if the dll's cannot be found and print the path
#Otherwise, import the dll
for path in [thermo_data_path, thermo_rawfilereader_path]:
    if not os.path.isfile(path):
        raise Exception("The file: " + path + " does not exist!")
    else:
        print("PATH")
        print(path)
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
#Make folder for parquet files
############

#Global variable for file output path
parquet_out = abspath(args.parquet_out)
if not os.path.exists(parquet_out):
    os.makedirs(parquet_out)


############
#Convert Raw Files
############
import struct
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv as csv
import pyarrow.dataset as dataset
import gc
import numpy as np

#Scan table holds all scans for a .raw file. 
#For example the key 'masses' has a value, that is a 
#list of lists. That is, one mass list for each scan in the .raw file
#Give a copy of scan table to each instance of "readRawFile"
class Scans(object):
    """
    Scans object represents scans from a raw file

    Attributes:
        __Dict__        'dict': Representation of scans form a raw file. Has a key for each column. 
                        Values are lists with an entry for each row(scan). 

        __PaSchema__    'pyarrow.schema': A collection of fields that corresponds to self.__Dict__
                        Used for writing parquet files

        __PaTable__     'pyarrow.Table': Columnar representation of scans data from self.__Dict__

    Methods:
        addScan():      Appends one value to each field in self.__Dict__ according to supplied arguments       
        __toPaTable__():Converts self.__Dict__ to a pyarrow.Table object self.__PaTable__ 
        writeParquet:   Writes self.__PaTable__ to a parquet file
    """
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
        self.__Dict__['fileName']  += [file_name]
        self.__Dict__['basePeakMass']  += [scan_stats.BasePeakMass]
        self.__Dict__['scanType']  += [scan_stats.ScanType]
        self.__Dict__['basePeakIntensity']  += [scan_stats.BasePeakIntensity]
        self.__Dict__['packetType']  += [scan_stats.PacketType]
        self.__Dict__['scanNumber']  += [scan_number]
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
        return self
    
    def __toPaTable__(self, f_out):
        try:
            self.__PaTable__ = pa.Table.from_pydict(self.__Dict__, schema = self.__PaSchema__)
        except:
            print("Conversion of self.__Dict__ to pyarrow.Table failed!")
        return self 

    def writeParquet(self, f_out):

        if self.__PaTable__ is None:
            self.__toPaTable__(f_out)

        #pq.write_table(self.__PaTable__, 
        #               where = f_out, #File name out
        #               compression = 'SNAPPY' 
        #               )
        from pyarrow import fs
        local = fs.LocalFileSystem()
        with local.open_output_stream(f_out) as file:
            with pa.RecordBatchFileWriter(file, self.__PaTable__.schema) as writer:
                writer.write_table(self.__PaTable__)
        #dataset.write_dataset(data = self.__PaTable__, 
        #       base_dir = '/Users/n.t.wamsley/Projects/RawToParquetPython/ThermoRawFileToParquetConverter/parquet_out/',#f_out,#, #File name out,
        #       format = "arrow",
        #       use_threads = True
        #       #compression = 'SNAPPY' 
        #       )
        return

def filterScan(scan_event_string, scan_filters):
    '''
    Assumptions:
        scan_event_string       str: A string describing the scan event. From GetScanEventStringForScanNumber()
                                method on ThermoFisher.CommonCore.Data.Interfaces.IRawDataExtended object
                                Examples: 
                                            "FTMS + p NSI Full ms [300.0000-1100.0000]"
                                            "FTMS + p NSI Full ms2 710.4059@hcd32.00 [150.0000-1700.0000]"
                                            "ITMS + p NSI t Full ms [300.0000-1100.0000]"

        scan_filters            list() of re.Pattern: List of re.Pattern objects. 

    Guarantees:
        Searches the 'scan_event_string' for matches to one or more of the scan_filters (re.Pattern objects). 
    Returns True if there is at least one match to at least one scan filter
    Returns false otherwise
    Could alternatively use "ThermoFisher.CommonCore.Data.FilterEnums"
    '''
    for scan_filter in scan_filters:
        if scan_filter.search(scan_event_string):
            return True
    return False

def GetCentroidedScan(rawFile, scan_filter, scan_number):
    '''
    Assumptions:
        rawFile                 'ThermoFisher.CommonCore.Data.Interfaces.IRawDataExtended':
                                A handle on a thermo raw file. 

        scan_filter             ThermoFisher.CommonCore.Data.Interfaces.IScanFilter:
                                Scan filter for the given scan number

        scan_number             int: Integer refereing to the scan number to get.              
        
    Guarantees:
        Gets and returns a centroid_stream for the given scan number from the rawFile.
        Returns a 'ThermoFisher.CommonCore.Data.Business.CentroidStream' object. 
        In the future will centroid profile mode scans. 
    '''
    centroid_stream = rawFile.GetCentroidStream(scan_number, False)

    #At present, very slow for ITMS (linar ion trap) scans because they are not centroided
    #by default. Ineficient in terms of memory, disk, space, and file-conversion time
    #Need a method that can quickly centroid scans profile data, perhaps using numpy
    #This block of code is avoided by passing 'ITMS' as a -sf argument. 
    if scan_filter.MassAnalyzer ==MassAnalyzerType.MassAnalyzerITMS:
        scan_data = rawFile.GetSimplifiedScan(scan_number)
        centroid_stream.Masses = np.array(scan_data.Masses)
        centroid_stream.Intensities = np.array(scan_data.Intensities)

    return centroid_stream


from tqdm import tqdm

def convertRawFile(raw_file_path, scan_filters, out_path):
    '''
    Assumptions:
        raw_file_paths          str: path to a Thermo '.raw' file

        scan_filters            list() of re.Pattern: regular expressions
                                used to filter match unwanted scans
                                that will not be processed

        out_path                str: path to folder where .parquet files
                                will be written

    Guarantees:
        Converts the a '.raw' file 'raw_file_path' to an Apache '.parquet' file in the 'out_path'
    but skips scans that match one or more of the 'scan_filters'.
    '''
    #Parse the file name out of the file path and remove the '.raw'
    #f_out = raw_file_path.split('/')[-1].split('.')[0]
    f_out = os.path.split(raw_file_path)[-1].split('.')[0]
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

    #Loop through all scans in the rawFile and read each that passes the scan_filters
    #Into a "Scans" object. 
    for scan_number in tqdm(range(first_scan_number, last_scan_number)): 

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
        try:
            if msOrder == MSOrderType.Ms:
                scans.addScan(scan_stats, centroid_stream, raw_file_path, msOrder, scan_number,
                              None, charge_state)
            else:
                scans.addScan(scan_stats, centroid_stream, raw_file_path, msOrder, scan_number,
                              rawFile.GetScanEventForScanNumber(scan_number).GetReaction(0).PrecursorMass,
                              charge_state)
        except:
            print("Could not add scan # " + str(scan_number))

    print("Writing scans to parquet for ", f_out, " ...")

    #Time converting "Scans" object to pyarrow.Table and writing to a .parquet file
    t0 = time.time()
    scans.writeParquet(join(out_path, f_out)+
        #'.parquet'
        '.arrow')
    print("Writing raw file ", f_out, " took ", str(time.time()-t0) , " seconds")

    #Dispose of rawFile and clear memory
    rawFile.Dispose()
    del scans #still error
    gc.collect()
    return 

def main():

    #Time for converting all files
    initial = time.time()

    printArguments()

    #Convert raw files in parallel
    from itertools import repeat
    import multiprocessing as mp
    with mp.Pool(args.num_workers) as pool: 
        arguments = zip(raw_file_paths, repeat(scan_filters), repeat(parquet_out))
        pool.starmap(convertRawFile, arguments)

    #Time for converting all files
    print("Converted " + str(len(raw_file_paths)) + " raw files in " + str((time.time() - initial)/60) + " minutes")

    print("Memory in use (MB): ")
    #############
    #WARNING! DO NOT REMOVE THE FOLLOWING LINE OF CODE!
    #ON MAC OS, REMOVING "import psutil" CAUSES A SEGV ERROR IN 
    #THE MONO RUNTIME. PART OF THE ERROR MESSAGE IS PASTED BELOW. 
    #############
    #PART OF THE ERROR RECIEVED IF YOU REMOVE "import psutil"
    #####
    #=================================================================
    #        Native Crash Reporting
    #=================================================================
    #Got a segv while executing native code. This usually indicates
    #a fatal error in the mono runtime or one of the native libraries 
    #used by your application.
    #####
    import psutil; print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)

    return 0

if __name__ == '__main__':
    sys.exit(main())

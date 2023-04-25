# ThermoRawFileToParquetConverter

### Cross platform conversion of thermo raw files to the parquet file format using pythonnet

 Uses Thermo RawFileReader .NET Assemblies (NetStandard20 .dlls) found in https://github.com/thermofisherlsms/RawFileReader
 to convert Thermo '.raw' files to Apache '.parquet files'. Features parallel processing. Given a directory containing Thermo *.raw files, converts each to
 a .arrow or .parquet format in a specified output folder. 
 
 The output files have the following fields with one entry per scan in the *.raw file. 
 |Name                |Type                |Description                    |
 |--------------------|--------------------|--------------------|
 |fileName            |String|Name of the original *.raw file with the suffix <br> removed and the .arrow suffix added
 |basePeakMass        |Float32|Mass of the most intense peak in the spectrum 
 |scanType            |String|FTMS or ITMS
 |basePeakIntensity   |Float32|Intensity of the most intense peak in the spectrum
 |packetType          |Int32|
 |scanNumber          |Int32|Scan index of i'th scan in the .*raw file in order of occurence
 |retentionTime       |Float32|Retention time recorded for the scan 
 |masses              |LargeList{Float32}|List of masses for peaks in the centroided <br> spectra in ascending order
 |intensities         |LargeList{Float32}|List of intensities for peaks in the centroiede <br> spectra. Corresponds with the peaks in `masses`
 |lowMass             |Float32|First mass in the scan
 |highMass            |Float32|Last mass in the scan
 |TIC                 |Float32|Summed intensity of all peaks in the scan
 |FileID              |String|Index of the *.raw file in order of processing
 |precursorMZ         |Float32|Precursor MZ or the center of the isolation window for an MSN scan <br> If no precursor was assigned. Missing for an MS1 scan.
 |precursorCharge     |Int32|Precursor charge
 |msOrder             |Int32|As in MS1, MS2, or MSN scan. Is "2" for an MS2 scan. 

 # Usage and Examples
 
 # Notes/Future Work
 At present there are at least two major limitations
 
 1) Will not centroid profile scans. Reading, converting, and writing, profile mode scans (ITMS scans)
    is slow and not efficient in terms of memory or disk space. It is recommended to only apply this tool to .raw files 
    where the data are centroided. 

 2) Does not put info about pressure traces into the .parquet file. Only retrieves and writes
    information from the first 'MS' device. This functionality is enabled by the Thermo RawFileReader but may not be possible to implement via pythonnet. 
    
   

# ThermoRawFileToParquetConverter

### Cross platform conversion of Thermo raw files to the parquet file format using pythonnet

 Uses Thermo RawFileReader .NET Assemblies (NetStandard20 .dlls) found in https://github.com/thermofisherlsms/RawFileReader
 to convert Thermo '.raw' files to Apache '.parquet files'. Features parallel processing. Given a directory containing Thermo *.raw files, converts each to
 a .arrow or .parquet format in a specified output folder. Multi-threading is supported to convert many files in parallel. Typical conversion times are 1-2 minutes per *.raw file on a single thread. Currently configured for conversion on Mac and linux. Need to run pip install pythonnet and have Mono runtime. 
 
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
 |scanHeader | MapArray{(utf8,utf8)} | Map Array of key-values of each scan header field for each scan.

 # Usage
 
 ## Convert *.raw files using command line options

1) <b>-d</b> flag sets directory to search for Thermo RawFileReader dll's
2) <b>-n</b> flag specifies the number of threads to use
3) <b>-sf</b> list of terms to search for in each scan header/filter. Scans that contain these words are ommited
4) <b>-o</b> path to folder where the converted files will be saved
5) **-sh** optional flag specified whether to parse `scanHeader` for each scan. Default is `False`, or use `--no-scan-header`.

###### Install Dependencies
```
pip install -r requirements.txt
```
###### POSIX
```
python3 raw_to_parquet.py ../raw -d ./libs -n 12 -sf ITMS kazoo -o ./parquet_out
```
###### WINDOWS
```
python .\raw_to_parquet.py ..\raw -d .\libs\ -n 12 -sf ITMS -o .\parquet_out
```
Note: For windows, downloading the repository as a .zip file will cause execution of the .dll files to be blocked. Instead, clone the repository as follows
```
git clone https://github.com/nwamsley1/ThermoRawFileToParquetConverter.git
```

This command converts all raw files in the "../raw" folder to a .arrow format, but excludes scans with "ITMS" or "kazoo" in the scan filter. An example scan filter is:  "ITMS + p NSI t Full ms [300.0000-1100.0000]". The new files are saved into the "./parquet_out" folder. In the current directory, an "args.json" file is also generated. 

```
import json
json_args_f = open('args.json')
json_args = json.load(json_args_f)
display(json_args)
{'raw_dir': '../raw',
 'thermo_dlls': '../libs',
 'scan_filter_regex_list': ['ITMS', 'kazoo'],
 'num_workers': 12,
 'parquet_out': './parquet_out',
 'scan_header_used': true}
```

## Convert *.raw files using .json arguments
In the previous example running the raw file converter using command line options generated a "args.json"
file specifying the options. Rather than supplying options through flags in the command line, a simple .json
file can be specified as follows

```
python3 raw_to_parquet.py args.json
```


 # Notes/Future Work
 At present there are at least two major limitations
 
 1) Will not centroid profile scans. Reading, converting, and writing, profile mode scans (ITMS scans)
    is slow and not efficient in terms of memory or disk space. It is recommended to only apply this tool to .raw files 
    where the data are centroided. 

 2) Does not put info about pressure traces into the .parquet file. Only retrieves and writes
    information from the first 'MS' device. This functionality is enabled by the Thermo RawFileReader but may not be possible to implement via pythonnet. 
    
 3) Currently files are converted to ".arrow". It is possible to convert to .csv or .parquet as well. See below. The desired output format should be specifiable as a command line argument in a future version.
 
 
 Parquet format
 ```
        pq.write_table(self.__PaTable__, 
                      where = f_out, #File name out
                       compression = 'SNAPPY' 
                       )
```
Arrow Format
```
        from pyarrow import fs
        local = fs.LocalFileSystem()
        with local.open_output_stream(f_out) as file:
            with pa.RecordBatchFileWriter(file, self.__PaTable__.schema) as writer:
                writer.write_table(self.__PaTable__)
 ```
4) Needs to be tested on Windows. The only errors at present should have to do with file path compatability. Should be possible to fix this without too much trouble. 

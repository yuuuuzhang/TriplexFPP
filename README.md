# TriplexFPP: A machine learning program for DNA:RNA triplex potential prediction
TriplexFPP is an integrated program for DNA:RNA triplex prediction. It contains three models, 1) triplex lncRNA prediction model, to predict the most likely triplex forming lncRNA in practical, 2) triplex DNA site prediction model, to predict if a DNA site can form triplex in practical, and 3) in cis / in trans lncRNA prediction model, to predict the interaction position of a triplex lncRNA. It will give researchers usegul guidelines to study the DNA:RNA triplex formation.

## SYSTEM REQUIREMENTS

### Hardware requirements
TriplexFPP requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements

#### OS Requirements

TriplexFPP has been tested on the following systems:

* macOS (10.14.6)
* Windows10
 with tensorflow 1 (>=1.13)

#### Python Dependencies

Based on python3.  
Python modules:  
```
numpy  
pandas  
csv  
tensorflow 
keras
Bio
math
pickle
```

## EXPLANATION
This repository contains four folders: code, embed, input_example and output.

### Code folder:
contains the python codes.  
```
triplex_util.py -- functions will be used.  
TriplexFPP.ipynb -- user interface.  
```
### embed folder:
This folder contains the models and files that will be used. The models are from one of the cross-fold validation in three prediction tasks.

### input_example folder:
This folder contains example input files. These example files are from one of the cross-fold validation test data for 3 models. The file must be in .fa or .fasta format.

### output_files folder:
The output result files will be put into this folder.


## USAGE:
  
Download all the files firstly, this takes about mins, open TriplexFPP.ipynb, change data path:  
```
datapath = '.../...'
```
to where TriplexFPP dwonloaded,
change input file path and name
```
inputfile = '.fasta' 
```
to where input file is,
change output file name
```
outputname = '.csv'
```
to the output file name, the output file will be at the datapath/output.

* To run triplex lncRNA prediction:
- use the triplex lncRNA section in TriplexFPP.ipynb
- change corresponding data path and file names
- run
the default name is the demo, it takes about mins to run, where the feature extraction can be slow.

* To run triplex DNA site prediction:
- use the triplex DNA site section in TriplexFPP.ipynb
- change corresponding data path and file names
- run
the default name is the demo, it takes about mins to run, where the feature extraction can be slow.

* To run in cis / in trans lncRNA prediction:
- use the triplex lncRNA section in TriplexFPP.ipynb
- change corresponding data path and file names
- run
the default name is the demo, it takes about mins to run, where the feature extraction can be slow.

More details can be found from [1]

## REFERANCE
[1] 
## CONTACT
If you have any inqueries, please contact YU007@e.ntu.edu.sg

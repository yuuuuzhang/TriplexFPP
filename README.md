# TriplexFPP: A machine learning program for DNA:RNA triplex potential prediction
TriplexFPP is an integrated program for DNA:RNA triplex prediction. It contains three models, 1) triplex lncRNA prediction model, to predict the most likely triplex forming lncRNA in practical, 2) triplex DNA site prediction model, to predict if a DNA site can form triplex in practical, and 3) in cis / in trans lncRNA prediction model, to predict the interaction position of a triplex lncRNA. It will give researchers usegul guidelines to study the DNA:RNA triplex formation.

## PUBLICATION
Please cite this paper if using the codes here: 

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
Based on python3.  
Python modules:  
```
numpy  
pandas  
csv  
keras
Bio
math
pickle
```
will be used. keras is backened on tensorflow.  
Download all the files firstly, open TriplexFPP.ipynb, change data path:  
```
datapath = '.../...'
inputfile = '.fasta' 
outputname = '.csv'
```

More details can be found from [1]

## REFERANCE
[1] 
## CONTACT
If you have any inqueries, please contact YU007@e.ntu.edu.sg

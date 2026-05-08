# **CDPred**
A computational approach to predict the celiac disease causing peptides using the sequence information.
## Introduction
CDPred is a bioinformatic-ware with the ability to discriminate peptides with the ability of causing celiac disease. It is developed to predict, scan, and, design the celiac disease causing peptides using sequence information only. In the standalone version, we have provided two models to make prediction such as "Main model" which is developed using experimentally-validated peptides whcih can cause celiac disease and random peptides generated using Swiss-Prot database; another model is also implemented in the backend called as "Alternate model" developed using experimentally-validated peptides whcih can cause celiac disease and experimentally-validated peptides with no ability to cause celiac disease. This method is available as five major modules such as:
- [ ] (i)   Predict     : It used machine learning model in the backend developed on amino acid composition as input features.
- [ ] (ii)  Design      : This module generates all the possible mutants of a peptide by altering a single amino acid at a time and used the model to predict its potential to cause disease.
- [ ] (iii) Scan        : This modile generates overlapping patterns from the submitted sequence of specificed length and use them to make predictions.
- [ ] (iv)  PQ Density  : In this module, density of amino acid P and Q are used to make predictions.
- [ ] (v)   Ensemble    : This module joins motif search and machine learning model to predict the potential of the submitted peptides to cause disease.

CDPred is also available as web-server at https://webs.iiitd.edu.in/raghava/cdpred. Please read/cite the content about the CDPred for complete information including algorithm behind the approach.

## Standalone
The Standalone version of transfacpred is written in python3 and following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy

## Minimum USAGE
To know about the available option for the stanadlone, type the following command:
```
python cdpred.py -h
```
To run the example, type the following command:
```
python3 cdpred.py -i example_input.fa
```
This will predict if the submitted sequences are Binders or Non-binder. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma seperated variables).

## Full Usage
```
usage: transfacpred.py [-h] 
                       [-i INPUT 
                       [-o OUTPUT]
		       [-m {1,2}]
		       [-j {1,2,3}]
		       [-t THRESHOLD]
                       [-w {9,10,11,12,13,14,15,16,17,18,19,20}]
                       [-p {3,4,5,6,7,8,9}]
		       [-d {1,2}]
```
```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -m {1,2}, --method {1,2}
                        Method Type: 1:Main (Disease Vs Random), 2:Alternate
                        (Disease Vs Non-Disease), by default 1
  -j {1,2,3,4,5}, --job {1,2,3,4,5}
                        Job Type: 1:Predict, 2:Design, 3:Scan, 4:PQ Density,
                        5: Ensemble Method, by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.32
  -w {9,10,11,12,13,14,15,16,17,18,19,20}, --winleng {9,10,11,12,13,14,15,16,17,18,19,20}
                        Window Length: 9 to 20 (scan mode only), by default 9
  -p {3,4,5,6,7,8,9}, --patleng {3,4,5,6,7,8,9}
                        Pattern Length for PQ Density Module: 3 to 9, by
                        default 5
  -d {1,2}, --display {1,2}
                        Display: 1:Only Disease-Causing Peptides, 2: All
                        peptides, by default 1
```

**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Method:** Two methods are provided to make predictions, where 1 is using Main model (Disease Vs Random), and 2 is for Alternate (Disease Vs Non-Disease) model.

**Job:** User is allowed to choose between three different modules, such as, 1 for prediction, 2 for Designing, 3 for scanning, 4 for PQ density, and 5 for ensemble method, by default its 1.

**Threshold:** User should provide threshold between 0 and 1, by default its 0.32 for other jobs except PQ Density for which the default threshold is 0.41.

**Window length**: User can choose any pattern length between 9 and 20 in long sequences. This option is available for only scanning module.

**Pattern length**: User can choose any pattern length between 3 and 9 in case of implementing PQ density module.

**Display type:** This option allow users to fetch either only celiac disease causing peptides by choosing option 1 or prediction against all peptides by choosing option 2.

CDPred Package Files
=======================
It contains following files, brief descript of these files given below

INSTALLATION                        : Installations instructions

LICENSE                             : License information

README.md                           : This file provide information about this package

model.zip                           : This zipped file contains the compressed version of model

envfile                             : This file compeises of paths for the database and blastp executable

cdpred.py                           : Main python program

example_input.fa                    : Example file contain peptide sequenaces in FASTA format

example_Predict_output.csv          : Example output file for predict module

example_Scan_output.csv             : Example output file for scan module

example_Design_output.csv           : Example output file for design module

example_PQ_Density_output.csv       : Example output file for PQ Density module

example_Ensemble_output.csv 	    : Example output file for Ensemble method module

# Reference
Tomar R., Patiyal S., Dhall A. and Raghava G. P. S. (2023) Prediction of celiac disease-associated epitopes and motifs in a protein. <a href="https://www.frontiersin.org/articles/10.3389/fimmu.2023.1056101/">Front. Immunol. 10.3389/fimmu.2023.1056101</a>

# Lieberman--Naranjo--Mihov--Velikov
 Code used to create results in Lieberman, Naranjo, Mihov, and Velikov (2024), Show Me the Receipts: B2B Payment Timeliness and Expected Returns

This repository contains code used to create the results in Lieberman, Naranjo, Mihov, and Velikov (2024), Show Me the Receipts: B2B Payment Timeliness and Expected Returns. This code is to be used in conjunction with the MATLAB asset pricing package that accompanies Novy-Marx and Velikov (2024), Assaying Anomalies. 

The order of operations to replicate the results in Lieberman, Naranjo, Mihov, and Velikov (2024) is:

1. Download and follow the instructions for setting up the MATLAB Toolkit from https://github.com/velikov-mihail/AssayingAnomalies.git
	* The results in Lieberman, Naranjo, Mihov, and Velikov (2024) use the pre-release v0.4 of the MATLAB Toolkit.
3. Download the code in this repository.
4. Add the  CSAD file from Dunn & Bradstreet to the Input folder.
5. Run lnmv.m. The script requires setting up the directories for the MATLAB asset pricing package repository and this repository. It starts a log file and calls multiple other scripts which perform the following functions:  
	* make_data.m creates some preliminary data to be used in creating the tables and figures
	* make_tables.m prints all tables in ~/Results/LNMVTablesOutput.txt
	* make_figures.m creates and stores all figures in ~/Figures/
6. Run lnmv.do. The program reads the file stata_data.csv created by make_data.m above and creates Figure 2 and Tables 4 and 6 
   

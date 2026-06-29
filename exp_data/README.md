**General Description**
This dataset contains single-cell tracking data of HeLa cervical cancer cells obtained via long-term time-lapse microscopy using a custom shear-free, diffusive microfluidic platform. The data was used for parameter estimation and testing a cellular automata (CA) model to analyze cell cycle variability and inheritance.
**File Structure**
The dataset consists of 15 individual files, each representing data from one of the 15 independent microfluidic chambers analyzed in the study.
**Data Format and Column Descriptions**
Each XLSX file contains the following columns:
•	Position X [µm]: The horizontal (X) coordinate of the cell's center of mass, measured in micrometers.
•	Position Y [µm]: The vertical (Y) coordinate of the cell's center of mass, measured in micrometers.
•	Time [h]: The elapsed time from the start of the experiment, recorded in hours.
•	Lineage ID: A unique hierarchical identifier (pedigree) encoding the cell’s division history:
o	Founder cells (cells present at the start of the observation) are assigned single-digit IDs (e.g., 1, 2, 3).
o	Daughter cells inherit their mother’s ID with an appended 1 or 2 upon division (e.g., cell 1 divides into 11 and 12; cell 11 subsequently divides into 111 and 112).
o	The total number of digits in the ID directly corresponds to the generation index of the cell.
**Dataset Partitioning**
For the purposes of the iterative parameter-testing and model training described in the manuscript, the chambers were divided as follows:
•	Training Set: Data from chambers 9, 10, and 11.
•	Test Set: Data from all other chambers (1–8 and 12–15).
**Contact Information**
For any questions regarding the dataset structure, please contact the corresponding author as indicated in the main manuscript.

This repository contains results and models for the project:

## Estimation of Greenhouse Gas Emissions (GHG) Using Big Traffic Data for Sustainable Traffic Control and Management
A joint effort between Miovision, ODX, and the Universit of Waterloo

To run all of the processes in this repository, addition data files are required. Due to Github's file size constraints, these files are shared at
https://www.dropbox.com/sh/y3vz8w00syzrmnj/AADUDsC1_1hr38RlKjGAsOWSa?dl=0
They must be placed in the repository directory.

### Prerequisites
Python 3.6 and corresponding libraries.

### About the Data Files
#### Vissim Outputs
The VissimOutputs folder contains output files from Vissim simulations, including vehicle records, signal timings, simulated detector data and simulated travel time measurements.
##### Vissim Model (Not Required)
The Vissim model, which is not required for running any part of this repository, is also shared in the above link in the folder OptionalFiles-VissimModel. This can be downloaded to run further simulations.
###### Prerequisite to Run 
PTV Vissim 7 or higher
#### MOVES Inputs
The files used as inputs for running emission simulations in MOVES include converted output data from the Vissim simulations, fuel data, meteorological data, and other specificiations.
##### Converted Vissim Outputs
The MOVESinputs folder contains files converted from Vissim vehicle records into operating mode distributions, links properties and link sources.
##### Other MOVES Inputs
Other MOVES inputs are either in the MOVESinputs folder or in the main directory, depending on their spatial-temporal specificity. 

#### MOVES Outputs 
The MOVESoutputs folder contains emission estimates from MOVES exported from its SQL database into csv files. Version: MOVES2014a

#### Files for Tableau
The ForTableau folder contains user-friendly formats of the MOVES outputs. MovesOutputForTableau.csv is the file used to generate the publicly shared Tableau workbook and visualizations at https://public.tableau.com/views/VehicleEmissionsonHespeler-DunbartoHwy401/Dashboard1

# REPLICATION REPOSITORY FOR LEPPERT (2023)
Repository to reproduce "No fences make bad neighbors, but markets make better ones: Cap and trade reduces cross-border SO2 in a natural experiment", published in Enverionmental Economics and Policy Studies, 25, pp. 407–433 (2023). Link: https://link.springer.com/article/10.1007/s10018-023-00367-z 

This repository generates the data panel for the working paper "No fences make bad neighbors, but markets make better ones: Cap and trade reduces cross-border SO2 in a natural experiment". The dataset is made up of publicly available survey- and monitoring products by the U.S. Environmental Protection Agency, the U.S. Energy Information Administration, and National Oceanic and Atmospheric Administration. Complete citations to these data are listed in the working paper bibliography. 

#### Before running this program, ensure that all associated folders and files are included as specified in the documentation and that no changes have been made to file names or data. Running the file directly from the master folder [EEPS_cross-border_SO2] available to download from [https://github.com/DanielLeppert] ensures this. 

To replicate the "data_final.csv" file, run "main.py":

## Steps:

Section 1: Loads data on plant characteristics, sulfur dioxide emissions and sulfur emissions control installations. Converts relevant variables to SI units. Restricts the data to coal-fired electric utilities 1997-2020.

Section 2: Loads weather data from the Global Historical Climatology Network which includes typical hourly observations from 460 U.S. weather stations. Computes bi-daily weighted averages for daytime and nightime.

Section 3: Links each electric utility to its most proximate weather station using a Haversine function of distance between geographic point coordinates.

Section 4: Initializes a GaussMod input matrices init_sulfur and init_weather and populates with monthly emissions-, weather- and utility variables. See inputs_readme.txt for input data documentation.

Section 5: Creates a spatial dataframe and coordinate reference system so that dispersion model output can be mapped onto geographic space centered on the plant.

Section 6: Runs GaussMod (Leppert 2022) using input matrices and computes the cross-state pollution from each utility in each month using the GaussMod output. GaussMod is a Gaussian air pollution dispersion model adapted from Abdel-Rahman (2008) and EPA. GaussMod is imported from 'gaussian_model.py' and relies on 'stability.py' to calculate atmospheric stability. The user can change the spatial resolution of GaussMod. The default is 1,000 meters. A higher resolution will result in a slower runtime.

Section 7: Adds annual variables and treatment indicators.

## Data cleaning and analysis

Run "data_cleaning&plots.R"


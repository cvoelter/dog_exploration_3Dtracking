# dog_exploration_3Dtracking
Data and R scripts for the study entitles: "Using machine learning to track dogs’ exploratory behaviour in presence and absence of their caregiver"

## 01_dog_tracking_data_processing.Rmd
* Data processing (filtering, interpolation, smoothing) and calculation of the summary stats for each trajectory (e.g., interest area durations, area covered, etc.)
* Plotting of individual key point trajectories

## 02_dog_tracking_analysis.rmd:
* Analysis of the aggregated response variables (distance travelled, tail movements, tail angles, area covered, visual angle to owner chair, duration in owner IA, duration in stranger IA, duration in door IA)
* Correlation between behavioral and survey (c-BARQ) data
* Correlation between machine-learning-based data and manual scorings

## 03_dog_manual_scorings_analysis.rmd
* Analysis of the manually scored response variables 

## saves folder
* saved analysis outputs including tables for the proportion tracked data and the model outputs

## graphics folder
* all data visualisations (Figure 2-5)
* interpolated subfolder: all preprocessed individual trajectories (for each subject, key point, and condition)

## data folder
Includes all data files supporting the analyses.
* distance_travelled_summary_data.csv: aggregated tracking data for analysis
* dog_exploration_3Dtracking_scoringdata.csv: manually scored data (raw data) of durations spent in interest areas
* dog_exploration_3Dtracking_agg_scoringdata.csv: manually scored data, aggregated by subject and condition
* 3D_tracking_object_locations_chairs_plates: location of chairs and plates for visualisation
* 3D_tracking_object_locations_chairs_plates_centers: location of chair and plate centers for analysis
* dog_tracking_project_data: demographics data and assignment to conditions
* survey subfolder: c-BARQ data
* detector2 subfolder: raw data (loopy output), one file per trajectory (for each dog and condition)

## functions folder
Functions kindly provided by Roger Mundry. 

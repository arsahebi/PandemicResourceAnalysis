# Pandemic Resource Analysis
Analysis of in-hospital code status changes during COVID-19 using R. Focuses on patient care trends and resource allocation in a Level I Trauma hospital. Provides valuable insights for healthcare data analysis


## Description
This project involves a comprehensive analysis of in-hospital code status updates, focusing on trends over time and the impact of the COVID-19 pandemic. The study, published in the "American Journal of Hospice and Palliative Medicine," utilizes data from a Level I Trauma hospital to evaluate the effectiveness of daily updates in code status documentation and its implications on patient care and resource allocation.

## Installation
To run the scripts in this repository, you'll need R installed on your system. You can download and install R from [CRAN](https://cran.r-project.org/). Once installed, you will need to install the following R packages:

```R
install.packages(c("data.table", "dplyr", "tidyr", "ggplot2", "stringr", "car", "lubridate", "emmeans", "semTools", "gridExtra", "scales", "grid", "MASS", "boot"))
```

## Data
The analysis is based on a comprehensive dataset detailing ICU and ED encounters, code status, and ventilator use across seven hospitals from March 2019 to December 2022. Note: The dataset used in this study is not publicly shareable due to privacy concerns.

## Usage
To use the script, set the `pathread` and `pathwrite` variables to the appropriate file paths where your dataset is located and where you want the output to be saved, respectively. The main script, `ProjectCodeStatus.R`, can then be run in any standard R environment.

## Function Descriptions
The script includes several custom functions for data analysis and visualization. For example, `train_sec` is used for creating secondary axis transforms in data visualizations.

## Contributors
- Amirreza Sahebi (Author)

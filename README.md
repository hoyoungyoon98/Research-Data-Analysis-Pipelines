# Research Data Analysis Pipelines

Python-based automated data analysis pipelines for experimental research data.

This project contains three analysis pipelines developed to automate the processing, statistical analysis, visualization, and export of experimental data:

- **qPCR Analysis**
  - Relative quantification
  - Absolute quantification

- **Flow Cytometry Analysis**
  - Standard analysis
  - Time-course analysis

- **Plaque Assay Analysis**
  - Viral growth curve analysis

Example datasets are included to demonstrate each analysis workflow without exposing original research data.

## Overview

Experimental data analysis often involves repetitive steps such as data preprocessing, normalization, statistical analysis, visualization, and result organization.

This project was developed to automate these workflows for qPCR, flow cytometry, and plaque assay data using Python.

Each pipeline is designed to:

- Process experimental raw data automatically
- Handle biological and technical replicates
- Perform experiment-specific calculations and normalization
- Conduct statistical analysis
- Generate publication-ready graphs
- Export analysis results to Excel files
- Process multiple experiment folders using a consistent workflow

The pipelines were developed and validated using anonymized example datasets based on experimental data structures.

## Analysis Pipelines

### 1. qPCR Analysis

The qPCR pipeline supports both relative and absolute quantification workflows.

#### Relative Quantification

The relative quantification workflow automates:

- Target and reference gene data loading
- Technical replicate processing
- ΔCt calculation
- ΔΔCt calculation using a control sample
- Relative expression calculation using 2^(-ΔΔCt)
- Paired statistical comparison against the control
- Mean and standard deviation calculation
- Bar graph generation with statistical significance
- Excel result export

#### Absolute Quantification

The absolute quantification workflow automates:

- Standard and sample data processing
- Standard curve generation
- Standard curve quality control
- Absolute copy number calculation
- Conversion to copies/mL
- Summary statistics by sample and time point
- Two-way ANOVA
- Dunnett's multiple comparison test against the control
- Absolute quantification graph generation
- Standard curve visualization
- Excel result export

### 2. Flow Cytometry Analysis

The flow cytometry pipeline processes gMFI data extracted from RTF files and supports both standard and time-course analysis workflows.

#### Standard Analysis

The standard analysis workflow automates:

- gMFI extraction from RTF files
- Sample and technical replicate identification
- Technical replicate averaging
- Normalization to a control sample
- Summary statistics across biological replicates
- Paired statistical comparison against the control
- Bar graph generation with statistical significance
- Excel result export

#### Time-course Analysis

The time-course analysis workflow automates:

- gMFI extraction from RTF files
- Sample, time point, and biological replicate identification
- Time-course data organization
- Biological replicate summary
- Optional normalization
- Statistical analysis across samples and time points
- Time-course graph generation with statistical significance
- Excel result export

### 3. Plaque Assay Analysis

The plaque assay pipeline automates viral growth curve analysis from plaque count data.

The workflow automates:

- Plaque count and dilution data loading
- Sample and time point identification
- PFU/mL calculation
- Summary statistics across biological replicates
- Log10 transformation for statistical analysis
- Two-way ANOVA
- Dunnett's multiple comparison test against the control
- Viral growth curve generation
- Limit of detection (LOD) visualization
- Statistical significance annotation
- Excel result export

## Repository Structure

```text
Research_Data_Analysis_Pipelines/
│
├── qPCR/
│   ├── qpcr_analysis.py
│   └── example_data/
│       ├── relative/
│       │   └── Experiment_001/
│       │       └── Final_Results/
│       └── absolute/
│           └── Experiment_001/
│               └── Final_Results/
│
├── Flow_Cytometry/
│   ├── flow_cytometry_analysis.py
│   └── example_data/
│       ├── standard/
│       │   └── Experiment_001/
│       │       └── Final_Results/
│       └── timecourse/
│           └── Experiment_001/
│               └── Final_Results/
│
├── Plaque_Assay/
│   ├── plaque_assay_analysis.py
│   └── example_data/
│       └── Experiment_001/
│           └── Final_Results/
│
├── README.md
├── requirements.txt
└── .gitignore
```
Each pipeline is organized independently with its analysis script and corresponding example dataset.

The example datasets reproduce the structure required by each pipeline and can be used to test the complete analysis workflow.

Example analysis outputs are also included in the `Final_Results` directories to demonstrate the expected results of each pipeline.

## Usage

### 1. Install Dependencies

Install the required Python packages using:

```bash
pip install -r requirements.txt
```

### 2. Select an Analysis

Open the corresponding Python script and configure the analysis settings in the `USER CONFIGURATION` section.

For qPCR analysis:

```python
ANALYSIS_TYPE = "relative"
```

or

```python
ANALYSIS_TYPE = "absolute"
```

For flow cytometry analysis:

```python
ANALYSIS_TYPE = "standard"
```

or

```python
ANALYSIS_TYPE = "timecourse"
```

The plaque assay pipeline runs as a single analysis workflow and does not require analysis type selection.

### 3. Prepare Input Data

Example datasets are provided in the `example_data` directory of each pipeline.

Each dataset is organized in an `Experiment_001` folder and follows the input structure required by the corresponding analysis pipeline.

The included example datasets can be used directly to test the complete analysis workflow.

### 4. Run the Pipeline

Run the corresponding Python script from the terminal.

For qPCR analysis:

```bash
python qPCR/qpcr_analysis.py
```

For flow cytometry analysis:

```bash
python Flow_Cytometry/flow_cytometry_analysis.py
```

For plaque assay analysis:

```bash
python Plaque_Assay/plaque_assay_analysis.py
```

The pipeline automatically processes the input data and generates the analysis results in the corresponding experiment folder.

## Output

Each pipeline automatically creates a `Final_Results` directory for the analyzed experiment.

The generated outputs include:

- Processed and summarized experimental data
- Statistical analysis results
- Excel files containing analysis results
- Graphs generated from the analyzed data

The specific outputs depend on the analysis pipeline.

### qPCR

- Relative or absolute quantification results
- Statistical analysis results
- Relative expression or absolute quantification graphs
- Standard curve results and visualization for absolute quantification

### Flow Cytometry

- Processed gMFI data
- Replicate and summary results
- Normalized results when applicable
- Statistical analysis results
- Standard or time-course graphs

### Plaque Assay

- Calculated PFU/mL results
- Growth curve summary
- Two-way ANOVA results
- Dunnett's multiple comparison results
- Viral growth curve graphs with LOD

## Requirements

The pipelines were developed using Python and require the following packages:

- NumPy
- pandas
- Matplotlib
- SciPy
- statsmodels
- openpyxl
- striprtf

All required packages can be installed using:

```bash
pip install -r requirements.txt
```

## Notes

These pipelines were developed based on experimental data analysis workflows used during biological research.

To protect research data, the original experimental datasets are not included in this repository. Instead, anonymized example datasets are provided to demonstrate and validate the functionality of each analysis pipeline.

The example datasets are intended for demonstration and portfolio purposes.

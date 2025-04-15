# RNA-Sequencing Set Up Tool

This Python application provides a graphical user interface (GUI) for processing RNA sample data from a CSV file. It helps automate RNA normalization calculations and generates worklists for liquid handling robots. It generates sample QC plots for run monitoring.

## Features

- Import RNA sample CSV files
- Validate and preprocess sample data
- Automatically calculate RNA and water volumes based on user input
- Assign destination wells across 384-well plates
- Generate worklists for:
  - RNA Dispensing
  - Water Dispensing
  - First Strand Synthetis Reagent addition
  - ERCC spike-ins
- Create and export summary plots:
  - RNA concentration and RIN distributions
  - Comparisons by `Fermentation Experiment` and `Library Prep Plate ID`

## Requirements

- Python 3.7+
- Libraries:
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `plotly`
  - `tkinter` 

## How to Use

1. Run the script:  
   ```bash
   python rna_seq_set_up_tool.py

# Repository for extended triplex formation analysis

## Before starting

You will want to download the rgt-genomics software.

Link here: https://reg-gen.readthedocs.io/en/latest/rgt/installation.html

Importantly, both the rgt-data and associated functions should be included in your bashrc path.

## Running the scripts 

Note about running: The scripts are designed to run on the BU cluster. As a result, most of the initially heavy lifting is performed with two qsub files. In addition, the qsub files take advantage of several modules present on the BU cluster. 

### Step 1 - Input data

### Step 2 - Setting variables/parameters

#### Step 2A - User variables for triplex identification

Four values can be set be set for the initial triplex identification step. These values are modified directly within the 'qsub' script here:

'src/00_TTS_calculation.qsub'

These include:

1) The directory - 'DIR' - The path to your version of the '2022_TriplexAnalysis' directory
2) The input folder - 'I_PATH' - Folder in 'data' with input files
3) The output folder - 'O_path' - Folder in 'output' where generated files will be placed
4) The complexity filter (boolean - on or off) - 'DEFAULT' - Determines whether the low complexity filter is off (true) or on (false)

#### Step 2B - Setting up parameter file for DBR and DBD determination


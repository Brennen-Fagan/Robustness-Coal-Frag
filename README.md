# Robustness-Coal-Frag
Code accompanying manuscript entitled "Robustness of steady state and stochastic cyclicity in generalized coalescence-fragmentation models".

# Folders
Three folder types. Data folders have raw data (created from Python code and then sparsified via MATLAB). Analysis folders have these data files summarised to relevant statistics (using R). The Code-Model folder contains the relevant scripts that were used to define, run, and analyse the simulations, as well as to organise the information for the tables and figures.

## Data Folders
Folders are organised roughly by date and theme. Each file is a MATLAB .mat file, with those that have the tag "_S" converted in MATLAB from those without as sparse matrices (see "MAT_to_Sparse" family of MATLAB scripts). Those that do not have the tag "_S" were instead created in Python using SciPy.io's savemat function (see "1_Generate" family of Python scripts). Each row of the matrix represents the population after a time step, with column i representing the number of clusters of size i. By convention, to reduce the amount of data stored, and to reduce autocorrelation (but see plots of summary statistics through time in the text), data was only stored every 10 time steps. As standard computational practice the simulations were burnt in, usually with 1,000,000 events, and were started from random partitions of the population. To facilitate saving and storing the simulation results in Python during batch jobs, we began saving (time) "slices" of the results matrix. These were then "stitched" together when the MATLAB sparsification script was run.

Note that not all data generated was used in the course of preparation of this manuscript.

The other tags used refer to

1. RandomPartition: our standard is to start the simulations from a random partition of the population size.
2. MixedKernel / AttrAccr / ImpComb / Imperfect / ImpFrag: our standard function at project completion was to run a master script that could accomodate varying coalescence and fragmentation kernels and forms (alt. perturbations). The other names are associated with intermediate functions that were specialised for varying coalescence and fragmentation forms.
3. MM / CM / MC / CC: Position indicates Coalescence and Fragmentation respectively. M indicates a multiplicative kernel. C indicates a constant kernel.
4. Bar _X_: Rarely used, but indicates a barrier at X below which fragmentation cannot occur.
5. v0 _X_: Indicates the probability of fragmentation. Here, v0 _X_ is equivalent to a P(Fragmentation) = 0.X, e.g., v001 implies P(Fragmentation) = 0.01.
6. Total / PL _X_ / BBa _X_ b _Y_ / BB _X_ AND _Y_ / Fragn _X_ d _Y_ / _X_ OF _Y_ Frag / Splitn _X_ d _Y_ / ImpComb / ImpAll / ImpCoFr / Accr _X_ Unq _Y_ / AtAc _X_ Unq _Y_ / Attr _X_ Unq _Y_ : These indicate the various simulations run.
    * Total refers to shattering fragmentation (i.e., the total breakdown of the cluster).
    * PL _X_ refers to a Chinese Restaurant Process Fragmentation that produces a Power Law with exponent _X_ / 100 (to avoid the period in the middle of a file name).
    * BBa _X_ b _Y_ / BB _X_ AND _Y_ / ImpFrag refer to beta binomial fragmentation with coefficients _X_ and _Y_. In the case of Imp(erfect) Frag(mentation), the coefficients are both set to 1 (a uniform distribution).
    * Fragn _X_ d _Y_ / _X_ OF _Y_ Frag / Splitn _X_ d _Y_ refer to fragmentation of _X_/_Y_ of the cluster chosen, with the remainder remaining clustered. When Frag(mentation) is used, the portion fragmented is shattered. When split is used, its cluster is otherwise preserved (binary fragmentation).
    * ImpComb / ImpAll / ImpCoFr refer to "imperfect" (BB 1 and 1) processes, but applied for either coalescence ("Combination") or both coalescence and fragmentation ("All" or "CoFr").
    * Accr _X_ Unq _Y_ / AtAc _X_ Unq _Y_ / Attr _X_ Unq _Y_ refer to accretion and erosion processes in addition to the existing binary coalescence and shattering fragmentation processes. Accr _X_ refers to a random accretion of _X_ individuals to clusters. Attr refers to random erosion (attrition) of _X_ individuals from clusters. AtAc _X_ combines both processes at the same rate of _X_. Unq _Y_ refers to whether there is a uniqueness constraint _Y_ = 1, wherein the same cluster can not be picked twice by accretion or erosion, or not when _Y_ = 0.
8. A number: when preceding the "Run" tag, this indicates population size. When following the "Run" tag, this indicates the number we used to keep track of our simulations.

## Analysis Folders
Analysis folders generally mirror the structure of the associated data folders, including file naming conventions. Additional tags refer to the starting and stopping samples from the files, but generally all samples were used.

Files are .RData files, loaded with the base R load function, and whose contents are primarily "result_storage_subset". This object has the following columns:

1. time: the row number of the results matrix, which can be translated by our convention established above to _10 * time_ time steps.
2. xmin: the fitted smallest cluster size. Note that for the basic MLE, this is the smallest cluster in the population, but can vary when a KS variant is used.
3. xmax: the "fitted" largest cluster size. This should be the largest cluster in the system (but we have unsupported additional functionality that tries to determine the largest cluster using a KS approach.
4. N or clusters: the number of clusters in the system. This is not a fitted quantity, but refers to the row sum of the original results matrix. Note that this is an additional column reported after an update to the pipeline, so analyses generated from code that predate this column are sometimes combined again with the data to generate this column for plots.
5. n_1: the number of clusters of size 1. This is also an additional column reported after an update to the pipeline, so analyses generated from code that predate this column are sometimes combined again with the data to generate this column for plots.
6. exponent: the fitted power-law exponent. This assumes that the power law is a good fit to the system, which we do not guarantee.
7. method_type: two options in the analyses provided. "Basic Truncated MLE" uses the smallest cluster size, the largest cluster size, and the MLE for the exponent. "CSN + ZXX MLE and KS: Lower Bound" follows the [Clauset, Shalizi, and Newman (2009)](https://www.jstor.org/stable/25662336) approach to fitting power laws. (Our implementation is meant to allow for variation in which distributional paramters are fitted via MLE or KS approach, according to [Zhu, Xie, and Xu (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1111/anzs.12162), but this functionality is essentially unused in the main text, which hews extremely closely to the CSN approach.)
8. error: generally FALSE, but occasionally used to report issues captured by try-catch guards during the fitting process. These issues can usually be rectified upon inspection during summarisation (i.e., in a table, plot, interval).
9. kernel: contains pertinent information from tags. Sometimes NA due to coding error (see the change in R's default stringsAsFactors argument).
10. file: contains the associated filename. Sometimes NA, as in kernel above.

## Coding and Modelling Folder
The files in this folder collect the various tools that were necessary to generate the data, analyses, tables, figures, and statistics discussed above and in the main text. We describe them in the rough flow of the pipeline below. We additionally note the usage of a .yaml file to describe our Python environment on the now decommissioned [Viking computing cluster](https://wiki.york.ac.uk/display/RCS/Viking+-+University+of+York+Research+Computing+Cluster).

### Python -- Simulation
Python code can be divided into code defining the process and code that runs the simulation. To assist in understanding our terminology, we internally referred to our original model as Johnson's model after [Johnson et al. (2016)](http://science.sciencemag.org/content/352/6292/1459), in which coalescence and fragmentation was used to model terrorist support networks. Defining the process, we have:

1. JohnsonAggregator.py: various coalescence functions.
2. JohnsonDisaggregator.py: various fragmentation functions.
3. picker.py: functions that characterise the kernels.
4. counting.py: functions that broadly a) detect numbers and positions of elements, b) partition numbers, c) summarise the population.
5. ChineseRestaurant.py: functions that implement various attempts at the Chinese Restaurant Process and related processes.
6. JohnsonOnePopModels.py: various stages of functions that bring together all of the above to implement the full coalescence and fragmentation process. MixedKernelModel is our end result.

The following files run the simulations:

1. 1_Generate_Populations_2.py
2. 1_Generate_Populations_3.py
3. 1_Generate_Populations_4.py
4. 1_Generate_Populations_5.py
5. 1_Generate_Populations_ImperfectComparisons.py
6. 1_Generate_Populations_KernalComparisonPopulations.py

These files all have similar formatting and expect command line arguments indicating which simulations to run and which parameters to set. They then generate directories, start the simulation from a random partition burn-in the simulation, and then run sections of the simulation, writing and then appending to a .mat file using the SciPy.io interface.

### MATLAB -- Sparsification
The following files load the .mat results matrices, combine the results into a single matrix per file, and then saves the single matrix as a sparse matrix.

1. MAT_to_Sparse.m
2. MAT_to_Sparse_NoFlux.m

Both functions serve to convert the .mat files generating in Python into a consistent sparse version of the .mat files, which proved easier to load into R. Differences between the two files reflect a beginnings of an analysis of flux across cluster sizes, and changing directory structure and software syntax.

### R -- Analysis
The following files are used to generate the analyses.

1. run_fitting_embarassingly_parallel.R: workhorse script for analysing .mat files. As noted above, the exact format of the output has changed over time. The script primarily reads in the sparse .mat file, breaks it into rows, and then applies power-law fitting methods to each row before binding the results together and saving. As the name suggests, this script is formatted to facilitate switching between parallel and serial fitting.
2. rexpl.R: Derived from run_fitting_embarassingly_parallel, this script constructs partitions using the equation for the steady-state and then fits them with the same power-law fitting methods. Generates figure 3.
3. Summary_CalculateRunLengths.R: an earlier version of Tables.R, this script takes the sparse .mat, organises them into runs in which the largest cluster is not fragmented, interpolates them slightly, and then converts and summarises the lengths of runs into the cyclicity parameter.
4. Tables.R: As in 3., but set to provide the cyclicities for both of the tables in the text, rather than base model rows.

The following files are used to generate elements of the text (statistics and figures).

1. Figure_CreateSteadyStateExample.R: Generates figures 1 and 2.
2. Figure_CreateSteadyStateExample2.R: Generates figure 4.
3. Figure_CreateSteadyStateExample5.R: Generates figure 5.
4. Figure_CreateSteadyStateExample6.R: Generates figures 6, 7, and 8.
5. Figure_CreateGelShatterExample9.R: Generates figures 9, 10, 11, and 12.
6. Figure_CreateTablePlots.R: Generates figure 13.
7. Figure_SI_LoadAndAverageTargets.R: A supplemental figure generator that lets one inspect PDF and CCDF plots.

### Pipeline
In principle, the pipeline should be straightforward. 

1. Select the kernels one is interested in, 
2. Select the coalescence and fragmentation functions one is interested in, 
3. Write a script to run MixedKernelModel with the desired combination of functions,
4. Save the results as .mat matrices,
5. Sparsify the matrices in MATLAB with something like MAT_to_Sparse.m,
6. Run the fitting and summarising script run_fitting_embarassingly_parallel (with some edits to make sure your specific .mats are identified),
7. Summarise the cyclicity (e.g. with some edits to Tables.R or Summary_CalculateRunLengths.R),
8. Plot the results.

Edits should mostly need to be made due to changes in directory structure, but there have been several changes in standards over the lifetime of this project, so no long term guarantee can be provided.

# Robustness-Coal-Frag
Code accompanying manuscript entitled "Robustness of steady state and stochastic cyclicity in generalized coalescence-fragmentation models".

# Folders
Three folder types. Data folders have raw data (created from Python code and then sparsified via MATLAB). Analysis folders have these data files summarised to relevant statistics (using R). The Code-Model folder contains the relevant scripts that were used to define, run, and analyse the simulations, as well as to organise the information for the tables and figures.

## Data Folders
Folders are organised roughly by date and theme. Each file is a MATLAB .mat file, with those that have the tag "_S" converted in MATLAB from those without as sparse matrices (see "MAT_to_Sparse" family of MATLAB scripts). Those that do not have the tag "_S" were instead created in Python using SciPy.io's savemat function (see "1_Generate" family of Python scripts). Each row of the matrix represents the population after a time step, with column i representing the number of clusters of size i. By convention and to reduce the amount of data stored, data was only stored every 10 time steps. As standard computational practice the simulations were burnt in, usually with 1,000,000 events, and were started from random partitions of the population. To facilitate saving and storing the simulation results in Python during batch jobs, we began saving (time) "slices" of the results matrix. These were then "stitched" together when the MATLAB sparsification script was run.

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

## Coding and Modelling Folder

### Pipeline

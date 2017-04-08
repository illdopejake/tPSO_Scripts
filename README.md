# tPSO_Scripts

Include in this repository is the code use to perform many of the analyses in the paper "Proximity to parental symptom onset is associated with brain network abnormalities in sporadic Alzheimer’s disease"

Most notably, a pipeline exists within this repository that was used to execute the spatial correlation analysis for the structural-functional section of the paper. This readme will contain detailed information on what the pipeline does, and some notes about running it.

# Notes about the spatial correlation pipeline and the contents of the repo

All of the actual code is written in Python and the individual commands in the pipeline, as well as full documentation of their use, can be found in bootstrap_spatial_correlations.py. There are a few other functions that are also borrowed from a different group of python utilities I wrote to interface with NIAK. These can be found in propagation_correlations.py and jake_niak_utils.py

The analysis for the paper involved a large amount of computing power, as 5000 voxelwise analysis permutations and 1227 regression analyses (* 5000 permutations) were performed. Therefore, the pipeline was parallelized and controlled using the Pipeline System for Octave and Matlab: https://github.com/SIMEXP/psom. 

To interface with PSOM, a few things had to be done: 
1) First each python function was wrapped into a command-line executable function. These are the cmd_*.py functions in the repo
2) Matlab functions were generated to take matlab arguments and use them to execute the command-line functions from 1). These are the [function].m functions
3) Finally, a PSOM-friendly pipeline was created to take user inputs and execute the various functions for the spatial correlation pipeline in parallel on the Guillimin super computer, an instance of Compute Canada. That's PIPELINE_spatial_correlation.m. An example wrapper script for the pipeline can be found in example_wrapper_script_PSOM2.m

Any of the individual commands can be run outside of the pipeline, but most of the functions accept outputs from other functions as their inputs, so its probably better to run the whole pipeline at once. 

At this time, the pipeline and commands are very specifically and exclusively written to handle NIAK output, and NIAK _glm.mat files in particular. The script can however be adapted to handle inputs from other rsfMRI pipelines. Please contact me if you're interested.

# Dependencies
The various functions require several nonstandard Python libraries, including numpy, pandas, randomstate, statsmodel and nibabel

In addition, the voxelwise_analysis function uses FSL if the norm argument is not set to False 

Finally, the actual pipeline requires PSOM, NIAK and minc_toolkit

# Methods 

Space limitations precluded a thorough and detailed account of the methods behind the spatial correlation pipeline. Here is a fully detailed explanation of the methods, as if it were appearing in the actual paper:

We ran analyses comparing the pattern of GMV reduction associated with spEYO to resting state networks generated at different resolutions. Briefly, we gathered functional connectivity atlases at nine different resolutions (see paper) and derived average functional connectivity maps for each region within each atlas using a group of healthy young adults. Modified spatial correlation analysis (Buckner et al., 2005; Zeighami et al., 2015) was used to characterize the topographical similarity between each connectivity map and an effect map representing the voxelwise relationship between spEYO and GMV (see paper: Figure 4A for visual aid). Stronger correlations represent greater spatial correspondence between patterns of functional connectivity and structural features related to spEYO, providing evidence that the observed structural features may be partially informed or constrained by the network topography of the brain.  

1227 analyses were performed in total – one for each parcel (brain region) within each of the nine functional atlases. First, the unthresholded voxelwise GMV t-value map (GMV-t) from the morphometric analysis (paper: Figure 1A) was resampled to 3mm isotropic, and then moved to the same resolution as each of the functional brain atlases by averaging GMV-t voxels within each of the parcels in the respective atlas. Next, for each of the 1227 parcels, the following procedures were performed: A covariate-adjusted functional connectivity (FC-t) map seeded from the nth parcel (parceln) was generated using young healthy controls (see above), and Fisher’s z-transformed correlations from this map were converted to t-values. Next, values in the GMV-t and FC-t maps were vectorized so that each vector contained k values – where k is the number of parcels in the atlas. Parcelwise spatial Spearman’s rank correlations between the GMV-t and FC-t vectors were computed, where each parcel represents one “observation”. The spatial correlation compared the effect of tPSO on GMV in each parcel, to the connectivity of that parcel to parceln. Higher positive correlations indicate increasing topographical similarity between the observed pattern of neurodegeneration specifically associated with spEYO and the resting state functional network seeded from parceln. 

Because of factors relating to the resolution and structure of the data, p-values associated with the Spearman's ρ values resulting from the spatial correlations will be meaningless. In addition, the statistics will be impossible to compare across resolutions because the degrees of freedom differ for spatial correlations at different resolutionss. To properly assess the likelihood of the relationship occurring by chance, the analysis was repeated using Monte Carlo simulations over 5000 permutations to derive exact p-values. For each of the 5000 samples, labels for the covariate-adjusted GMV-t volumes were permuted without replacement using the mrg32k3a pseudorandom number generator (L’Ecuyer, 1999), which offers 2^51 parallel streams of 2^71 numbers. A seed was entered to control the stream, and for each successive sample, a “jump” of 2^121 positions was instituted. The number at the next consecutive position was used as the seed for random number generation, which was used to permute the sample. The permuted sample then underwent voxelwise correlation with spEYO to create a newly randomized GMV-t map, and the methods described in the previous paragraph were implemented once again to attain the spatial correlation statistics. Using this method, a null distribution of Spearman's ρ values was generated, respectively. Exact p-values were derived by determining where the observed test value fell within the null distribution. These methods were implemented through in-house scripts (from this repo) using Python 2.7.9. For Monte Carlo simulations, these scripts were wrapped into a pipeline and executed in parallel using PSOM on the Guillimin super computer, an implementation of Compute Canada




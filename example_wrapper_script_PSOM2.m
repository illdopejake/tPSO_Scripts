clear all

%RUN jake_bootstrap_spatial_correlation
%   For explanation of function and inputs, see python scripts by doing one of the following:
%   1. For any of the cmd_* functions, into your shell, type 'python [func_name]'
%   2. Open python or ipython. Import boostrap_spatial_correlation. Type help(func_name) 

%addpath(genpath('/gs/project/gsf-624-aa/quarantaine/niak-boss-0.12.18/'));
addpath(genpath('/gs/project/gsf-624-aa/quarantaine/niak-dev-0.14.0'))
addpath(genpath('/gs/scratch/jvogel44/bsc/scripts/'));

%%%SET INPUT FILES%%%

% spreadsheet with data
files_in.ss = '/gs/scratch/jvogel44/bsc/stopad_updated_csf_fMRI_edited.csv';
% path to subject scans
files_in.subpath = '/gs/scratch/jvogel44/bsc/gm_files_raw_real/';
% path to directory containing other inputs (e.g. atlases and such)
files_in.indir = '/gs/scratch/jvogel44/bsc/indir/';
% path to file to be used for coregistration
files_in.norm_fl = '/gs/scratch/jvogel44/bsc/rMNI152_T1_2mm_brain.nii';
% path to transformation to be used for coregistration
files_in.xfm = '/gs/scratch/jvogel44/bsc/indir/bsc_tfm'
%files_in.mss=

%%%SET PIPELINE OPTIONS%%%

opt.folder_out = '/gs/scratch/jvogel44/bsc/results/bsc_perm_1000_yc_raw_20160908';
opt.label_out = 'raw'; % label for all outputs
opt.rand_seed = 1; % set the seed for randomization, so results can be reproduced
opt.nb_samp = 1000; % number of bootstrap samples
opt.flag_test = false; % Generate all the inputs without executing pipeline

%%%SET INPUTS AND OPTIONS FOR DEFINE BOOTSTRAP SAMPLE%%%

%INPUTS
opt.dbs.subcol = 'scan'; % name of column label with subject ID
opt.dbs.pv = 'expected_age_onset_no_sibs_min_actual'; % column label for predictor variable

%OPTIONS

opt.dbs.resamp = 'permute'; % run permutation testing
opt.dbs.par = 'True'; % run permutation samples in parallel
%opt.dbs.perc = 0.50; 
%opt.dbs.gps = 3;

%%%SET OPTIONS FOR VOXELWISE ANALYSIS%%%

opt.va.outmap = 't';  % generate an output tmap or rmap
opt.va.nonpar = 'False'; % run voxelwise spearmans instead of pearsons correlations
opt.va.par = 'True'; % Run voxelwise analyses in parallel
opt.va.inter = 'False'; % Save intermediate 4D map for voxelwise analysis

%%%SET INPUTS AND OPTIONS FOR SEARCHLIGHT%%%

%INPUTS
opt.scs.contrast= 'avg'; % Contrast to be used from NIAK glm_file
opt.scs.templ_str= 'masked'; % string to help find atlas templates

%OPTIONS
opt.scs.sclstr = 'scale';  % string to help determine scale of a given atlas
opt.scs.eff = 't'; % take information from eff or ttest structures in glm_file
opt.scs.poly = 2; % also run n-order polynomial regression where n=poly

%%%SET OPTIONS FOR ID SIG RESULTS%%%

% The following inputs will save *all* results
opt.isr.perc = 'perc'; 
opt.isr.type = 'r';
opt.isr.thresh = '1'; 
opt.isr.outtype = 'samp';
opt.isr.res = 'all' ;
%opt.isr.mthresh = 0.01
opt.isr.save = 'True';

%%%SET OPTIONS FOR COLLECT RESULTS

opt.cr.thr = 'r';
opt.cr.resamp = '0.05';
opt.cr.test = 'True';
opt.cr.sum = 'False';

%%%SET INPUTS AND OPTIONS FOR CREATE OUTPUT IMAGES

opt.coi.sclstr = 'masked_scale';
opt.coi.outtp = 'resamp';

%%%%%%%%%%%%%%%%%%%%
%%%RUN THE PIEPLINE
%%%%%%%%%%%%%%%%%%%%
%opt.psom.qsub_options = '-A yai-974-aa -q lm -l nodes=1:ppn=2 -l pmem=5700 -l walltime=2:00:00';
%opt.psom.max_queued = 15;
[pipeline,opt] = jake_bootstrap_spatial_correlation(files_in,opt);

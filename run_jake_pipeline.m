clear all

%RUN jake_bootstrap_spatial_correlation
%   For explanation of function and inputs, see python scripts by doing one of the following:
%   1. For any of the cmd_* functions, into your shell, type 'python [func_name]'
%   2. Open python or ipython. Import boostrap_spatial_correlation. Type help(func_name) 

addpath(genpath('/gs/project/gsf-624-aa/quarantaine/niak-boss-0.12.18/'));
addpath(genpath('/gs/scratch/jvogel44/bsc/scripts/'));

%%%SET INPUTS%%%

files_in.ss = '/gs/scratch/jvogel44/bsc/stopad_updated_csf_fMRI_edited.csv';
files_in.subcol = 'scan';
files_in.subpath = '/gs/scratch/jvogel44/bsc/gm_files/';
files_in.pv = 'expected_age_onset_no_sibs_min_actual';
files_in.indir = '/gs/scratch/jvogel44/bsc/indir/';
files_in.contrast= 'avg';
files_in.templ_str= 'masked';
files_in.norm_fl = '/gs/scratch/jvogel44/bsc/rMNI152_T1_2mm_brain.nii';
%files_in.mss=

%%%SET PIPELINE OPTIONS%%%

opt.folder_out = '/gs/scratch/jvogel44/bsc/results';
opt.label_out = 'pipetest';
opt.nb_samp = 5;
opt.flag_test = true;

%%%SET OPTIONS FOR DEFINE BOOTSTRAP SAMPLE%%%

opt.dbs.perc = 0.50;
opt.dbs.gps = 3;
opt.dbs.par = 'True';

%%%SET OPTIONS FOR VOXELWISE ANALYSIS%%%

opt.va.outmap = 't';
opt.va.nonpar = 'False';
opt.va.par = 'True';
opt.va.inter = 'False';

%%%SET OPTIONS FOR SEARCHLIGHT%%%

opt.scs.sclstr = 'scale'; 
opt.scs.eff = 't';
opt.scs.poly = 2;

%%%SET OPTIONS FOR ID SIG RESULTS%%%

opt.isr.perc = 'top';
opt.isr.type = 'r';
%opt.isr.thresh = 
opt.isr.outtype = 'samp';
opt.isr.res = 'all' ;
%opt.isr.mthresh = 
opt.isr.save = 'True';

%%%%%%%%%%%%%%%%%%%%
%%%RUN THE PIEPLINE
%%%%%%%%%%%%%%%%%%%%
opt.psom.qsub_options = '-A yai-974-aa -q sw -l nodes=1:ppn=2 -l walltime=1:00:00';
opt.psom.max_queued = 15;
[pipeline,opt] = jake_bootstrap_spatial_correlation(files_in,opt);

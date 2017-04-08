function [pipe,opt] = jake_bootstrap_spatial_correlation(files_in,opt)
niak_gb_vars
warning on backtrace
%UNTITLED Summary of this function goes here
%   For explanation of function and inputs, see python scripts by doing one of the following:
%   1. For any of the cmd_* functions, into your shell, type 'python [func_name]'
%   2. Open python or ipython. Import boostrap_spatial_correlation. Type help(func_name) 


%% set up default arguments

if ~exist('files_in','var')||~exist('opt','var')
    error('Input arguments required.')
end

% files_in
files_in = psom_struct_defaults(files_in,...
           { 'ss' , 'subpath', 'indir', 'norm_fl'         , 'xfm'            , 'mss'            },...
           { NaN   , NaN     ,  NaN   , 'gb_niak_omitted' , 'gb_niak_omitted', 'gb_niak_omitted'});
% options
opt = psom_struct_defaults(opt,...
           { 'folder_out' , 'nb_samp' , 'label_out' , 'rand_seed' , 'flag_test', 'psom'    , 'dbs'     , 'va'     , 'scs'    , 'isr'    , 'cr'     , 'coi'   },...
           { NaN          , NaN       ,  'outpt'    ,  false      ,  false     ,  struct() ,  struct() , struct() , struct() , struct() , struct() , struct()}); 


pipe = struct;


%% options for the pipeline

% Psom options
opt.psom = psom_struct_defaults(opt.psom,...
           { 'path_logs'                     },...
           { [opt.folder_out filesep 'logs'] });

% dbs options
opt.dbs = psom_struct_defaults(opt.dbs,...
            { 'subcol' , 'pv' , 'resamp' , 'par'    , 'perc' , 'gps' },...
            { NaN      , NaN  ,  'False' , 'False'  ,  0.5   ,  3    });

opt.va = psom_struct_defaults(opt.va,...
            {'outmap' , 'nonpar' , 'par'  , 'inter' , },...
            { 't'     ,  'False' , 'True' , 'False' , });

opt.scs = psom_struct_defaults(opt.scs,...
            {'contrast' , 'templ_str' , 'sclstr' , 'eff' , 'poly'},...
            { NaN       , NaN         , 'scale'  ,  't'  ,   1   });

opt.isr = psom_struct_defaults(opt.isr,...
            {'perc' , 'type' , 'thresh'           , 'outtype' , 'res' , 'mthresh'          , 'save' },...
            { 'top' , 'r'    ,  'gb_niak_omitted' ,  'samp'   , 'all' ,  'gb_niak_omitted' , 'False'});

opt.cr = psom_struct_defaults(opt.cr,...
            {'thr' , 'resamp' , 'test' , 'sum'  },...
            { 'r'  ,  'False' ,  false , 'False'});

opt.coi = psom_struct_defaults(opt.coi,...
            {'sclstr' , 'outtp'  },...
            { NaN     , 'resamp'});


%% building the pipeline

cr_in = files_in;

for taskid = 1:opt.nb_samp

    % set up the options for function dbs
    opt.dbs.tid = taskid;
    opt.dbs.outdir = opt.folder_out;
    opt.dbs.flag_test = opt.flag_test;
    opt.dbs.rand = opt.rand_seed;
    dbs_name = sprintf('define_boostrap_samples_%d',taskid);
    dbs_in = files_in;
    dbs_opt = opt.dbs;
    dbs_out = [opt.folder_out filesep sprintf('dbc_out%d.xls',taskid)];
    
    pipe = psom_add_job(pipe, dbs_name, 'define_bootstrap_samples', dbs_in, dbs_out, dbs_opt);
    
    % set up options for va
    opt.va.tid = taskid;
    opt.va.outdir = opt.folder_out;
    opt.va.outstr = opt.label_out;
    opt.va.flag_test = opt.flag_test;
    va_name = sprintf('voxelwise_analysis_%d',taskid);
    va_in = files_in;
    va_in.parin = dbs_out;
    va_opt = opt.va;
    va_out = [opt.folder_out filesep sprintf('%s%d.nii',opt.label_out,taskid)]; 

    pipe = psom_add_job(pipe,va_name, 'voxelwise_analysis', va_in, va_out, va_opt);

    % set up the options for function scs 
    opt.scs.tid = taskid;
    opt.scs.outdir = opt.folder_out;
    opt.scs.outstr = opt.label_out;
    opt.scs.flag_test = opt.flag_test;
    scs_name = sprintf('searchlight_%d',taskid);
    scs_in = files_in;
    scs_in.cov_img = va_out;
    scs_opt = opt.scs;
    %this line is a placeholder for psom. Also needs to be fixed so the 7 is a real number from inputs
    scs_out = [opt.folder_out filesep sprintf('%s_scl7_res%d.xls',scs_opt.outstr,taskid)];

    pipe = psom_add_job(pipe,scs_name, 'searchlight', scs_in, scs_out, scs_opt);

    % set up options for function isr
    opt.isr.tid = taskid;
    opt.isr.outdir = opt.folder_out;
    opt.isr.outstr = opt.label_out;
    opt.isr.flag_test = opt.flag_test;
    isr_name = sprintf('id_significant_results_%d',taskid);
    isr_in = files_in;
    %this line is a placeholder for psom
    isr_in.dud = scs_out;
    isr_opt = opt.isr;
    isr_out = struct();
    %isr_out.outdir = opt.folder_out;
    isr_out.samp = [opt.folder_out filesep sprintf('%s_thrres%d.xls',opt.label_out,taskid)];
    if isr_opt.outtype ~= 'samp'
        if size(files_in.mss)(1) > 0
            isr_out.mast = files_in.mss;
        else 
            isr_out.mast = [opt.folder_out filesep sprintf('master_ss.csv')];
        end
    end
    pipe = psom_add_job(pipe,isr_name,'id_sig_results',isr_in, isr_out, isr_opt);

    % load outputs into job cr files_in structure (actually just placeholders for psom dependencies)
    cr_in.(int2str(taskid)) = isr_out.samp;
end

% if necessary, create test (original) data and results 
if opt.cr.test

    ntid = opt.nb_samp + 1

    opt.dbs.tid = ntid;
    opt.dbs.outdir = opt.folder_out;
    opt.dbs.flag_test = opt.flag_test;
    opt.dbs.rand = opt.rand_seed;
    dbs_name = 'define_boostrap_samples_test';
    dbs_in = files_in;
    dbs_opt = opt.dbs;
    dbs_opt.resamp = 'False';
    dbs_out = [opt.folder_out filesep sprintf('dbc_out%d.xls',ntid)];
    pipe = psom_add_job(pipe, dbs_name, 'define_bootstrap_samples', dbs_in, dbs_out, dbs_opt);

    opt.va.tid = ntid;
    opt.va.outdir = opt.folder_out;
    opt.va.outstr = 'TEST';
    opt.va.flag_test = opt.flag_test;
    va_name = 'voxelwise_analysis_test';
    va_in = files_in;
    va_in.parin = dbs_out;
    va_opt = opt.va;
    va_out = [opt.folder_out filesep sprintf('%s%d.nii',va_opt.outstr,ntid)];     
    pipe = psom_add_job(pipe,va_name, 'voxelwise_analysis', va_in, va_out, va_opt);

    % set up the options for function scs
    opt.scs.tid = ntid;
    opt.scs.outdir = opt.folder_out;
    opt.scs.outstr = 'TEST';
    opt.scs.flag_test = opt.flag_test;
    scs_name = 'searchlight_test';
    scs_in = files_in;
    scs_in.cov_img = va_out;
    scs_opt = opt.scs;
    %this line is a placeholder for psom. Also needs to be fixed so the 7 is a real number from inputs
    scs_out = [opt.folder_out filesep sprintf('%s_scl7_res%d.xls',scs_opt.outstr,ntid)];
    pipe = psom_add_job(pipe,scs_name, 'searchlight', scs_in, scs_out, scs_opt);

   % set up options for function isr
    opt.isr.tid = ntid;
    %opt.isr.outdir = opt.folder_out;
    opt.isr.outstr = 'TEST';
    opt.isr.flag_test = opt.flag_test;
    isr_name ='id_significant_results_test';
    isr_in = files_in;
    %this line is a placeholder for psom
    isr_in.dud = scs_out;
    isr_opt = opt.isr;
    isr_out = struct();
    %isr_out.outdir = opt.folder_out;
    isr_opt.outtype = 'samp';
    isr_out.samp = [opt.folder_out filesep sprintf('%s_thrres%d.xls',isr_opt.outstr,ntid)];
    
    pipe = psom_add_job(pipe,isr_name,'id_sig_results',isr_in, isr_out, isr_opt);

else
    opt.cr.test = 'False';
end

% set up options for function cr

cr_name = 'collect_results';
cr_in.ssdir = opt.folder_out;
if ~(strncmp(opt.cr.test,'False',4))
    cr_in.test = isr_out.samp;
end
opt.cr.ssstr = opt.label_out;
opt.cr.ssext = 'xls';
opt.cr.outdir = opt.folder_out;
opt.cr.outfl = opt.label_out;
opt.cr.flag_test = opt.flag_test;
cr_opt = opt.cr;
cr_out = struct();
cr_out.pearson = [opt.folder_out filesep sprintf('%s_finalres%s.xls',cr_opt.outfl,'r')];
cr_out.spearman = [opt.folder_out filesep sprintf('%s_finalres%s.xls',cr_opt.outfl,'rho')];
if opt.scs.poly > 1
    cr_out.poly = [opt.folder_out filesep sprintf('%s_finalres%s.xls',cr_opt.outfl,'poly')]
end
pipe = psom_add_job(pipe,cr_name,'collect_results',cr_in, cr_out, cr_opt);


% set up options for functions coi

coi_in = files_in;
opt.coi.outdir = opt.folder_out;
opt.coi.outstr = opt.label_out;
opt.coi.flag_test = opt.flag_test;
coi_opt = opt.coi;

coi_name = 'create_output_images_pearson';
coi_in.resdf = cr_out.pearson;
coi_opt.input = 'r';
%FOR TESTING! THIS NEEDS TO CHANGE!
coi_out = [opt.folder_out filesep sprintf('%s_resample_p_%s7.nii.gz',coi_opt.outstr,coi_opt.input)];
pipe = psom_add_job(pipe,coi_name,'create_output_images',coi_in,coi_out,coi_opt);

coi_name = 'create_output_images_spearman';
coi_in.resdf = cr_out.spearman;
coi_opt.input = 'rho';
%FOR TESTING! THIS NEEDS TO CHANGE!
coi_out = [opt.folder_out filesep sprintf('%s_resample_p_%s7.nii.gz',coi_opt.outstr,coi_opt.input)];
pipe = psom_add_job(pipe,coi_name,'create_output_images',coi_in,coi_out,coi_opt);

if opt.scs.poly > 1
    coi_name = 'create_output_images_poly';
    coi_in.resdf = cr_out.poly;
    coi_opt.input = 'poly';
    %FOR TESTING! THIS NEEDS TO CHANGE!
    coi_out = [opt.folder_out filesep sprintf('%s_resample_p_%s7.nii.gz',coi_opt.outstr,coi_opt.input)];
    pipe = psom_add_job(pipe,coi_name,'create_output_images',coi_in,coi_out,coi_opt);
end

%% run the pipeline

if ~opt.flag_test
    opt.psom.qsub_options = '-A gsf-624-aa -q sw -l nodes=1:ppn=12 -l pmem=2700m -l walltime=12:00:00';
    opt.psom.max_queued = 50;
    psom_run_pipeline(pipe,opt.psom)
end



function [pipe,opt] = jake_bootstrap_spatial_correlation(files_in,opt)
niak_gb_vars

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
           { 'ss' , 'subcol', 'subpath', 'pv', 'indir', 'contrast', 'templ_str', 'norm_fl',         'mss'            },...
           { NaN    , NaN   ,  NaN     ,  NaN,  NaN   ,  NaN      ,  NaN       , 'gb_niak_omitted', 'gb_niak_omitted'});
% options
opt = psom_struct_defaults(opt,...
           { 'folder_out' , 'nb_samp' , 'label_out' , 'flag_test'},...
           { NaN          , NaN       ,  'outpt'    ,  false    });


pipe = struct;


%% options for the pipeline

% Psom options
opt.psom = psom_struct_defaults(opt.psom,...
           { 'path_logs'                     },...
           { [opt.folder_out filesep 'logs'] });

% dbs options
opt.dbs = psom_struct_defaults(opt.dbs,...
            { 'perc' , 'gps' , 'par' , },...
            { 0.5    , 3     , 'True'  });

opt.va = psom_struct_defaults(opt.va,...
            {'outmap' , 'nonpar' , 'par'  , 'inter' , },...
            { 't'     ,  'False' , 'True' , 'False' , });

opt.scs = psom_struct_defaults(opt.scs,...
            {'sclstr' , 'eff' , 'poly'},...
            {'scale'  ,  't'  ,   1   });

opt.isr = psom_struct_defaults(opt.isr,...
            {'perc' , 'type' , 'thresh'           , 'outtype' , 'res' , 'mthresh'          , 'save'},...
            { 'top' , 'r'    ,  'gb_niak_omitted' ,  'samp'   , 'all' ,  'gb_niak_omitted' , 'False'});
      
%% building the pipeline

for taskid = 1:opt.nb_samp

    % set up the options for function dbs
    opt.dbs.tid = taskid;
    opt.dbs.outdir = opt.folder_out;
    dbs_name = sprintf('dbs_%d',taskid);
    dbs_in = files_in;
    dbs_opt = opt.dbs;
    dbs_out = [opt.folder_out filesep sprintf('dbc_out%d.xls',taskid)];
    
    pipe = psom_add_job(pipe, dbs_name, 'define_bootstrap_sample', dbs_in, dbs_out, dbs_opt);
    
    % set up options for va
    opt.va.tid = taskid;
    opt.va.outdir = opt.folder_out;
    opt.va.outstr = opt.label_out;
    va_name = sprintf('va_%d',taskid);
    va_in = files_in;
    va_in.parin = dbs_out;
    va_opt = opt.va;
    va_out = [opt.folder_out filesep sprintf('%s%d.nii',opt.label_out,taskid)]; 

    pipe = psom_add_job(pipe,va_name, 'voxelwise_analysis', va_in, va_out, va_opt);

    % set up the options for function scs 
    opt.scs.tid = taskid;
    opt.scs.outdir = opt.folder_out;
    opt.scs.outstr = opt.label_out;
    scs_name = sprintf('scs_%d',taskid);
    scs_in = files_in;
    scs_in.cov_img = va_out;
    scs_opt = opt.scs;
    scs_out = scs_opt.outstr;

    pipe = psom_add_job(pipe,va_name, 'searchlight', scs_in, scs_out, scs_opt);

    % set up options for function isr
    opt.isr.tid = taskid;
    opt.isr.outdir = opt.folder_out;
    isr_name = sprintf('isr_%d',taskid);
    isr_in = files_in;
    isr_in.outstr = scs_out;
    isr_opt = opt.isr;
    if isr.opt.outtype ~= 'mast'
        isr_out.samp = [opt.folder_out filesep sprintf('%s%s.csv',opt.label_out,taskid)];
    end
    if isr.opt.outtype ~= 'samp'
        if size(files_in.mss)(1) > 0
            isr.out.mast = files_in.mss;
        else 
            isr.out.mast = [opt.folder_out filesep sprintf('master_ss.csv')];
        end
    end

    pipe = psom_add_job(pipe,va_name,'id_sig_results',isr_in, isr_opt, isr_out);
end

%% run the pipeline

if ~opt.flag_test
    psom_run_pipeline(pipe,opt.psom)
end



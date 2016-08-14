function [files_in,files_out,opt] = define_bootstrap_samples(files_in,files_out,opt)
warning on backtrace
%DEFINE BOOTSTRAP SAMPLE
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_define_bootstrap_sample.py
%	2. Open python or ipython, import bootstrap_spatial_correlation,
%	type help(define_bootstrap_sample)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';
if ~opt.flag_test
    system(sprintf('python %s/cmd_define_boostrap_sample.py %s %s %s %s %s -perc %0.2f -gps %d -par %s -tid %d',...
				    script_pth, files_in.ss, opt.subcol, files_in.subpath, opt.pv, opt.outdir,...
				    opt.perc, opt.gps, opt.par, opt.tid ))
end

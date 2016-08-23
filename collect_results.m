function [files_in,files_out,opt] = collect_results(files_in,files_out,opt)
warning on backtrace
%COLLECT RESULTS
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_collect_results.py
%	2. Open python or ipython, import bootstrap_spatial_correlation,
%	type help(collect_results)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';
if ~opt.flag_test
    system(sprintf('python %s/cmd_define_boostrap_sample.py %s %s %s %s -outfl %s -thr %s -resamp %s -test %s -sum %s',...
				    script_pth, files_in.ssdir, opt.ssstr, opt.ssext, opt.outdir,...
				    opt.outfl, opt.thr, opt.resamp, opt.test, opt.sum))
end

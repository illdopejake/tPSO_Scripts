function [files_in,files_out,opt] = searchlight(files_in,files_out,opt)
%Searchlight
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_searchlight.py
%	2. Open python or ipython, import bootstrap_spatial_correlation,
%	type help(spatial_correlation_searchlight_from_NIAK_GLMs)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';
if ~opt.flag_test
    system(sprintf('python %s/cmd_searchlight.py %s %s %s %s %s %s -sclstr %s -eff %s -poly %d -norm %s -tid %d',...
				    script_pth, files_in.indir, files_in.cov_img, opt.outstr, opt.outdir, files_in.contrast, files_in.templ_str,...
				    opt.sclstr, opt.eff, opt.poly, files_in.norm_fl, opt.tid))
end

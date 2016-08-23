function [files_in,files_out,opt] = create_output_images(files_in,files_out,opt)
warning on backtrace
%COLLECT RESULTS
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_create_output_images.py
%	2. Open python or ipython, import bootstrap_spatial_correlations,
%	type help(create_output_images)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';
if ~opt.flag_test
    system(sprintf('python %s/cmd_define_boostrap_sample.py %s %s %s %s %s %s -outtp %s -par %s',...
				    script_pth, 'jnk', files_in.indir, opt.sclstr, opt.outdir, opt.outstr, opt.input,...
				    opt.outtp, files_in.resdf))
end

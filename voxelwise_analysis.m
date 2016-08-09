function [files_in,files_out,opt] = voxelwise_analysis(files_in,files_out,opt)
%DEFINE BOOTSTRAP SAMPLE
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_voxelwise_analysis.py
%	2. Open python or ipython, import bootstrap_spatial_correlation,
%	type help(voxelwise_analysis)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';
if ~opt.flag_test
    system(sprintf('python %s/cmd_voxelwise_analysis.py %s %s %s -out %s -nonpar %s -par %s -tid %d -inter %s',...
				    script_pth, files_in.parin, opt.outstr, opt.outdir,...
				    opt.outmap, opt.nonpar, opt.par, opt.tid, opt.inter))
end
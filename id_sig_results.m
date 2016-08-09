function [files_in,files_out,opt] = id_sig_results(files_in,files_out,opt)
%id_sig_results
%   For detailed explanation of function and inputs/options, you can:
%	1. Into your shell, type python cmd_id_sig_results.py
%	2. Open python or ipython, import bootstrap_spatial_correlation,
%	type help(id_sig_results)

script_pth = '/gs/scratch/jvogel44/bsc/scripts';

if opt.thresh
    if ~ischar(opt.thresh)
	    opt.thresh = mat2str(opt.thresh)
	end
end

if opt.mthresh
    if ~ischar(opt.mthresh)
	    opt.mthresh = mat2str(opt.mthresh)
    end
end

subcmd = sprintf('python %s/cmd_id_sig_results.py',script_pth);
cmd = sprintf(' %s %s -perc %s -tp %s -thr %s -outtp %s -res %s -mss %s -mthr %s -save %s -tid %d',...
				files_in.outstr, opt.outdir,...
				opt.perc, opt.type, opt.thresh, opt.outtype, opt.res, files_in.mss, opt.mthresh, opt.save, opt.tid);
if ~opt.flag_test
    system(sprintf('%s%s',subcmd, cmd));
end
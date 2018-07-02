% EVAL_VOL Evaluate volume estimation accuracy
%
% Usage
%    vol_perf = eval_vol(vol_true, vol_est);
%
% Input
%    vol_true: The true volumes in the form of an L-by-L-by-L-by-K array.
%    vol_est: The estimated volumes in the same form.
%
% Output
%    vol_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the volumes in a vector of size K.
%          - corr: The correlations of the volumes in a vector of size K.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol_perf = eval_vol(vol_true, vol_est)
    if ndims(vol_true) ~= ndims(vol_est) || ...
        any(size(vol_true) ~= size(vol_est))

        error('Volume input arrays must be the same shape.');
    end

    vol_perf = struct();

    err = anorm(vol_true-vol_est, 1:3);
    norm_true = anorm(vol_true, 1:3);

    err = permute(err, [4:ndims(err) 1:3]);
    norm_true = permute(norm_true, [4:ndims(err) 1:3]);

    rel_err = err./norm_true;
    corr = acorr(vol_true, vol_est, 1:3);

    corr = permute(corr, [4:ndims(corr) 1:3]);

    vol_perf.err = err;
    vol_perf.rel_err = rel_err;
    vol_perf.corr = corr;
end

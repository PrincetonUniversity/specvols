% EVAL_VOLMAT Evaluate volume matrix estimation accuracy
%
% Usage
%    volmat_perf = eval_volmat(volmat_true, volmat_est);
%
% Input
%    volmat_true: The true volume matrices in the form of an
%       L-by-L-by-L-by-L-by-L-by-L-by-K array.
%    volmat_est: The estimated volume matrices in the same form.
%
% Output
%    volmat_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the volume matrices in a vector of
%            size K.
%          - corr: The correlations of the volume matrices in a vector of size
%            K.
%
% See also
%    eval_vol

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function volmat_perf = eval_volmat(volmat_true, volmat_est)
    if ndims(volmat_true) ~= ndims(volmat_est) || ...
        any(size(volmat_true) ~= size(volmat_est))

        error('Volume matrices must be the same shape.');
    end

    volmat_perf = struct();

    rel_err = anorm(volmat_true-volmat_est, 1:6)./anorm(volmat_true, 1:6);
    corr = acorr(volmat_true, volmat_est, 1:6);

    volmat_perf.rel_err = permute(rel_err, [7:ndims(rel_err) 1:6]);
    volmat_perf.corr = permute(corr, [7:ndims(corr) 1:6]);
end

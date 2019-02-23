% SRC_PARTITION_MEAN_BACKWARD Apply adjoint mappings to a partition
%
% Usage
%    mean_b_coeffs = src_partition_mean_backward(src, partition_idx, ...
%       basis, mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    partition_idx: A partition in the form of a cell array obtained from
%       `create_partition` containing n_folds folds.
%    basis: A basis object used for representing the volumes.
%    mean_est_opt: A struct of mean estimation options. For more details, see
%       `src_mean_backward`.
%
% Output
%    mean_b_coeffs: The adjoint mapping applied to the different subsets of
%       the source corresponding to folds of the partition in the form of a
%       basis.count-by-n_folds array. For more details, see
%       `src_mean_backward`.
%
% See also
%    src_partition_mean_kernel, param_search_mean, create_partition

function mean_b_coeffs = src_partition_mean_backward(src, partition_idx, ...
    basis, mean_est_opt)

    if nargin < 4 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    n_folds = numel(partition_idx);

    partition_srcs = cellfun(@(idx)(subset_src(src, idx)), partition_idx, ...
        'uniformoutput', false);

    mean_b_fun = @(src_sub)(src_mean_backward(src_sub, basis, mean_est_opt));

    mean_b_coeffs = cellfun(mean_b_fun, partition_srcs, ...
        'uniformoutput', false);

    mean_b_coeffs = cat(2, mean_b_coeffs{:});
end

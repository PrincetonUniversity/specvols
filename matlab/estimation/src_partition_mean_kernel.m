% SRC_PARTITION_MEAN_KERNEL Calculate mean estimation kernels for a partition
%
% Usage
%    mean_kernels_f = src_partition_mean_kernel(src, partition_idx, ...
%       mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    partition_idx: A partition in the form of a cell array obtained from
%       `create_partition` containing n_folds folds.
%    mean_est_opt: A struct of mean estimation options. For more details, see
%       `src_mean_kernel`.
%
% Output
%    mean_kernels_f: A 2*L-by-2*L-by-2*L-by-n_folds array containing the
%       non-centered Fourier transforms of the mean least-squares estimator
%       kernels. One kernel is calculated for each fold in the partition. For
%       more detials, see `src_mean_kernel`.
%
% See also
%    src_partition_mean_backward, param_search_mean, create_partition

function mean_kernels_f = src_partition_mean_kernel(src, partition_idx, ...
    mean_est_opt)

    if nargin < 3 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    partition_srcs = cellfun(@(idx)(subset_src(src, idx)), partition_idx, ...
        'uniformoutput', false);

    mean_kernel_fun = @(src_sub)(src_mean_kernel(src_sub, mean_est_opt));

    mean_kernels_f = cellfun(mean_kernel_fun, partition_srcs, ...
        'uniformoutput', false);

    mean_kernels_f = cat(4, mean_kernels_f{:});
end

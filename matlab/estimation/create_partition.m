% CREATE_PARTITION Create partition for cross validation
%
% Usage
%    partition_idx = create_partition(n_images, n_folds);
%
% Input
%    n_images: The number of images to partition.
%    n_folds: The number of folds in the partition.
%
% Output
%    partition_idx: A cell array of indices forming the partition.

function partition_idx = create_partition(n_images, n_folds)
    partition_idx = cell(n_folds, 1);

    n_per_fold = floor(n_images/n_folds);
    n_large_folds = n_images - n_per_fold * n_folds;

    current_ind = 1;

    for k = 1:n_folds
        last_ind = current_ind + n_per_fold - 1;

        if k <= n_large_folds
            last_ind = last_ind + 1;
        end

        partition_idx{k} = current_ind:last_ind;

        current_ind = last_ind + 1;
    end
end

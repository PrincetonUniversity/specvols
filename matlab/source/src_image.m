% SRC_IMAGE Extract images from source
%
% Usage
%    im = src_image(src, start, num);
%
% Input
%    src: Source object.
%    start: First index of the images to extract.
%    num: The number of images to extract.
%
% Output
%    im: The source images from start to start+num-1.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = src_image(src, s, n)
    if ~isfield(src, 'type')
        error('''src'' must have ''type'' field');
    end

    if src.type == src_type_sim()
        im = sim_image(src.sim, s, n);
    elseif src.type == src_type_file()
        im = zeros([src.L*ones(1, 2) n], src.precision);

        all_idx = s:s+n-1;

        file_idx = src.file_idx(all_idx);
        image_idx = src.image_idx(all_idx);

        unique_file_idx = unique(file_idx);
        file_image_idx = arrayfun( ...
            @(file_id)(image_idx(file_idx == file_id)), ...
            unique_file_idx, 'uniformoutput', false);

        for file_id = unique_file_idx(:)'
            file_name = src.file_names{file_id};
            file_s = min(file_image_idx{file_id});
            file_n = max(file_image_idx{file_id}) - file_s + 1;

            file_im = load_image(file_name, file_s, file_n);
            file_im = file_im(:,:,file_image_idx{file_id} - file_s + 1);

            im(:,:,file_idx == file_id) = file_im;
        end
    elseif src.type == src_type_downsampled()
        im_original = src_image(src.original_src, s, n);

        im = im_downsample(im_original, src.L);
    elseif src.type == src_type_array()
        im = src.images(:,:,s:s+n-1);
    elseif src.type == src_type_subset()
        im = zeros([src.L*ones(1, 2) n], src.precision);

        % TODO: Speed this up by getting images in contiguous blocks.

        for k = 1:n
            im(:,:,k) = src_image(src.original_src, src.subset_idx(s+k-1), 1);
        end
    elseif src.type == src_type_filtered()
        im = src_image(src.original_src, s, n);

        filter_idx = src.filter_idx(s:s+n-1);

        unique_filter_idx = unique(filter_idx);

        for k = unique_filter_idx
            mask = find(filter_idx == k);

            im(:,:,mask) = im_filter(im, src.filters(k));
        end
    elseif src.type == src_type_normalized_bg()
        im = src_image(src.original_src, s, n);

        im = im_normalize_bg(im, src.bg_radius, src.do_ramp);
    end
end

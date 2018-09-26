% SRC_VOLS_WT_BACKWARD Apply weighted adjoint mapping to source with
% weights
%
% Usage
%    vols_wt_b_coeff = src_vols_wt_backward(src, basis, wts, mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    basis: A basis object used for representing the volumes.
%    wts: a vector of weights to apply to the respective images when
%    finding the mean of backprojected images
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%
% Output
%    vols_wt_b_coeff: A set of volumes, each of which is a copy of the adjoint
%       mapping applied to the images, averaged with weights given by the
%       respective columns of wts, over the whole dataset and expressed as
%       coefficients of "basis".  A basis.count x r matrix
%
% See also
%    im_backward_wt, src_vols_wt_kermat

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function vols_wt_b_coeff = src_vols_wt_backward(src, basis, wts, vols_wt_est_opt)
    if nargin < 3 || isempty(vols_wt_est_opt)
        vols_wt_est_opt = struct();
    end

    vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    params = src.params;

    L = src.L;
    n = src.n;
    r = size(wts,1);
    
    if n ~= size(params.rots, 3)
        error('Number of images in source and parameters do not agree.');
    end

    batch_ct = ceil(n/vols_wt_est_opt.batch_size);

    vols_wt_b = zeros([L L L r], vols_wt_est_opt.precision);

    tic
    disp('coef_b')
    for batch = 1:batch_ct
        disp(['batch ' num2str(batch) ' t = ' num2str(toc)]);
        batch_s = (batch-1)*vols_wt_est_opt.batch_size+1;
        batch_n = min(batch*vols_wt_est_opt.batch_size, n)-batch_s+1;

        batch_idx = batch_s:(batch_s+batch_n-1);
        
        local_slices = zeros([L L L batch_n]);
        
        parfor k = 1:batch_n
            idx = k+batch_s-1;
%             local_slices = im_backward(src,src_image(src,k+batch_s-1,1), k+batch_s-1);
            im = src_image(src,idx,1);
            scaled_im = params.amplitudes(idx) * im;
            % I should incorporated the shifts here directly, since it's
            % just a phase shift in Fourier...
            shifted_im = im_translate(scaled_im,-params.offsets(:,idx));
            filtered_im = im_filter(shifted_im,params.filters(params.filter_idx(idx)));
            %local_slices(:,:,:,k) = imbackproject(filtered_im,params.rots(:,:,idx);
            pts_rot = rotated_grids(L,params.rots(:,:,idx));
            pts_rot = reshape(pts_rot,[3 L^2]);
            im_f = 1/L^2 * centered_fft2(filtered_im);
            if mod(L, 2) == 0
                im_f(1,:,:) = 0;
                im_f(:,1,:) = 0;
            end

            im_f_flat = reshape(im_f, [L^2 1]);
            nufft_opt = struct();
            nufft_opt.num_threads = 1;
            vol = 1/L*anufft3(im_f_flat, pts_rot, [L L L],nufft_opt);
            local_slices(:,:,:,k) = real(vol);
        end

        %note: I can do the below using a single bsxfun instead of a for
        %loop and a bsxfun, but I like how this is similar to the kernel
        %calculation this way (and that one isn't at present implemented as
        %a bsxfun, though it could be with an extra step).
        %It seems like I should be able to get rid of the second permute by
        %just taking the sum along the first coordinate, but then I'd still
        %have to permute to get right of the first dimension being a
        %singleton!  Life is pain.
        for l = 1:r
            wted_sum = sum(permute(bsxfun(@times,wts(l,batch_idx)',permute(local_slices,[4 1 2 3])),[2 3 4 1]),4);
            vols_wt_b(:,:,:,l) = vols_wt_b(:,:,:,l) + cast(wted_sum,vols_wt_est_opt.precision);
        end
    end

    for l = 1:r
        vols_wt_b_coeff(:,l) = basis_evaluate_t(basis,vols_wt_b(:,:,:,l));
    end
%     vols_wt_b_coeff = vols_wt_b;
end

% apply_vols_wt_kernel Applies the matrix of weighted mean kernels
% represented by sums of convolutions
%
% Usage
%    vol_coeff = vol_coeff = apply_vols_wt_kernel(vol_coeff, kermat_f, ...
%       basis, mean_est_opt)
%
% Input
%    vol_coeff: The volume to be convolved, stored in the basis coefficients.
%    kermat_f: A "matrix" of non-centered FTs of convolution kernels, each
%       representing a weighted-mean projection/backprojection operator (as
%       obtained from `src_vols_wt_kernel`). Note that each kernel is applied
%       separately and the sums taken.
%    basis: A basis object corresponding to the basis used to store
%       `vol_coeff`.
%    mean_est_opt: An options structure. Currently no options are used.
%
% Output
%    vols_coeff: The result of evaluating `vols_coeff` in the given basis,
%       convolving with the kernel given by `kernel_f`, and backprojecting into
%       the basis.

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function vols_coeff = apply_vols_wt_kernel(vols_coeff, kermat_f, basis, ...
    vols_wt_est_opt)

    unflatten_flag = 0;

    if(size(vols_coeff,1) ~= basis.count)
        r = size(vols_coeff,1) / basis.count;
        if( rem(r,1) == 0)
            unflatten_flag = 1;
            vols_coeff = reshape(vols_coeff,[basis.count r]);
        else
            error('Dimensions must be basis.count times number of vols');
        end
    end

    vols_in = basis_evaluate(basis, vols_coeff);

    vols_out = zeros(size(vols_in), class(vols_in));

    for k = 1:size(kermat_f,4)
        vols_out(:,:,:,k) = sum(conv_vols_wt(vols_in, permute(kermat_f(:,:,:,k,:),[1 2 3 5 4])), 4);
    end

    vols_coeff = basis_evaluate_t(basis, vols_out);
    
    if( unflatten_flag )
        vols_coeff = reshape(vols_coeff,[basis.count * r 1]);
    end
        
end

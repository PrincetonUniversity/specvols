%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the alternating minimization method on experimental data (here
% `frank70s_10k`). We do this by initializing with random weights and
% iterating a certain number of times. The performance is evaluated by
% clustering the coordinates and comparing the clusters with the
% “ground-truth” partition provided with the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 64;

num_eigs = 2;

C = 2;

star_opt = struct();
star_opt.pixel_size = 2.82;
star_opt.B = 0;

noise_est_opt = struct();
noise_est_opt.bg_radius = 60/65;
noise_est_opt.noise_type = 'anisotropic';

% Control the datatype and what tolerance we should have CG converge to.

precision = 'single';
rel_tolerance = 1e-6;

% Set the number of alternating minimization iterations.

altmin_iter = 10;

% Initialize basis and options structs.

basis = ffb_basis(L*ones(1, 3));

vols_wt_est_opt = struct();
vols_wt_est_opt.rel_tolerance = rel_tolerance;
vols_wt_est_opt.max_iter = 250;
vols_wt_est_opt.precision = precision;
vols_wt_est_opt.batch_size = 2^13;
vols_wt_est_opt.regularizer = 1e-6;

% Locate and load the data.

data_root = fullfile(pkg_root(), 'data', 'frank70s_10k');

filename = fullfile(data_root, 'run_001_it025_data.star');

star = load_star(filename);

src = star_to_src(star, data_root, star_opt);

% For this particular data, we know the “ground truth” class assignment, and
% so can use this to compare results.

true_idx = [1*ones(1, 5000) 2*ones(1, 5000)];
true_idx = true_idx(src.image_idx);

% Process the data (downsample, whiten).

src = downsample_src(src, L);

noise_psd_est = estimate_noise(src, noise_est_opt);

src = whiten_src(src, noise_psd_est);

% The data are flipped in sign, so incorporate this in the forward model.

src.params.filters = mult_filter(src.params.filters, scalar_filter(-1));

% Cache all these operations so we don't have to go back to disk.

src = cache_src(src);

% Initialize with random weights and normalize them. Then get the initial
% volume estimate.

wts = zeros(1+num_eigs, src.n);
wts(1,:) = ones(1, src.n);
wts(2:end,:) = randn(size(wts, 1)-1, size(wts, 2));

nrm = sqrt(sum(abs(wts).^2, 2));
wts = bsxfun(@times, wts, 1./nrm);

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% Run the alternating minimization from these random weights.

for k = 1:altmin_iter
    [Q, ~] = qr(vol_to_vec(vols_wt_est(:,:,:,2:end)), 0);
    vols_wt_est(:,:,:,2:end) = vec_to_vol(Q);

    coords = src_wiener_coords(src, ...
        vols_wt_est(:,:,:,1), vols_wt_est(:,:,:,2:end));

    wts(2:end,:) = coords;

    nrm = sqrt(sum(abs(wts).^2, 2));
    wts = bsxfun(@times, wts, 1./nrm);

    vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

    fprintf('Finished %d.\n', k);
end

% Cluster the resulting coordinates so that we can compare with the “ground
% truth.”

max_l_likelihood = -Inf;
for k = 1:10
    [~, idx_k, ~, ~, l_likelihood] = cluster_gmm(coords, C);

    if l_likelihood > max_l_likelihood
        max_l_likelihood = l_likelihood;
        idx = idx_k;
    end
end

accuracy = feval(@(p)(max(p, 1-p)), sum(idx == true_idx)/numel(idx));

print_value('Accuracy', corr(end,3));

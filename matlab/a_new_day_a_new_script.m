% A new day, a new script!
% By Amit Halevi
% New script to familiarize with Joakim's new, theoretically self-contained
% covariance estimation code.  Hooray!

%addpath!

% tic
addpath_cov3d;

%N is resolution (per side), n is number of projected images
% N = 8;
% n = 2^9-1;

if ~exist('vols') || N~= size(vols,1)
%     uncomputed_basis = fb_basis(N*ones(1,3),[]);
    uncomputed_basis = dirac_basis(N*ones(1,3),[]);
    basis = precompute_basis(uncomputed_basis);


    % number of eigenvolumes to extract from covariance matrix
    % num_eigs = 32;
    % number of diff ap coordinates
    % r = 5;
    % Diffusion time parameter for diffusion maps - large t makes things get
    % more "uniform" - will place increasing weight on the first map/subspace
    % of maps
    dmap_t = 0;

    %% let's generate volumes!

    num_angles_1 = 16;
    max_angle_1 = pi/2;
    num_angles_2 = 1;
    max_angle_2 = 0;

    % %  "Real way" to generate volumes: we have a rod and we attach spinning
    % %  things to it
    fixed_map = 'FakeKvMapAlphaOne.mrc';
    moving_map1 = 'FakeKvMapAlphaTwo.mrc';
    moving_map2 = 'FakeKvMapBeta.mrc';

    vol_fixed = load_mrc(fixed_map);
    vol_moving1 = load_mrc(moving_map1);
    vol_moving2 = load_mrc(moving_map2);

    vol_fixed = max(0, vol_fixed);
    vol_moving1 = max(0, vol_moving1);
    vol_moving2 = 0.85*max(0, vol_moving2);

    vol_moving1_rot = vol_rotate_z(vol_moving1, linspace(0, max_angle_1, num_angles_1));
    vol_moving2_rot = vol_rotate_z(vol_moving2, linspace(0, max_angle_2, num_angles_2));

    vol_moving_rot = bsxfun(@plus, vol_moving1_rot, permute(vol_moving2_rot, [1 2 3 5 4]));
    vol_moving_rot = reshape(vol_moving_rot, [size(vol_fixed) size(vol_moving1_rot, 4)*size(vol_moving2_rot, 4)]);

    vol_fixed = vol_downsample(vol_fixed, N);
    vol_moving_rot = vol_downsample(vol_moving_rot, N);

    vols = bsxfun(@plus, vol_fixed, vol_moving_rot);

    vols = basis_project(basis, vols);

    %we'll normalize anyway!
    avg_vol_energy = sum(vols(:).^2);
    if(avg_vol_energy ~= 0)
        normed_vols = vols/sqrt(avg_vol_energy) * sqrt(size(vols,4));
    else
        error('Vols have zero energy!');
    end
end
% normed_vols = vols / max(vols(:));

% Current implementation uses vols as an NxNxNx[num_vols] 4D array, not a
% cell of 3D arrays
% vols = num2cell(vols, 1:3);
% vols = vols(:);

%% First, we generate our sim!

% Output
%    sim: A simulation structure containing the fields:
%       - n: the number of images in the problem (default 1024),
%       - L: the resolution, that is L-by-L images and L-by-L-by-L volumes
%          (default 8),
%       - vols: An L-by-L-by-L-by-C array of volumes representing the C
%          volumes in the simulation (default `gaussian_blob_vols`),
%       - states: a 1-by-n array containing the hvolume states for each
%          images (default randomly sampled uniformly between 1 and K),
%       - rots: a 3-by-3-by-n array of rotation matrices corresponding to
%          viewing directions (default generated using 'rand_rots'),
%       - filters: a struct array of filter F objects (see `eval_filter`)
%          (default `identity_filter()`),
%       - filter_idx: a 1-by-n array containing the filter function assigments
%          for each of the images (default randomly sampled uniformly between
%          1 and F),
%       - offsets: a 2-by-n array specifying the shifts of the images (default
%          generated from a Gaussian distribution of standard deviation
%          L/16),
%       - amplitudes: a 1-by-n array specifying the amplitude multipliers of
%          the images (default uniformly sampled between 2/3 and 3/2),
%       - noise_seed: the random seed for generating the noise in the images
%           (default 0), and
%       - noise_var: the variance of the noise in the images (default 1).

sim_params = struct();
sim_params = fill_struct(sim_params, ...
        'n', n, ...
        'L', N,  ...
        'vols', normed_vols, ...
        'states', [], ...
        'rots', [], ...
        'filters', identity_filter(), ...
        'filter_idx', [], ...
        'offsets',zeros(2,n), ...
        'amplitudes', [], ...
        'noise_seed', 0, ...
        'noise_var', 0); 
    
% sim_params.filters = filter;
    
%fix states for practice!
sim_params.states = ceil(cumsum(ones(1,n))/n * num_angles_1 * num_angles_2);
    
% sim_params.states = ones(1,n);

sim = create_sim(sim_params);

% Normalize noise_var so that the number we put in is normalized by the
% amount of noise in each pixel

clean_images = sim_clean_image(sim,1,n);
sim_params.noise_var = sim_params.noise_var * ...
    (( sum(clean_images(:).^2)) / (sim_params.n * sim_params.L^2)) ;
sim.noise_var = sim_params.noise_var;


%Create the simulation!
src = sim_to_src(sim);

%% Now it's time to estimate covariance!

%Let's cheat for covariance!

mean_true = sim_mean(sim);
covar_true = sim_covar(sim);
[eigs_true, lambdas_true] = sim_eigs(sim);
[coords_true_states, residuals] = sim_vol_coords(sim, mean_true, eigs_true);
coords_true = coords_true_states(:,sim.states);
coords = coords_true;


%%  Diffusion maps!

n2 = sum(abs(coords).^2, 1);

d2 = bsxfun(@plus, n2, n2') - 2*coords'*coords;

dists = sqrt(max(0, d2));

%  As written, dists_to_dmap_coords assumes a Gaussian kernel on the
%  distances.  "Epsilon" is the 2 \sigma^2 term
%  Let's estimate a good value!  We use a subset in the event that n gets
%  very large.
% dists_sorted = sort(dists, 1);
% dists_epsilon = mean(mean(dists_sorted([1:128]+1,:)));
dists_epsilon = 0.2;   %empirically determined when using perfect

% Reminder: r and t are set at the top of the script
[dmap_coords dmap_evals]= dists_to_dmap_coords(dists, dists_epsilon, r, dmap_t);
% 
% disp('Done with calculating diffusion maps!')
% disp([num2str(toc,6) ' seconds elapsed.']);
% 
% %%  Diffusion map fitting!
% 
vols_wt_est_opt = struct();
vols_wt_est_opt = fill_struct(vols_wt_est_opt ,...
        'precision','single',...
        'batch_size',512, ...
        'preconditioner','none');
vols_wt_est_opt.max_iter = 500;
vols_wt_est_opt.verbose = 0;
vols_wt_est_opt.iter_callback = [];
vols_wt_est_opt.preconditioner = [];
vols_wt_est_opt.rel_tolerance = 1e-5;
vols_wt_est_opt.store_iterates = true;
% 
% % dmap_coords = coords_true;
% % r = size(dmap_coords,1);
% 
% wts = dmap_coords;
% 
% % [vols_wt_est cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
% % wts = ones(1,n);
% % "UNWRAPPED FUNCTION" estimate_vols_wt
% 
% %     wts = dmap_coords;
% % 
% % 
% % 
%      L = src.L;
%      n = src.n;
%  
%  %     vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
%  %         'preconditioner', 'circulant', ...
%  %         'precision', 'single');
%      vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
%          'preconditioner', 'none', ...
%          'precision', 'single');
%  %     if isempty(basis)
%  %         basis = dirac_basis(L*ones(1, 3));
%  %     end
%  
%      kermat_f = src_vols_wt_kermat(src, wts, vols_wt_est_opt);
%  
%      precond_kermat_f = [];
%  
%      if ischar(vols_wt_est_opt.preconditioner)
%          if strcmp(vols_wt_est_opt.preconditioner, 'none')
%              precond_kermat_f = [];
%          elseif strcmp(vols_wt_est_opt.preconditioner, 'circulant')
%              precond_kermat_f = circularize_kernel(kermat_f, 3);
%          else
%              error('Invalid preconditioner type.');
%          end
%  
%          % Reset so this is not used by the `conj_grad` function.
%          vols_wt_est_opt.preconditioner = @(x)(x);
%      end
%  
%      vols_wt_b_coeff = src_vols_wt_backward(src, basis, wts, vols_wt_est_opt);
%  
%      
%      [vols_wt_est_coeff, cg_info] = conj_grad_vols_wt(kermat_f, vols_wt_b_coeff, ...
%          basis, precond_kermat_f, vols_wt_est_opt);
%  
%      vols_wt_est = basis_evaluate(basis, vols_wt_est_coeff);
%      
% % END "UNWRAPPED FUNCTION"
% 
% toc

% dists = distances(coords');
% 
% dists_sorted = sort(dists, 1);
% % dists_epsilon = mean(mean(dists_sorted([1:128]+1,:)));

% covariance
% [dmap_coords, dmap_dists] = calc_dmap_coords(dists, dists_epsilon^2, 1, 32);
% dmap_dists_sorted = sort(dmap_dists, 1);
% dmap_dists_epsilon = mean(mean(dmap_dists_sorted([1:4]+1,:)));
% coords_ave = average_coords(coords, dmap_dists, ...
%     dmap_dists_epsilon, 'gaussian');
% 
% mu_true = simulation_mean(sim);
% [eigs_true, lambdas_true] = simulation_covariance(sim);
% eigs_true = eigs_true(:,:,:,1:num_eigs);
% lambdas_true = diag(lambdas_true);
% lambdas_true = lambdas_true(1:num_eigs);
% 
% mu_true = vol_downsample(mu_true, N_ds);
% eigs_true = vol_downsample(eigs_true, N_ds);
% 
% src_wi_true = wiener_filter_source(src, mu_true, eigs_true, lambdas_true, ...
%     sim.sigma);
% 
% coords_true = calc_coords(src_wi_true, mu_true, eigs_true);
% 
% rand('state', 0);
% rand(1, 1);
% [Q, ~] = qr(rand(num_eigs));
% 
% coords_proj = coords*Q;
% dmap_coords_proj = dmap_coords(:,1:num_eigs)*Q;
% coords_ave_proj = coords_ave*Q;
% coords_true_proj = coords_true*Q;
% 
% numIter = 25;
% 
% old_vol_norms = zeros(1,numIter);
% new_vol_norms = zeros(1,numIter);
% diff_vol_norms = zeros(1,numIter);
% 
% all_vols = zeros([N N N r numIter]);
% tic
% for iterate_index = 1:numIter
%     iterate_index
%     vols = zeros([N N N n]);
%     for k = 1:n
%         for l = 1:r
%         vols(:, :,:,k) = vols(:,:,:,k) + dmap_coords(l,n) * vols_wt_est(:,:,:,l);
%         end
%     end
%     normed_vols = vols/sqrt(sum(vols(:).^2)) * sqrt(size(vols,4));
% %     normed_vols = vols;
% %     sim_params = struct();
% %     sim_params = fill_struct(sim_params, ...
% %     'n', n, ...
% %     'L', N,  ...
% %     'vols', normed_vols, ...
% %     'states', [], ...
% %     'rots', [], ...
% %     'filters', identity_filter(), ...
% %     'filter_idx', [], ...
% %     'offsets', [], ...
% %     'amplitudes', [], ...
% %     'noise_seed', 0, ...
% %     'noise_var', 0);
%     %fix states for practice!
% %     sim_params.states = ceil(cumsum(ones(1,n))/n * num_angles_1);
%     sim = create_sim(sim_params);
%     sim.vols = normed_vols;
%     %Create the simulation!
%     src = sim_to_src(sim);
% %     [new_vols_wt_est] = estimate_vols_wt(src, basis, dmap_coords, vols_wt_est_opt);
% 
%     toc
%     
%     all_vols(:,:,:,:,iterate_index) = new_vols_wt_est;
%     
%     old_vol_norms(iterate_index) = sum(sum(sum(sum((vols_wt_est).^2))));
%     new_vol_norms(iterate_index) =  sum(sum(sum(sum((new_vols_wt_est).^2))));
%     diff_vol_norms(iterate_index) = sum(sum(sum(sum((new_vols_wt_est - vols_wt_est).^2))));
%     
%     vols_wt_est = new_vols_wt_est;
% 
% end
% 
% plot(1:numIter,old_vol_norms,'r',1:numIter,new_vol_norms,'b',1:numIter,diff_vol_norms,'k');

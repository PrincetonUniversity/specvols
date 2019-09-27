%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  __   __       __                 __      
% |__) |_  |\ | /   |__| |\/|  /\  |__) |_/ 
% |__) |__ | \| \__ |  | |  | /--\ | \  | \ 
% 
%          __  __  __     __  ___ 
%         (_  /   |__) | |__)  |  
%         __) \__ | \  | |     |  
% 
% A script to generate data for Channel_Stretch simulation
%
% Outline of script:
%   1. Initialize variables and settings
%   2. Set up test volumes
%   3. Set up sim object
%   4. Compute mean volume, covariance, covariance coordinates
%   5. Computer diffusion map (laplacian eigenvectors)
%   6. Estimate diffusion volumes
%   7. Check results
%   8. Outputs
%
%   Amit Halevi (ahalevi@princeton.edu)
%   September 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables and settings
% if this script was a function, these would be the inputs
%
% First, we addpath!

addpath_specvols;

% install_fftw;
% install_finufft;
install_mex;


disp('Let''s get started!');
disp(datetime('now'));
tic

max_threads = omp_get_max_threads;
threads_to_use = min(max_threads,16);
omp_set_num_threads(threads_to_use);

noise_levels = [30];
r = 15;

vols_wt_est_all = cell(length(noise_levels),r);
cov_eigs_all = cell(length(noise_levels));
coords_all = cell(length(noise_levels));
dmap_coords_all = cell(length(noise_levels));
lambdas_all = cell(length(noise_levels));
covs_all = cell(length(noise_levels));

% Now, the parameters relating to the volume
% 
num_angles_1 = 1;
max_angle_1 = 0;
num_angles_2 = 1;
max_angle_2 = 0;

num_shifts_x = 11;
min_shifts_x = -5;
max_shifts_x = 5;
num_shifts_y = 11;
min_shifts_y = -5;
max_shifts_y = 5;

shift_mask = ones(num_shifts_x,num_shifts_y);                                                 

[x_locs y_locs] = find(shift_mask);
shifts_x = x_locs - ceil(num_shifts_x/2);
shifts_y = y_locs - ceil(num_shifts_y/2);

shifts = [shifts_x(:) shifts_y(:)]';
  
% parameters related to the simulation

L = 108;
down_L = 16;
n = 1e5;
noise_var = 0;    %Frankly, not entirely certain how this is normalized
noise_seed = 0;     %for reproducibility
offsets = zeros(2,n); %currently unused
rots = [];          %[] means uniformly drawn from SO(3)
amplitudes = [];    %for image filters
%states = ceil(cumsum(ones(1,n))/n * num_angles_1 * num_angles_2);
                   %which volumes to use for each image.  They're placed 
                   %together here for ease of debugging
states = ceil(rand(1,n) * num_angles_1 * num_shifts_x * num_shifts_y);
%filters = identity_filter();
ctf_params = gen_ctf_params();
filters = ctf_filter(ctf_params);
filter_idx = [];    %draw which filter to use uniformly at random


% Parameters related to covariance stuff
% Note! You can cheat on any one of the estimations indivually, allowing
% you to e.g. run a real coordinate estimation with cheating covariance
num_cov_coords = 16;    %number of coordinates to extract from the
                       %covariance PCA thing.  These are used to calculate
                       %the adjacency matrix for the diffusion map
num_coords_used_for_real = 6;
mean_cheat = 0;
cov_cheat = 0;
eigs_cheat = 0;
coords_cheat = 0;

% Parameters related to diffusion maps
dmap_t = 0;              %time parameter for diffusion maps
%dists_epsilon = 0.2;      %width parameter for kernel used on the distances
knn_graph_k = 10;

% Parameters related to estimation of diffusion volumes
%r = 4;  %also counts for above...


% Parameters related to result checking

% for right now, we check several

disp(['Finished initializing, t = ' num2str(toc)]);
                   
%% Set up test volumes
uncomputed_basis = dirac_basis(L*ones(1,3),[]);
if basis_precompute
   basis = precompute_basis(uncomputed_basis);
else
   basis = uncomputed_basis;
end
down_basis = fb_basis(down_L * ones(1,3));

angles = zeros(1,size(shifts,2));

vol_tops = gen_fakekv_top_rotations(L,angles);
vol_bottoms = gen_fakekv_mid_bottom_shifts(L,shifts);
vols = vol_tops + vol_bottoms;                                                                

disp(['Finished with volumes, t = ' num2str(toc)]);
%% Set up sim object



sim_params = struct();
sim_params = fill_struct(sim_params, ...
       'n', n, ...
       'L', L,  ...
       'vols', vols, ...
       'states', states, ...
       'rots', rots, ...
       'filters', filters, ...
       'filter_idx', filter_idx, ...
       'offsets',offsets, ...
       'amplitudes', amplitudes, ...
       'noise_psd', scalar_filter(noise_var / L^3), ...
       'noise_seed', noise_seed);
sim = create_sim(sim_params);

uncached_src = sim_to_src(sim);

src = cache_src(uncached_src);
clean_src = src;
image_energy = sum(sum(sum(src.images.^2,3),2),1);
image_noise = randn([L L n]) * sqrt(image_energy / (L * L *n));

noisy_src = src;

down_src = (downsample_src(uncached_src,down_L));

clean_images = src.images;
noise_images = image_noise;

time_keeper = zeros(1,10000);
time_keeper_desc = cell(1,10000);
time_idx = 1;

noise_idx_min = 1;
noise_idx_max = 1;
r_used_min = 1;
r_used_max = min(r,16);

for noise_idx = noise_idx_min:noise_idx_max

noisy_src.images = clean_images + sqrt(noise_levels(noise_idx)) * noise_images;
src = noisy_src;

disp(['Finished setting up sim, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['Sim, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

%% We do the covariance thing

mean_est_opt = struct();
mean_est_opt = fill_struct(mean_est_opt ,...
       'precision','single',...
       'batch_size',2^16, ...
       'preconditioner','none');
mean_est_opt.max_iter = 500;
mean_est_opt.verbose = 0;
mean_est_opt.iter_callback = [];
mean_est_opt.preconditioner = [];
mean_est_opt.rel_tolerance = 1e-3;
mean_est_opt.store_iterates = true;

mean_vol = estimate_mean(down_src, down_basis, mean_est_opt);

disp(['Finished with mean, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['Mean, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

noise_var_est = estimate_noise_power(down_src);

cov_est_opt = struct();
cov_est_opt = fill_struct(cov_est_opt ,...
       'precision','single',...
       'batch_size',2^16, ...
       'preconditioner','none');
cov_est_opt.max_iter = 50;
cov_est_opt.verbose = 0;
cov_est_opt.iter_callback = [];
cov_est_opt.preconditioner = [];
cov_est_opt.rel_tolerance = 1e-3;
cov_est_opt.store_iterates = true;

covar = estimate_covar(down_src, mean_vol, noise_var_est, down_basis, cov_est_opt);

disp(['Finished with covariance, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['Cov, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

cov_flat = reshape(covar,[down_L^3 down_L^3]);
cov_flat_sym = 1/2 * (cov_flat + cov_flat');
cov_sym_unflat = reshape(cov_flat_sym, size(covar));
[cov_eigs, lambdas] = mdim_eigs(cov_sym_unflat, num_cov_coords, 'la');

cov_eigs_all{noise_idx} = cov_eigs;
lambdas_all{noise_idx} = lambdas;
covs_all{noise_idx} = covar;

disp(['Finished with covar eigs, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['Cov eigs, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

noise_var_est = estimate_noise_power(down_src);
coords = src_wiener_coords(down_src, mean_vol, cov_eigs, lambdas, ...
    noise_var_est);

coords_trunc = coords(1:num_coords_used_for_real,:);
coords_all{noise_idx} = coords;

disp(['Finished with covar coords, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['Cov coords, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

disp(['Finished with all covar stuff, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['All cov stuff, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

%% Diffusion map calculation

graph_weights = graph_knn(coords', knn_graph_k);
graph_laplacian = laplacian(graph_weights, 'normalized');
[dmap_coords, dmap_evals] = eigs(graph_laplacian, r, 'smallestabs');
dmap_coords = dmap_coords';
dmap_coords_all{noise_idx} = dmap_coords;

disp(['Finished with dmap coords, t = ' num2str(toc)]);
time_keeper_desc{time_idx} = ['dmap_coords, noise_idx = ' num2str(noise_idx)];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

% Time to estimate vols!

vols_wt_est_opt = struct();
vols_wt_est_opt = fill_struct(vols_wt_est_opt ,...
        'precision','single',...
        'batch_size',3e4, ...
        'preconditioner','none');
vols_wt_est_opt.max_iter = 50;
vols_wt_est_opt.verbose = 0;
vols_wt_est_opt.iter_callback = [];
vols_wt_est_opt.preconditioner = [];
vols_wt_est_opt.rel_tolerance = 1e-5;
vols_wt_est_opt.store_iterates = true;

src = src;
basis = basis;
wts = dmap_coords;
vols_wt_est_opt = vols_wt_est_opt;

    L = src.L;
    n = src.n;
    vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
        'preconditioner', 'none', ...
        'precision', 'single');
    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kermat_f = sqrt(n^2) * src_vols_wt_kermat(src, wts, vols_wt_est_opt);

    disp(['Finished with kermat, t = ' num2str(toc)]);
    time_keeper_desc{time_idx} = ['Kermat, noise_idx = ' num2str(noise_idx)];
    time_keeper(time_idx) = toc;time_idx = time_idx + 1;

    precond_kermat_f = [];

    if ischar(vols_wt_est_opt.preconditioner)
        if strcmp(vols_wt_est_opt.preconditioner, 'none')
            precond_kermat_f = [];
        elseif strcmp(vols_wt_est_opt.preconditioner, 'circulant')
            precond_kermat_f = circularize_kernel(kermat_f, 3);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        vols_wt_est_opt.preconditioner = @(x)(x);
    end

    vols_wt_b_coeff = sqrt(n) * src_vols_wt_backward(src, basis, wts, vols_wt_est_opt);
    disp(['Finished with backproj, t = ' num2str(toc)]);
    time_keeper_desc{time_idx} = ['backproj, noise_idx = ' num2str(noise_idx)];
    time_keeper(time_idx) = toc;time_idx = time_idx + 1;

        for r_used = r_used_min:r_used_max
            [vols_wt_est_coeff,cg_info]=conj_grad_vols_wt(kermat_f(:,:,:,1:r_used,1:r_used)...
                , vols_wt_b_coeff(:,1:r_used),basis, precond_kermat_f, vols_wt_est_opt);
            vols_wt_est = basis_evaluate(basis, vols_wt_est_coeff);
            vols_wt_est_all{noise_idx,r_used} = vols_wt_est;
            disp(['Done with estimation! Noise level = ' num2str(noise_levels(noise_idx)) ...
                ', r_used = ' num2str(r_used) ', t = ' num2str(toc)]);
            time_keeper_desc{time_idx} = ['Est, noise_idx = ' num2str(noise_idx) ', r_used = ' num2str(r_used)];
            time_keeper(time_idx) = toc;time_idx = time_idx + 1;
        end

end %this "end" closes the for loop over noise levels
%%End vols_wt_est_unwrapped

%%
disp(['Finished estimating vols, t = ' num2str(toc)]);
%

time_keeper_desc{time_idx} = ['Finished estimating vols!'];
time_keeper(time_idx) = toc;time_idx = time_idx + 1;

disp('All done!');
disp(datetime('now'));


% errs_mse = zeros(length(rs),n);
% errs_corr = zeros(length(rs),n);
errs_mse = cell(length(script_params),1);
errs_corr = cell(length(script_params),1);

diff_vols = cell(length(script_params),5);
diff_cg_info = cell(length(script_params),5);
diff_srcs = cell(length(script_params),5);
diff_wts = cell(length(script_params),5);

disp(datetime('now'))

%add pool startup

tic;
for k = 1:length(script_params)
    N = script_params(k).N;
    n = script_params(k).n;
    r = script_params(k).r;
    a_new_day_a_new_script;

    
    wts = dmap_coords;
    [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
    diff_wts{k,1} = wts;
    diff_vols{k,1} = vols_wt_est;
    diff_cg_info{k,1} = cg_info;
    diff_srcs{k,1} = src;
    
    disp(['k = ' num2str(k) ', first time, t = ' num2str(toc)]);
    
%     recon_vols = reshape( reshape(vols_wt_est,[N^3 r]) * wts,[N N N n]);
%     sim.states = 1:n;
%     sim.vols = recon_vols;
%     src = sim_to_src(sim);
%     [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
%     diff_wts{k,2} = wts;
%     diff_vols{k,2} = vols_wt_est;
%     diff_cg_info{k,2} = cg_info;
%     diff_srcs{k,2} = src;
%     disp(['k = ' num2str(k) ', second time, t = ' num2str(toc)]);
    
%     recon_vols = reshape( reshape(vols_wt_est,[N^3 r]) * wts,[N N N n]);
%     sim.states = 1:n;
%     sim.vols = recon_vols;
%     src = sim_to_src(sim);
%     [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
%     diff_wts{k,3} = wts;
%     diff_vols{k,3} = vols_wt_est;
%     diff_cg_info{k,3} = cg_info;
%     diff_srcs{k,3} = src;
%     disp(['k = ' num2str(k) ', third time, t = ' num2str(toc)]);
% 
%     recon_vols = reshape( reshape(vols_wt_est,[N^3 r]) * wts,[N N N n]);
%     sim.states = 1:n;
%     sim.vols = recon_vols;
%     src = sim_to_src(sim);
%     [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
%     diff_wts{k,4} = wts;
%     diff_vols{k,4} = vols_wt_est;
%     diff_cg_info{k,4} = cg_info;
%     diff_srcs{k,4} = src;
%     disp(['k = ' num2str(k) ', fourth time, t = ' num2str(toc)]);
% 
%     recon_vols = reshape( reshape(vols_wt_est,[N^3 r]) * wts,[N N N n]);
%     sim.states = 1:n;
%     sim.vols = recon_vols;
%     src = sim_to_src(sim);
%     [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);
%     diff_wts{k,5} = wts;
%     diff_vols{k,5} = vols_wt_est;
%     diff_cg_info{k,5} = cg_info;
%     diff_srcs{k,5} = src;
%     disp(['k = ' num2str(k) ', fifth time, t = ' num2str(toc)]);

    %     r = rs(k);
%     [dmap_coords dmap_evals]= dists_to_dmap_coords(dists, dists_epsilon, r, dmap_t);
%     wts = dmap_coords;
%     
%     kermat_f = src_vols_wt_kermat(src, wts, vols_wt_est_opt);
% 
%     vols_wt_est_opt = struct();
%     vols_wt_est_opt = fill_struct(vols_wt_est_opt ,...
%             'precision','double',...
%             'batch_size',512, ...
%             'preconditioner','none');
%     vols_wt_est_opt.max_iter = 50;
%     vols_wt_est_opt.verbose = 0;
%     vols_wt_est_opt.iter_callback = [];
%     vols_wt_est_opt.preconditioner = [];
%     vols_wt_est_opt.rel_tolerance = 1e-5;
%     vols_wt_est_opt.store_iterates = true;
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
%     errs_corr{k} = check_recon_wts(src,wts,vols_wt_est,'corr');
%     errs_mse{k} = check_recon_wts(src,wts,vols_wt_est,'mse');
%     
%     cg_infos{k} = cg_info;
% %     mean_errs{k} = check_recon_wts(src,wts,vols_wt_est,'square');
%     disp(['r = ' num2str(r) ' time = ' num2str(toc)]);
end
% 
% subplot_width = ceil(sqrt(2*length(script_params)));
% subplot_height = ceil(length(script_params) / subplot_width);
% 
% figure
% for k = 1:length(script_params)
%     subplot_loc = k + ((fix((k-1)/subplot_width)) * subplot_width);
%     subplot(2*subplot_height,subplot_width, subplot_loc )
%     plot(errs_corr{k})
%     ylabel('Corr');
%     title([num2str(script_params(k).n) ' n, mean corr = ' num2str(mean(errs_corr{k}))]);
%     subplot(2*subplot_height,subplot_width,subplot_loc+subplot_width)
%     plot(errs_mse{k})
%     ylabel('MSE');
%     title([num2str(script_params(k).n) ' n, MSE = ' num2str(mean(errs_mse{k}))]);
% end
% suptitle(['N = ' num2str(N) ', n = ' num2str(n) ', SNR = ' num2str(1 / ...
%     sim_params.noise_var * ...
%     (( sum(clean_images(:).^2)) / (sim_params.n * sim_params.L^2)) )]);

for k = 1:length(script_params)
    
    N = script_params(k).N;
    n = script_params(k).n;
    r = script_params(k).r;

    figure(k)
    subplot(4,5,1)
    plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,1},'mse'));
    title('mse vs orig')
    subplot(4,5,6)
    plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,1},'corr'));
    title('corr vs orig')

%     subplot(4,5,2)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,2},'mse'));
%     title('mse vs orig')
%     subplot(4,5,7)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,2},'corr'));
%     title('corr vs orig')
%     subplot(4,5,12)
%     plot(check_recon_wts(diff_srcs{k,2},diff_wts{k,1},diff_vols{k,2},'mse'));
%     title('mse vs this round')
%     subplot(4,5,17)
%     plot(check_recon_wts(diff_srcs{k,2},diff_wts{k,1},diff_vols{k,2},'corr'));
%     title('corr vs this round')

%     subplot(4,5,3)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,3},'mse'));
%     title('mse vs orig')
%     subplot(4,5,8)
%     title('corr vs orig')
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,3},'corr'));
%     subplot(4,5,13)
%     plot(check_recon_wts(diff_srcs{k,3},diff_wts{k,1},diff_vols{k,3},'mse'));
%     title('mse vs this round')
%     subplot(4,5,18)
%     plot(check_recon_wts(diff_srcs{k,3},diff_wts{k,1},diff_vols{k,3},'corr'));
%     title('corr vs this round')
% 
%     subplot(4,5,4)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,4},'mse'));
%     title('mse vs orig')
%     subplot(4,5,9)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,4},'corr'));
%     title('corr vs orig')
%     subplot(4,5,14)
%     plot(check_recon_wts(diff_srcs{k,4},diff_wts{k,1},diff_vols{k,4},'mse'));
%     title('mse vs this round')
%     subplot(4,5,19)
%     plot(check_recon_wts(diff_srcs{k,4},diff_wts{k,1},diff_vols{k,4},'corr'));
%     title('corr vs this round')
% 
%     subplot(4,5,5)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,5},'mse'));
%     title('mse vs orig')
%     subplot(4,5,10)
%     plot(check_recon_wts(diff_srcs{k,1},diff_wts{k,1},diff_vols{k,5},'corr'));
%     title('corr vs orig')
%     subplot(4,5,15)
%     plot(check_recon_wts(diff_srcs{k,5},diff_wts{k,1},diff_vols{k,5},'mse'));
%     title('mse vs this round')
%     subplot(4,5,20)
%     plot(check_recon_wts(diff_srcs{k,5},diff_wts{k,1},diff_vols{k,5},'corr'));
%     title('corr vs this round')

    suptitle(['N = ' num2str(N) ', n = ' num2str(n) ', r = ' num2str(r)]);
    
    %add pool shutdown
    p = gcp('nocreate');
    delete(p);
end
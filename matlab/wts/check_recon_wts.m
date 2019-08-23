% CHECK_RECON_WTS Check errors for wted vol reconstruction
%
% Usage
%    errs = check_recon_wts(src, wts, vols_wt_est, err_type)
%
% Input
%    src: a src object from `sim_to_src`
%    wts: an n x r matrix specifying the weights for each volume at each
%    respective image (i.e. coordinates)
%    vols_wts_est: an N x N x N x r array of weighted basis/"diffusion"
%    volumes
%    err_type: the type of error to calculate (default squared error)
%    batch_size: the size of the batch to use; default n
%
% Output
%    errs: a vector of length n giving the error between each 
%
% See also
%    estimate_vols_wt

% Author
%    Amit Halevi <ahalevi@princeton.edu>


function [ errs ] = check_recon_wts( src, wts, vols_wt_est, err_type, batch_size )

    N = src.L;
    n = src.n;
    r = size(wts,1);
    
    if nargin < 4
        err_type = 'mse';
    end
    %Add batch size?
    if nargin < 5
        batch_size = n;
    end

    if strcmp(err_type,'fsc_hollow')
        mean_vol = sim_mean(src.sim);
        hollow_orig_vols = zeros([N N N size(src.sim.vols,4)]);
        for vol_idx = 1:size(hollow_orig_vols,4)
            hollow_orig_vols(:,:,:,vol_idx) = src.sim.vols(:,:,:,vol_idx) - mean_vol;
        end
    end
    
    
    if(strcmp(err_type,'fsc') || strcmp(err_type,'fsc_hollow'))
        errs = zeros(floor(N/2),n);    %Floor for odd N
    else
        errs = zeros(1,n);
    end

    
%    unique_states = unique(src.sim.states);
    
    num_batches = ceil(n / batch_size);

    for batch_idx = 1:num_batches
    batch_start = ((batch_idx - 1) * batch_size + 1);
    batch_end = min(batch_idx * batch_size,n);

%    for i = 1:length(unique_states)
%        unbatch_idxes = find(src.sim.states == unique_states(i));
        unbatch_idxes = 1:n;
        idxes = unbatch_idxes((batch_start <= unbatch_idxes) & (unbatch_idxes <= batch_end));
    
    %if(~isempty(idxes))
    
    if strcmp(err_type,'fsc_hollow')
        hollow_recon_vols_flat = reshape(vols_wt_est(:,:,:,2:end), [N^3 r-1]) * wts(2:end,idxes);
        hollow_recon_vols = reshape(hollow_recon_vols_flat,[N N N length(idxes)]);
    else
        recon_vols_flat = reshape(vols_wt_est,[N^3 r]) * wts(:,idxes);
        recon_vols = reshape(recon_vols_flat,[N N N length(idxes)]);
    end
    batch_errs = zeros([size(errs,1) length(idxes)]);
    for k = 1:length(idxes)    
%        if(strcmp(err_type,'mse'))
%            errs(idxes) = sum(sum(sum( (bsxfun(@minus,src.sim.vols(...
%                :,:,:,unique_states(i)),recon_vols).^2),1),2),3);
%        elseif(strcmp(err_type,'corr'))
%            errs(idxes) = acorr(repmat(src.sim.vols(...
%                :,:,:,unique_states(i)),[1 1 1 length(idxes)]),recon_vols,[1 2 3]);
%        elseif(strcmp(err_type,'ratios'))
%            errs(idxes) = [];
%        elseif(strcmp(err_type,'fsc'))
%            errs(:,idxes) = repmat(FSCorr(src.sim.vols(:,:,:, ...
%                unique_states(i)),recon_vols),[1, numel(idxes)]);
%        elseif(strcmp(err_type,'fsc_hollow'))
%            errs(:,idxes) = repmat( FSCorr( ...
%                src.sim.vols(:,:,:,unique_states(i)) - mean_vol, ...
%                hollow_recon_vols), [1, numel(idxes)]);
%        end
         if(strcmp(err_type,'fsc'))
             batch_errs(:,k) = FSCorr(recon_vols(:,:,:,k),src.sim.vols(:,:,:,src.sim.states(idxes(k))));
         elseif(strcmp(err_type,'fsc_hollow'))
             batch_errs(:,k) = FSCorr(hollow_recon_vols(:,:,:,k),hollow_orig_vols(:,:,:,src.sim.states(idxes(k))));
%         end
    end %for k = 1:length(idxes)
    errs(:,batch_start:batch_end) = batch_errs;
  %  end %if(~isempty(idxes))
%    end % for i = 1:length(unique_states)

    end % for batch_idx = 1:num_batches

end


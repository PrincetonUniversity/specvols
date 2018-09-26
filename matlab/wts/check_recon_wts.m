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
%
% Output
%    errs: a vector of length n giving the error between each 
%
% See also
%    estimate_vols_wt

% Author
%    Amit Halevi <ahalevi@princeton.edu>


function [ errs ] = check_recon_wts( src, wts, vols_wt_est, err_type )

    N = src.L;
    n = src.n;
    r = size(wts,1);
    
    if nargin < 4
        err_type = 'mse';
    end
    %Add batch size?
    
    recon_vols = reshape( reshape(vols_wt_est,[N^3 r]) * wts,[N N N n]);
    
    if(strcmp(err_type,'FSC'))
        errs = zeros(floor(N/2),n);    %Floor for odd N
    else
        errs = zeros(1,n);
    end

    
    unique_states = unique(src.sim.states);
    
    for i = 1:length(unique_states)
        idxes = find(src.sim.states == unique_states(i));
        
        if(strcmp(err_type,'mse'))
            errs(idxes) = sum(sum(sum( (bsxfun(@minus,src.sim.vols(...
                :,:,:,unique_states(i)),recon_vols(:,:,:,idxes)).^2),1),2),3);
        elseif(strcmp(err_type,'corr'))
            errs(idxes) = acorr(repmat(src.sim.vols(...
                :,:,:,unique_states(i)),[1 1 1 length(idxes)]),recon_vols(:,:,:,idxes),[1 2 3]);
        elseif(strcmp(err_type,'ratios'))
            errs(idxes) = [];
        elseif(strcmp(err_type,'FSC'))
            errs(:,idxes) = repmat(FSCorr(src.sim.vols(:,:,:, ...
                unique_states(i)),recon_vols(:,:,:,idxes)),[1, numel(idxes)]);
        end
    end

end


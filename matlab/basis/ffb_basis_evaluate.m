% FFB_BASIS_EVALUATE Evaluate coefficient vector in fast Fourier-Bessel basis
%
% Usage
%    x = ffb_basis_evaluate(basis, v);
%
% Input/Output
%    See documentation for `basis_evaluate`.
%
% See also
%    basis_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = ffb_basis_evaluate(basis, v)
    d = numel(basis.sz);

    if d == 2
        x = ffb_basis_evaluate_2d(basis, v);
    elseif d == 3
        x = ffb_basis_evaluate_3d(basis, v);
    end
end

function x = ffb_basis_evaluate_2d(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    n_theta = size(basis.precomp.freqs, 3);
    n_r = size(basis.precomp.freqs, 2);
    n_data = size(v, 2);

    % TODO: Rename. This is not actually the polar FT.
    pf = zeros([n_r 2*n_theta n_data], class(v));

    mask = (basis.indices.ells == 0);

    ind = 1;

    idx = ind + [0:basis.k_max(1)-1];

    pf(:,1,:) = basis.precomp.radial(:,idx)*v(mask,:);

    ind = ind + numel(idx);

    ind_pos = ind;

    for ell = 1:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        idx_pos = ind_pos + [0:basis.k_max(ell+1)-1];
        idx_neg = idx_pos + basis.k_max(ell+1);

        v_ell = (v(idx_pos,:) - 1i*v(idx_neg,:))/2;

        if mod(ell, 2) == 1
            v_ell = 1i*v_ell;
        end

        pf_ell = basis.precomp.radial(:,idx)*v_ell;
        pf(:,ell+1,:) = pf_ell;

        if mod(ell, 2) == 0
            pf(:,end-ell+1,:) = conj(pf_ell);
        else
            pf(:,end-ell+1,:) = -conj(pf_ell);
        end

        ind = ind + numel(idx);

        ind_pos = ind_pos + 2*basis.k_max(ell+1);
    end

    pf = 2*pi*ifft(pf, [], 2);

    % Only need "positive" frequencies.
    pf = pf(:,1:end/2,:);

    pf = bsxfun(@times, pf, basis.precomp.w.*basis.precomp.r);

    pf = reshape(pf, [n_r*n_theta n_data]);
    freqs = reshape(basis.precomp.freqs, [2 n_r*n_theta]);

    x = zeros([basis.sz n_data], class(pf));
    for ell = 1:n_data
        x(:,:,ell) = 2*real(anufft2(pf(:,ell), 2*pi*freqs, basis.sz));
    end

    x = roll_dim(x, sz_roll);
end

function x = ffb_basis_evaluate_3d(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    n_data = size(v, 2);

    n_r = size(basis.precomp.radial_wtd, 1);
    n_phi = size(basis.precomp.angular_phi_wtd_even{1}, 1);
    n_theta = size(basis.precomp.angular_theta_wtd, 1);

    u_even = zeros( ...
        [n_r 2*basis.ell_max+1 n_data floor(basis.ell_max/2)+1], class(v));
    u_odd = zeros( ...
        [n_r 2*basis.ell_max+1 n_data ceil(basis.ell_max/2)], class(v));

    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_wtd = basis.precomp.radial_wtd(:,1:k_max_ell,ell+1);

        % TODO: Fix this to avoid lookup each time.
        ind = (basis.indices.ells == ell);

        v_ell = reshape(v(ind,:), [k_max_ell (2*ell+1)*n_data]);

        v_ell = radial_wtd*v_ell;

        v_ell = reshape(v_ell, [n_r 2*ell+1 n_data]);

        if mod(ell, 2) == 0
            u_even(:,[-ell:ell]+basis.ell_max+1,:,ell/2+1) = v_ell;
        else
            u_odd(:,[-ell:ell]+basis.ell_max+1,:,(ell-1)/2+1) = v_ell;
        end
    end

    u_even = permute(u_even, [4 1 2 3]);
    u_odd = permute(u_odd, [4 1 2 3]);

    w_even = zeros([n_phi n_r n_data 2*basis.ell_max+1], class(v));
    w_odd = zeros([n_phi n_r n_data 2*basis.ell_max+1], class(v));

    for m = 0:basis.ell_max
        angular_phi_wtd_m_even = basis.precomp.angular_phi_wtd_even{m+1};
        angular_phi_wtd_m_odd = basis.precomp.angular_phi_wtd_odd{m+1};

        n_even_ell = size(angular_phi_wtd_m_even, 2);
        n_odd_ell = size(angular_phi_wtd_m_odd, 2);

        if m == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            u_m_even = u_even(end-n_even_ell+1:end,:,basis.ell_max+1+sgn*m,:);
            u_m_odd = u_odd(end-n_odd_ell+1:end,:,basis.ell_max+1+sgn*m,:);

            u_m_even = reshape(u_m_even, [n_even_ell n_r*n_data]);
            u_m_odd = reshape(u_m_odd, [n_odd_ell n_r*n_data]);

            w_m_even = angular_phi_wtd_m_even*u_m_even;
            w_m_odd = angular_phi_wtd_m_odd*u_m_odd;

            w_m_even = reshape(w_m_even, [n_phi n_r n_data]);
            w_m_odd = reshape(w_m_odd, [n_phi n_r n_data]);

            w_even(:,:,:,basis.ell_max+1+sgn*m) = w_m_even;
            w_odd(:,:,:,basis.ell_max+1+sgn*m) = w_m_odd;
        end
    end

    w_even = permute(w_even, [4 1 2 3]);
    w_odd = permute(w_odd, [4 1 2 3]);

    u_even = w_even;
    u_odd = w_odd;

    u_even = reshape(u_even, [2*basis.ell_max+1 n_phi*n_r*n_data]);
    u_odd = reshape(u_odd, [2*basis.ell_max+1 n_phi*n_r*n_data]);

    w_even = basis.precomp.angular_theta_wtd*u_even;
    w_odd = basis.precomp.angular_theta_wtd*u_odd;

    pf = w_even + 1i*w_odd;
    pf = reshape(pf, [n_theta*n_phi*n_r n_data]);

    x = zeros([basis.sz n_data], class(pf));
    for ell = 1:n_data
        x(:,:,:,ell) = ...
            real(anufft3(pf(:,ell), basis.precomp.fourier_pts, basis.sz));
    end

    x = cast(x, class(v));

    x = roll_dim(x, sz_roll);
end

% FFB_BASIS_EVALUATE_T Evaluate coefficient in dual fast Fourier-Bessel basis
%
% Usage
%    v = ffb_basis_evaluate_t(basis, x);
%
% Input/Output
%    See documentation for `basis_evaluate_t`.
%
% See also
%    basis_evaluate_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = ffb_basis_evaluate_t(basis, x)
    d = numel(basis.sz);

    if d == 2
        v = ffb_basis_evaluate_t_2d(basis, x);
    elseif d == 3
        v = ffb_basis_evaluate_t_3d(basis, x);
    end
end

function v = ffb_basis_evaluate_t_2d(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    n_theta = size(basis.precomp.freqs, 3);
    n_r = size(basis.precomp.freqs, 2);
    n_data = size(x, 3);

    freqs = reshape(basis.precomp.freqs, [2 n_r*n_theta]);
    pf = zeros([size(freqs, 2) n_data], class(x));
    for ell = 1:n_data
        pf(:,ell) = nufft2(x(:,:,ell), 2*pi*freqs);
    end
    pf = reshape(pf, [n_r n_theta n_data]);

    % Recover "negative" frequencies from "positive" half plane.
    pf = cat(2, pf, conj(pf));

    pf = bsxfun(@times, pf, basis.precomp.w.*basis.precomp.r);

    % TODO: Rename. This isn't actually the polar FT.
    pf = 2*pi/(2*n_theta)*fft(pf, [], 2);

    % This only makes it easier to slice the array later.
    pf = permute(pf, [1 3 2]);

    v = zeros([basis.count n_data], class(x));

    ind = 1;

    idx = ind + [0:basis.k_max(1)-1];

    mask = (basis.indices.ells == 0);

    v(mask,:) = basis.precomp.radial(:,idx)'*real(pf(:,:,1));

    ind = ind + numel(idx);

    ind_pos = ind;

    for ell = 1:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        idx_pos = ind_pos + [0:basis.k_max(ell+1)-1];
        idx_neg = idx_pos + basis.k_max(ell+1);

        v_ell = basis.precomp.radial(:,idx)'*pf(:,:,ell+1);

        if mod(ell, 2) == 0
            v_pos = real(v_ell);
            v_neg = -imag(v_ell);
        else
            v_pos = imag(v_ell);
            v_neg = real(v_ell);
        end

        v(idx_pos,:) = v_pos;
        v(idx_neg,:) = v_neg;

        ind = ind + numel(idx);

        ind_pos = ind_pos + 2*basis.k_max(ell+1);
    end

    v = roll_dim(v, sz_roll);
end

function v = ffb_basis_evaluate_t_3d(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    n_data = size(x, 4);

    n_r = size(basis.precomp.radial_wtd, 1);
    n_phi = size(basis.precomp.angular_phi_wtd_even{1}, 1);
    n_theta = size(basis.precomp.angular_theta_wtd, 1);

    pf = zeros([size(basis.precomp.fourier_pts, 2) n_data], class(x));
    for ell = 1:n_data
        pf(:,ell) = nufft3(x(:,:,:,ell), basis.precomp.fourier_pts);
    end

    pf = cast(pf, class(x));

    pf = reshape(pf, [n_theta n_phi*n_r*n_data]);

    u_even = basis.precomp.angular_theta_wtd'*real(pf);
    u_odd = basis.precomp.angular_theta_wtd'*imag(pf);

    u_even = reshape(u_even, [2*basis.ell_max+1 n_phi n_r n_data]);
    u_odd = reshape(u_odd, [2*basis.ell_max+1 n_phi n_r n_data]);

    u_even = permute(u_even, [2 3 4 1]);
    u_odd = permute(u_odd, [2 3 4 1]);

    w_even = zeros( ...
        [floor(basis.ell_max/2)+1 n_r 2*basis.ell_max+1 n_data], class(x));
    w_odd = zeros( ...
        [ceil(basis.ell_max/2) n_r 2*basis.ell_max+1 n_data], class(x));

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
            u_m_even = u_even(:,:,:,basis.ell_max+1+sgn*m);
            u_m_odd = u_odd(:,:,:,basis.ell_max+1+sgn*m);

            u_m_even = reshape(u_m_even, [n_phi n_r*n_data]);
            u_m_odd = reshape(u_m_odd, [n_phi n_r*n_data]);

            w_m_even = angular_phi_wtd_m_even'*u_m_even;
            w_m_odd = angular_phi_wtd_m_odd'*u_m_odd;

            w_m_even = reshape(w_m_even, [n_even_ell n_r n_data]);
            w_m_odd = reshape(w_m_odd, [n_odd_ell n_r n_data]);

            w_even(end-n_even_ell+1:end,:,basis.ell_max+1+sgn*m,:) = w_m_even;
            w_odd(end-n_odd_ell+1:end,:,basis.ell_max+1+sgn*m,:) = w_m_odd;
        end
    end

    w_even = permute(w_even, [2 3 4 1]);
    w_odd = permute(w_odd, [2 3 4 1]);

    v = zeros([basis.count n_data], class(x));

    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_wtd = basis.precomp.radial_wtd(:,1:k_max_ell,ell+1);

        if mod(ell, 2) == 0
            v_ell = w_even(:,[-ell:ell]+basis.ell_max+1,:,ell/2+1);
        else
            v_ell = w_odd(:,[-ell:ell]+basis.ell_max+1,:,(ell-1)/2+1);
        end

        v_ell = reshape(v_ell, [n_r (2*ell+1)*n_data]);

        v_ell = radial_wtd'*v_ell;

        v_ell = reshape(v_ell, [k_max_ell*(2*ell+1) n_data]);

        % TODO: Fix this to avoid lookup each time.
        ind = (basis.indices.ells == ell);

        v(ind,:) = v_ell;
    end

    v = roll_dim(v, sz_roll);
end

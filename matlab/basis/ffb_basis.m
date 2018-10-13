% FFB_BASIS Construct a fast Fourier-Bessel basis object
%
% Usage
%    basis = ffb_basis(sz, ell_max, domain);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only square images are supported.
%    ell_max: The maximum order ell of the basis elements. If set to Inf,
%       the basis includes all ell such that the resulting basis vectors
%       are concentrated below the Nyquist frequency (default Inf).
%    domain: Specifies whether the decomposition should be in the spatial (0)
%       or frequency (1) domain. Currently, only frequency-domain basis
%       vectors are supported (default 1).
%
% Output
%    basis: A fast Fourier-Bessel basis object corresponding to the parameters.
%
% Description
%    The returned basis object can be used to project a vector x in the
%    standard basis onto the basis or to map a coefficient vector v into
%    the standard basis. This is achieved using the basis_expand and
%    basis_evaluate functions, respectively. The fields basis.sz denotes
%    the size of the vectors in the standard basis while basis.count denotes
%    the number of vectors in the basis.
%
%    For example, we can generate a random vector in the standard basis and
%    project it onto the Fourier-Bessel basis through
%
%       sz = 16*ones(1, 2);
%       basis = ffb_basis(sz, floor(sz(1)/2));
%       x = randn(sz);
%       v = basis_expand(basis, x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis_evaluate(basis, v);
%
%    The adjoint of the evaluate function is obtained through
%    basis_evaluate_t.
%
%    This particular basis is a Fourier-Bessel basis in two or three
%    dimensions. It is a separable basis consisting of a radial part multiplied
%    by an angular part.
%
%    The radial part is a Bessel function in 2D or a spherical Bessel function
%    in 3D. These are dilated to obtain a given number of zeros on the
%    interval [0, 1]. Specifically, they are of the form
%
%       f_{ell,k}(r) = Z_ell(R_{ell,k} r),
%
%    where Z_ell is the ellth-order (spherical) Bessel function and R_{ell,k}
%    is its kth zero.
%
%    The angular part is given by a sinusoid in 2D or a real spherical harmonic
%    function in 3D. Specifically, for ell = 0, the angular part is constant,
%    while for ell > 0, it is either a sine or a cosine.
%
%    All of the basis vectors are normalized so that their 2-norm, as continuous
%    functions, is equal to 1. Additionally, the basis functions are orthogonal,
%    again in the continuous case. Due to the discrete sampling, these
%    properties only hold asymptotically as max(sz) goes to infinity.
%
%    The parameter ell_max determines the highest frequency of the angular
%    sinusoids in 2D and the highest degree of the real spherical harmonics in
%    3D. For each ell, the maximum number of zeros in the radial part f_{ell,k}
%    is determined by ensuring that the Fourier transform of the entire basis
%    function is mostly concentrated within the Nyquist disk or ball. This
%    maximum number is denoted by k_max(ell) and is stored in the basis object.
%
%    The coefficients of the basis are ordered by increasing ell, from 0 to
%    ell_max. For each ell, the basis functions are ordered by the radial index
%    k from 1 to k_max(ell). In 2D, for each radial index, we then have the
%    cosine part followed by the sine part (except for ell = 0, in which case we
%    only have the cosine part). In 3D, we instead have the coefficients
%    arranged by order m, ranging from -ell to +ell.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function basis = ffb_basis(sz, ell_max, domain)
    if nargin < 2 || isempty(ell_max)
        ell_max = Inf;
    end

    if nargin < 3 || isempty(domain)
        domain = 1;
    end

    if numel(sz) ~= 3 && numel(sz) ~= 2
        error('Only two- or three-dimesional basis functions are supported.');
    end

    if ~all(sz == sz(1))
        error('Only cubic domains are supported.');
    end

    d = numel(sz);

    N = sz(1);

    % 2*sz(1) is an upper bound for ell_max, so use it here.
    k_max = zeros(min(ell_max+1, 2*sz(1)+1), 1);

    r0 = cell(1, numel(k_max));

    for ell = 0:numel(k_max)-1
        [k_max(ell+1), r0{ell+1}] = num_besselj_zeros(ell+(d-2)/2, N*pi/2);

        if k_max(ell+1) == 0
            ell_max = ell-1;

            k_max = k_max(1:ell_max+1);
            r0 = r0(1:ell_max+1);
            break;
        end
    end

    for ell = 0:ell_max
        r0{ell+1} = cat(1, r0{ell+1}', NaN(max(k_max)-k_max(ell+1), 1));
    end

    r0 = cell2mat(r0);

    basis = struct();

    basis.type = basis_type_ffb();

    basis.sz = sz;

    basis.domain = domain;

    if d == 2
        basis.count = k_max(1) + sum(2*k_max(2:end));
    elseif d == 3
        basis.count = sum(k_max.*(2*[0:ell_max]'+1));
    end

    basis.ell_max = ell_max;
    basis.k_max = k_max;
    basis.r0 = r0;

    basis.indices = fb_indices(basis);

    basis.precomp = ffb_precomp(basis);
end

function [k, r0] = num_besselj_zeros(ell, r)
    k = 4;

    r0 = besselj_zeros(ell, k);
    while all(r0 < r)
        k = 2*k;

        r0 = besselj_zeros(ell, k);
    end

    r0 = r0(r0 < r);

    k = numel(r0);
end

function nrm = basis_norm_2d(basis, ell, k)
    nrm = abs(besselj(ell+1, basis.r0(k,ell+1)))*sqrt(pi/2)* ...
        sqrt(prod(basis.sz/2));

    if ell == 0
        nrm = nrm*sqrt(2);
    end
end

function precomp = ffb_precomp(basis)
    d = numel(basis.sz);

    if d == 2
        precomp = ffb_precomp_2d(basis);
    elseif d == 3
        precomp = ffb_precomp_3d(basis);
    end
end

function precomp = ffb_precomp_2d(basis)
    n_r = basis.sz(1);

    [precomp.r, precomp.w] = lgwt(n_r, 0, 1/2);

    % Second dimension below is not basis.count because we're not counting
    % signs, since radial part is the same for both cos and sin.
    precomp.radial = zeros(n_r, sum(basis.k_max));

    ind = 1;

    for ell = 0:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        besselj_zeros = basis.r0(1:basis.k_max(ell+1), ell+1)';
        radial_ell = besselj(ell, 2*precomp.r*besselj_zeros);

        % NOTE: We need to remove the factor due to the discretization here
        % since it is already included in our quadrature weights.
        nrms = 1/(sqrt(prod(basis.sz)))*basis_norm_2d(basis, ell, [1:basis.k_max(ell+1)])';

        radial_ell = bsxfun(@times, radial_ell, 1./nrms);

        precomp.radial(:,idx) = radial_ell;

        ind = ind + numel(idx);
    end

    n_theta = 2*basis.sz(1);

    % Only calculate "positive" frequencies in one half-plane.
    freqs_x = precomp.r*cos([0:n_theta-1]*2*pi/(2*n_theta));
    freqs_y = precomp.r*sin([0:n_theta-1]*2*pi/(2*n_theta));

    precomp.freqs = cat(1, permute(freqs_x, [3 1 2]), permute(freqs_y, [3 1 2]));
end

function precomp = ffb_precomp_3d(basis)
    precomp = struct();

    n_r = basis.ell_max+1;
    n_theta = 2*basis.sz(1);
    n_phi = basis.ell_max+1;

    [r, wt_r] = lgwt(n_r, 0, 1/2);
    [z, wt_z] = lgwt(n_phi, -1, 1);
    phi = acos(z);
    wt_phi = wt_z;
    theta = 2*pi*[0:n_theta-1]'/(2*n_theta);

    precomp.radial_wtd = zeros([n_r max(basis.k_max) basis.ell_max+1]);
    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_ell = sph_bessel(ell, 2*r*basis.r0(1:k_max_ell,ell+1)');
        nrm = abs(sph_bessel(ell+1, basis.r0(1:k_max_ell,ell+1)')/4);
        radial_ell = bsxfun(@times, radial_ell, 1./nrm);

        radial_ell_wtd = bsxfun(@times, r.^2.*wt_r, radial_ell);

        precomp.radial_wtd(:,1:k_max_ell,ell+1) = radial_ell_wtd;
    end

    precomp.angular_phi_wtd_even = {};
    precomp.angular_phi_wtd_odd = {};

    for m = 0:basis.ell_max
        n_even_ell = floor((basis.ell_max-m)/2)+1 ...
            -mod(basis.ell_max, 2)*mod(m, 2);
        n_odd_ell = basis.ell_max-m+1-n_even_ell;

        angular_phi_wtd_m_even = zeros(n_phi, n_even_ell);
        angular_phi_wtd_m_odd = zeros(n_phi, n_odd_ell);

        ind_even = 1;
        ind_odd = 1;

        for ell = m:basis.ell_max
            angular_phi_m_ell = norm_assoc_legendre(ell, m, cos(phi));

            nrm_inv = sqrt(1/(2*pi));
            angular_phi_m_ell = angular_phi_m_ell*nrm_inv;

            angular_phi_wtd_m_ell = wt_phi.*angular_phi_m_ell;

            if mod(ell, 2) == 0
                angular_phi_wtd_m_even(:,ind_even) = ...
                    angular_phi_wtd_m_ell;
                ind_even = ind_even+1;
            else
                angular_phi_wtd_m_odd(:,ind_odd) = ...
                    angular_phi_wtd_m_ell;
                ind_odd = ind_odd+1;
            end
        end

        precomp.angular_phi_wtd_even{m+1} = angular_phi_wtd_m_even;
        precomp.angular_phi_wtd_odd{m+1} = angular_phi_wtd_m_odd;
    end

    angular_theta = zeros(n_theta, 2*basis.ell_max+1);

    angular_theta(:,1:basis.ell_max) = ...
        sqrt(2)*sin(theta*[basis.ell_max:-1:1]);
    angular_theta(:,basis.ell_max+1) = ones(n_theta, 1);
    angular_theta(:,basis.ell_max+2:2*basis.ell_max+1) = ...
        sqrt(2)*cos(theta*[1:basis.ell_max]);

    angular_theta_wtd = 2*pi/n_theta*angular_theta;

    precomp.angular_theta_wtd = angular_theta_wtd;

    [theta_grid, phi_grid, r_grid] = ndgrid(theta, phi, r);

    fourier_x = r_grid .* cos(theta_grid) .* sin(phi_grid);
    fourier_y = r_grid .* sin(theta_grid) .* sin(phi_grid);
    fourier_z = r_grid .* cos(phi_grid);

    precomp.fourier_pts = 2*pi*[fourier_x(:) fourier_y(:) fourier_z(:)]';
end

function indices = fb_indices(basis)
    indices = struct();

    d = numel(basis.sz);

    if d == 2
        indices.ells = zeros(basis.count, 1);
        indices.ks = zeros(basis.count, 1);
        indices.sgns = zeros(basis.count, 1);

        ind = 1;

        for ell = 0:basis.ell_max
            if ell == 0
                sgns = 1;
            else
                sgns = [1 -1];
            end

            ks = 0:basis.k_max(ell+1)-1;

            for sgn = sgns
                rng = ind:ind+numel(ks)-1;

                indices.ells(rng) = ell;
                indices.ks(rng) = ks;
                indices.sgns(rng) = sgn;

                ind = ind + numel(rng);
            end
        end
    elseif d == 3
        indices.ells = zeros(basis.count, 1);
        indices.ms = zeros(basis.count, 1);
        indices.ks = zeros(basis.count, 1);

        ind = 1;

        for ell = 0:basis.ell_max
            ks = 0:basis.k_max(ell+1)-1;
            ms = -ell:ell;

            for m = -ell:ell
                rng = ind:ind+numel(ks)-1;

                indices.ells(rng) = ell;
                indices.ms(rng) = m;
                indices.ks(rng) = ks;

                ind = ind+numel(ks);
            end
        end
    end
end

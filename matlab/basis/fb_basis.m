% FB_BASIS Construct a Fourier-Bessel basis object
%
% Usage
%    basis = fb_basis(sz, ell_max);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only square images and cubic volumes are supported.
%    ell_max: The maximum order ell of the basis elements. If set to Inf,
%       the basis includes all ell such that the resulting basis vectors
%       are concentrated below the Nyquist frequency (default Inf).
%
% Output
%    basis: A Fourier-Bessel basis object corresponding to the parameters.
%
% Description
%    The returned basis object can be used to project a vector x in the
%    standard basis onto the basis or to map a coefficient vector v into
%    the standard basis. This is achieved using the basis.expand and
%    basis.evaluate functions, respectively. The fields basis.sz denotes
%    the size of the vectors in the standard basis while basis.count denotes
%    the number of vectors in the basis.
%
%    For example, we can generate a random vector in the standard basis and
%    project it onto the Fourier-Bessel basis through
%
%       sz = 16*ones(1, 2);
%       basis = fb_basis(sz, floor(sz(1)/2));
%       x = randn(sz);
%       v = basis_expand(basis, x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis_evaluate(basis, v);
%
%    The adjoints of the expand and evaluate functions are obtained through
%    basis.expand_t and basis.evaluate_t functions, respectively.
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
%    All of the basis vectors are normalized so that their 2-norm, as continu-
%    ous functions, is equal to 1. Additionally, the basis functions are ortho-
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
%    cosine part followed by the sine part (except for ell = 0, in which case
%    we only have the cosine part). In 3D, we instead have the coefficients
%    arranged by order m, ranging from -ell to +ell.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function basis = fb_basis(sz, ell_max)
    if nargin < 2 || isempty(ell_max)
        ell_max = Inf;
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

    basis.type = basis_type_fb();

    basis.sz = sz;

    if d == 2
        basis.count = k_max(1) + sum(2*k_max(2:end));
    elseif d == 3
        basis.count = sum(k_max.*(2*[0:ell_max]'+1));
    end

    basis.ell_max = ell_max;
    basis.k_max = k_max;
    basis.r0 = r0;

    if d == 2
        basis.coords = unique_coordinates_2d(N);
    elseif d == 3
        basis.coords = unique_coordinates_3d(N);
    end

    basis.indices = fb_indices(basis);

    basis.precomp = fb_precomp(basis);

    if isempty(basis.precomp)
        basis = rmfield(basis, 'precomp');
    end

    basis.norms = fb_norms(basis);
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

function coords = unique_coordinates_2d(N)
    grid2d = grid_2d(N);

    mask = grid2d.r<=1;

    r = grid2d.r(mask);
    phi = grid2d.phi(mask);

    [r_unique, ~, r_idx] = unique(r);
    [ang_unique, ~, ang_idx] = unique(phi);

    coords = struct();

    coords.r_unique = r_unique;
    coords.ang_unique = ang_unique;
    coords.r_idx = r_idx;
    coords.ang_idx = ang_idx;
    coords.mask = mask;
end

function coords = unique_coordinates_3d(N)
    grid3d = grid_3d(N);

    mask = grid3d.r<=1;

    r = grid3d.r(mask);
    theta = grid3d.theta(mask);
    phi = grid3d.phi(mask);

    [r_unique, ~, r_idx] = unique(r);
    [ang_unique, ~, ang_idx] = unique([theta phi], 'rows');

    coords = struct();

    coords.r_unique = r_unique;
    coords.ang_unique = ang_unique;
    coords.r_idx = r_idx;
    coords.ang_idx = ang_idx;
    coords.mask = mask;
end

function norms = fb_norms(basis)
    d = numel(basis.sz);

    norms = zeros(sum(basis.k_max), 1);

    ind = 1;
    for ell = 0:basis.ell_max
        for k = 1:basis.k_max(ell+1)
            if d == 2
                norms(ind) = basis_norm_2d(basis, ell, k);
            elseif d == 3
                norms(ind) = basis_norm_3d(basis, ell, k);
            end

            ind = ind+1;
        end
    end
end

function nrm = basis_norm_2d(basis, ell, k)
    nrm = abs(besselj(ell+1, basis.r0(k,ell+1)))*sqrt(pi/2)* ...
        sqrt(prod(basis.sz/2));

    if ell == 0
        nrm = nrm*sqrt(2);
    end
end

function nrm = basis_norm_3d(basis, ell, k)
    nrm = abs(sph_bessel(ell+1, basis.r0(k,ell+1)))/sqrt(2)* ...
        sqrt(prod(basis.sz/2));
end

function precomp = fb_precomp(basis)
    d = numel(basis.sz);

    if d == 2
        r_unique = basis.coords.r_unique;
        ang_unique = basis.coords.ang_unique;
        r_idx = basis.coords.r_idx;
        ang_idx = basis.coords.ang_idx;
        mask = basis.coords.mask;

        ind_radial = 1;
        ind_ang = 1;

        radial = zeros(size(r_unique, 1), sum(basis.k_max));
        ang = zeros(size(ang_unique, 1), 1+2*(basis.ell_max-1));

        for ell = 0:basis.ell_max
            for k = 1:basis.k_max(ell+1)
                radial(:,ind_radial) = ...
                    besselj(ell, basis.r0(k,ell+1)*r_unique);

                ind_radial = ind_radial+1;
            end

            if ell == 0
                sgns = 1;
            else
                sgns = [1 -1];
            end

            for sgn = sgns
                if sgn == 1
                    ang(:,ind_ang) = cos(ell*ang_unique);
                else
                    ang(:,ind_ang) = sin(ell*ang_unique);
                end

                ind_ang = ind_ang+1;
            end
        end
    elseif d == 3
        r_unique = basis.coords.r_unique;
        ang_unique = basis.coords.ang_unique;
        r_idx = basis.coords.r_idx;
        ang_idx = basis.coords.ang_idx;
        mask = basis.coords.mask;

        ind_radial = 1;
        ind_ang = 1;

        radial = zeros(size(r_unique, 1), sum(basis.k_max));
        ang = zeros(size(ang_unique, 1), sum(2:[0:basis.ell_max]+1));

        for ell = 0:basis.ell_max
            for k = 1:basis.k_max(ell+1)
                radial(:,ind_radial) = ...
                    sph_bessel(ell, basis.r0(k,ell+1)*r_unique);

                ind_radial = ind_radial+1;
            end

            for m = -ell:ell
                ang(:,ind_ang) = real_sph_harmonic(ell, m, ...
                    ang_unique(:,1), ang_unique(:,2));

                ind_ang = ind_ang+1;
            end
        end
    end

    precomp = struct();
    precomp.radial = radial;
    precomp.ang = ang;
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

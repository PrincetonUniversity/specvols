% FB_BASIS_EVALUATE_T Evaluate coefficient in dual Fourier-Bessel basis
%
% Usage
%    v = fb_basis_evaluate_t(basis, x);
%
% Input/Output
%    See documentation for `basis_evaluate_t`.
%
% See also
%    basis_evaluate_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = fb_basis_evaluate_t(basis, x)
    d = numel(basis.sz);

    if d == 2
        v = fb_basis_evaluate_t_2d(basis, x);
    elseif d == 3
        v = fb_basis_evaluate_t_3d(basis, x);
    end
end

function v = fb_basis_evaluate_t_2d(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    r_unique = basis.coords.r_unique;
    ang_unique = basis.coords.ang_unique;
    r_idx = basis.coords.r_idx;
    ang_idx = basis.coords.ang_idx;
    mask = basis.coords.mask;

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    v = zeros([basis.count size(x, 2)], class(x));

    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = basis.norms(idx_radial);

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = besselj(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        if ell == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            if ~is_precomp
                if sgn == 1
                    ang = cos(ell*ang_unique);
                else
                    ang = sin(ell*ang_unique);
                end
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            v(idx,:) = ang_radial'*x(mask,:);

            ind = ind + numel(idx);

            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    v = roll_dim(v, sz_roll);
end

function v = fb_basis_evaluate_t_3d(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    r_unique = basis.coords.r_unique;
    ang_unique = basis.coords.ang_unique;
    r_idx = basis.coords.r_idx;
    ang_idx = basis.coords.ang_idx;
    mask = basis.coords.mask;

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    v = zeros([basis.count size(x, 2)], class(x));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = basis.norms(idx_radial);

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = sph_bessel(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        for m = -ell:ell
            if ~is_precomp
                ang = real_sph_harmonic(ell, m, ang_unique(:,1), ...
                    ang_unique(:,2));
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            v(idx,:) = ang_radial'*x(mask,:);

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    v = roll_dim(v, sz_roll);
end

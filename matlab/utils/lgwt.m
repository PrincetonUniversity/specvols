% LGWT Compute Legendre-Gauss quadrature
%
% Usage
%    [x, w] = lgwt(N, a, b);
%
% Input
%    N: The truncation order, that is, the number of nodes.
%    a, b: The endpoints of the interval over which the quadrature is defined.
%
% Output
%    x: The quadrature nodes.
%    w: The quadrature weights.
%
% Description
%    This script is for computing definite integrals using Legendre-Gauss
%    quadrature. Computes the Legendre-Gauss nodes and weights on an interval
%    [a, b] with truncation order N.
%
%    Suppose you have a continuous function f(x) which is defined on [a, b]
%    which you can evaluate at any x in [a, b]. Simply evaluate it at all of
%    the values contained in the x vector to obtain a vector f, then compute
%    the definite integral using sum(f.*w);
%
%    The nodes are sorted in ascending order.

% Authors
%    Greg von Winckel
%    Zhizhen Zhao <zhizhenz@illinois.edu>
%    Joakim Anden <janden@flatironinstitute.org>

function [x, w] = lgwt(N, a, b)
    N = N-1;
    N1 = N+1;
    N2 = N+2;

    xu = linspace(-1, 1, N1)';

    % Initial guess.
    y = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

    % Legendre-Gauss Vandermonde matrix.
    L = zeros(N1, N2);

    % Derivative of the matrix.
    Lp = zeros(N1, N2);

    % Compute the zeros of the (N+1)th Legendre polynomial using the recursion
    % relation and the Newton-Raphson method. Iterate until new points are
    % uniformly within epsilon of old points.
    y0 = 2;
    while max(abs(y-y0)) > eps
        L(:,1) = 1;
        Lp(:,1) = 0;

        L(:,2) = y;
        Lp(:,2) = 1;

        for k = 2:N1
            L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
        end

        Lp = N2*(L(:,N1)-y.*L(:,N2))./(1-y.^2);

        y0 = y;
        y = y0-L(:,N2)./Lp;
    end

    % Linear map from [-1, 1] to [a, b].
    x = (a*(1-y)+b*(1+y))/2;

    % Compute the weights.
    w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

    % Sort x in ascending order.
    [x, id] = sort(x, 'ascend');
    w = w(id);
end

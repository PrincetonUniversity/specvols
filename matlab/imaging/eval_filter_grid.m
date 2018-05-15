% EVAL_FILTER_GRID Evaluates a filter on a grid
%
% Usage
%    h = eval_filter_grid(filter, L);
%
% Input
%    filter: A filter object (see `eval_filter` for more information). This
%       could be an array of filters, in which case all filters in the array
%       are evaluated.
%    L: The dimension of the grid.
%
% Output
%    h: An L-by-L-by-K array containing the values of all K filter on the fre-
%       quency grid.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function h = eval_filter_grid(filter, L)
    grid2d = grid_2d(L);

    omega = pi*[grid2d.x(:) grid2d.y(:)]';

    h = eval_filter(filter, omega);

    h = reshape(h, [size(grid2d.x) numel(filter)]);
end

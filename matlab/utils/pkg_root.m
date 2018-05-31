% PKG_ROOT Obtain root directory of the package
%
% Usage
%    root = pkg_root();
%
% Output
%    root: A string containing the path to the package.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function root = pkg_root()
    full_path = mfilename('fullpath');

    [root, ~, ~] = fileparts(full_path);

    ind = find(root == filesep, 2, 'last');
    ind = ind(1);

    root = root(1:ind-1);
end

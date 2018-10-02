% STAR_TO_SRC Create a source from a STAR file
%
% Usage
%    src = star_to_src(star, file_root, star_opt);
%
% Input
%    star: A STAR file struct obtained from `load_star`.
%    file_root: The root directory of any relative paths in the STAR file
%       (default '.').
%    star_opt: An options structure with the fields:
%          - strip_path: if true, strips the paths in the STAR files, only
%             keeping the filenames (default false).
%
% Output
%    src: A source object corresponding to the given STAR file.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = star_to_src(star, file_root, star_opt)
    if nargin < 2 || isempty(file_root)
        file_root = '.';
    end

    if nargin < 3 || isempty(star_opt)
        star_opt = struct();
    end

    star_opt = fill_struct(star_opt, ...
        'strip_path', false);

    src = struct();

    src.type = src_type_mrcs();

    if ~isstruct(star)
        error(['Input `star` must be a STAR file struct obtained from ' ...
            '`load_star`']);
    end

    if isfield(star, 'root')
        star = star.root;
    elseif isfield(star, 'images')
        star = star.images;
    elseif ~isfield(star, 'rlnImageName')
        error('Invalid STAR file format.');
    end

    src.L = [];
    src.n = numel(star);

    [file_names, file_idx, image_idx] = star_image_files(star, file_root, ...
        star_opt.strip_path);

    src.file_names = file_names;
    src.file_idx = file_idx;
    src.image_idx = image_idx;

    if ~exist(src.file_names{1}, 'file')
        warning(sprintf('File ''%s'' does not exist.', src.file_names{1}));

        src.L = [];
        src.precision = [];
    else
        mrc = mrc_open(src.file_names{1});

        if mrc.header.N(2) ~= mrc.header.N(1)
            error('Only square images are supported.');
        end

        src.L = mrc.header.N(1);

        src.precision = mrc.data_type;

        mrc_close(mrc);
    end

    src.params = star_to_params(star, star_opt);
end

function [file_names, file_idx, image_idx] = star_image_files(star, ...
    file_root, strip_path)

    file_names = cell(numel(star), 1);
    image_idx = zeros(numel(star), 1);

    for s = 1:numel(star)
        image_ref = star(s).rlnImageName;

        if ~isempty(strfind(image_ref, '@'))
            % Part of a stack.

            file_parts = fast_split(image_ref, {'@'});

            file_names{s} = file_parts{2};
            image_idx(s) = sscanf(file_parts{1}, '%d');
        else
            % Just an image file.

            file_names{s} = image_ref;
            image_idx(s) = 1;
        end

        if strip_path
            [~, name, ext] = fileparts(file_names{s});
            file_names{s} = [name ext];
        end

        if ~is_absolute_pathname(file_names{s})
            file_names{s} = fullfile(file_root, file_names{s});
        end
    end

    [file_names, ~, file_idx] = unique(file_names);
end

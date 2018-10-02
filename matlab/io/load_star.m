% LOAD_STAR Read STAR file
%
% Usage
%    data = load_star(filename);
%
% Input
%    filename: The filename of the STAR file to read.
%
% Output
%    data: A struct containing all the data blocks in the STAR file. Each data
%       block corresponds to a field. In turn, each data block is itself a
%       struct (in the case of a list) or a struct array (in the case of a
%       loop).

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function data = load_star(filename)
    fid = fopen(filename, 'r');

    % Name of the current data block.
    name = [];
    % The output struct containing all data blocks.
    data = struct();

    while true
        ch = fpeek(fid);
        if ~ischar(ch)
            % Reached EOF
            break;
        end

        if ch == '_'
            % Found a list. Try to read it and append to output struct.
            if isempty(name)
                error('list outside of data block');
            end

            star_list = load_star_struct(fid);

            data = setfield(data, name, star_list);
        end

        line = fgetl(fid);

        if ~isempty(strfind(line, 'data_'))
            % Found the start of a data block. Update current name.
            name = line(strfind(line, 'data_') + 5:end);
            if isempty(name)
                name = 'root';
            end
        end

        if ~isempty(strfind(line, 'loop_'))
            % Found a loop. Try to read it and append to output struct.
            if isempty(name)
                error('loop outside of data block');
            end

            star_loop = load_star_loop(fid);

            data = setfield(data, name, star_loop);
        end
    end

    fclose(fid);
end

function star_struct = load_star_struct(fid)
    [labels, values] = load_star_list(fid);

    % Interleave labels and values for call to struct.
    labels_values = [labels; values];
    labels_values = labels_values(:);

    star_struct = struct(labels_values{:});
end

function star_loop = load_star_loop(fid)
    % Note: load_star_list will only read labels here since it's a loop
    % header.
    labels = load_star_list(fid);
    values = load_star_data(fid, length(labels));

    % Interleave labels and values for call to struct.
    labels_values = {};
    for n = 1:length(labels)
        labels_values{2*(n-1)+1} = labels{n};
        labels_values{2*(n-1)+2} = values(:,n);
    end

    star_loop = struct(labels_values{:});
end

function [labels, values] = load_star_list(fid)
    n = 1;
    labels = {};
    values = {};

    while true
        ch = fpeek(fid);

        if isempty(ch)
            % Reached EOF (this shouldn't happen unless loop is empty).
            break;
        end

        if logical(ch == '#') || logical(ch == char(0)) || logical(ch == ';')
            % Comment line. Skip.
            fgetl(fid);
            continue;
        end

        if ch ~= '_'
            % This is tne end of the header.
            break;
        end

        line = fgetl(fid);

        % Strip out trailing comments.
        if ~isempty(strfind(line, '#'))
            line = line(1:strfind(line, '#')-1);
        end

        parts = strsplit(line, {' ', '\t'});

        % First character is '_' so skip it.
        labels{n} = parts{1}(2:end);
        values{n} = parts{2};

        n = n+1;
    end

    % Make sure all data has correct type.
    values = cellfun(@cast_data, values, 'UniformOutput', false);
end

function values = load_star_data(fid, num_fields)
    k = 1;
    values = {};

    while true
        line = fgetl(fid);

        if ~ischar(line)
            % Reached EOF
            break;
        end

        if ~isempty(find(line=='#', 1))
            ind = find(line=='#', 1);
            line = line(1:ind-1);
        end

        line = strtrim(line);

        if isempty(line)
            % Reached end of loop
            break;
        end

        %line_values = strsplit(line, {' ', '\t'});
        line_values = fast_split(line, {' ', sprintf('\t')});

        if length(line_values) > num_fields
            warning('too many values on line; truncating');
            line_values = line_values(1:num_fields);
        elseif length(line_values) < num_fields
            warning('too few values on line; padding with empty');
            [line_values{end+1:num_fields}] = deal([]);
        end

        % Cast data to correct type.
        line_values = cellfun(@cast_data, line_values, 'UniformOutput', false);

        values(k,:) = line_values;

        k = k+1;
    end
end

function data = cast_data(data)
    % Don't use ismember here. It's too slow.
    if ~isempty(data) && is_double(data)
        % Cast to double.
        data = sscanf(data, '%g');
        return;
    elseif strcmp(data, 'None')
        % Replace 'None' string with empty value.
        data = [];
        return;
    end
end

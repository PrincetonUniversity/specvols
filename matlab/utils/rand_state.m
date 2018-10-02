% RAND_STATE Gets and sets PRNG sates
%
% Usage
%    old_state = rand_state();
%    old_state = rand_state(new_state);
%
% Input
%    new_state: The new state vector to use for the PRNGs. This is either a
%       single column vector or an array with five columns, one for each
%       PRNG. If a single column is given, the same state vector is used
%       for all PRNGs. Each column is used as input for the rand*('state')
%       call, and so can be a state vector or an arbitrary vector (or scalar
%       from which a hash is calculated. If missing or empty, the state
%       vectors of the PRNGs are left unchanged.
%
% Output
%    old_state: The old state vectors arranged as columns in an array of size
%       625-by-5.
%
% Description
%    This function gets and/or sets the state vectors of the common PRNG
%    functions 'rand', 'randn', 'rande', 'randg', and 'randp', each of which
%    corresponds to a column in a 625-by-5 array, in that order. This can
%    be used to set all the PRNG states at once, save the old state, and then
%    restore that state at some later time. This can be done most easily using
%    the `rand_push` and `rand_pop` PRNG stack functions.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function old_state = rand_state(new_state)
    if nargin < 1
        new_state = [];
    end

    if isoctave()
        old_state = zeros(625, 5);

        old_state(:,1) = rand('state');
        old_state(:,2) = randn('state');
        old_state(:,3) = rande('state');
        old_state(:,4) = randg('state');
        old_state(:,5) = randp('state');

        if ~isempty(new_state)
            if size(new_state, 2) == 1
                new_state = repmat(new_state, 1, 5);
            elseif size(new_state, 2) ~= 5
                error('new_state must have either 1 or 5 columns');
            end

            rand('state', new_state(:,1));
            randn('state', new_state(:,2));
            rande('state', new_state(:,3));
            randg('state', new_state(:,4));
            randp('state', new_state(:,5));
        end
    else
        old_state = rng;

        if ~isempty(new_state)
            if isnumeric(new_state) && numel(new_state) == 1
                rng(new_state);
            elseif isstruct(new_state) && isfield(new_state, 'Type') && ...
                isfield(new_state, 'Seed') && isfield(new_state, 'State')
                rng(new_state);
            else
                error('new_state must be a scalar or a RNG state structure');
            end
        end
    end
end

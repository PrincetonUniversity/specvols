% RAND_STACK Utility function for managing PRNG state stack
%
% Usage
%    out = rand_stack(action, in);
%
% Input
%    action: An action code. One of
%          - 0: Push onto stack. `in` is the input to the `rand_state`
%             function.
%          - 1: Pop stack.
%          - 2: Extract stack. `out` returns the whole stack as a three-
%             dimensional array.
%    in: Depends on action code.
%
% Output
%    out: Depends on action code.
%
% Note
%    This function should not be called directly. Instead, use the `rand_push`
%    and `rand_pop` functions.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function out = rand_stack(action, in)
    persistent p_state_stack;

    if nargin < 2
        in = [];
    end

    if isempty(p_state_stack)
        p_state_stack = {};
    end

    if action == 0
        % Push onto stack.

        old_state = rand_state(in);

        p_state_stack = cat(1, p_state_stack, old_state);
    elseif action == 1
        % Pop stack.

        if numel(p_state_stack) == 0
            error('PRNG state stack is empty.');
        end

        rand_state(p_state_stack{end});

        p_state_stack = p_state_stack(1:end-1);
    elseif action == 2
        % Extract stack.

        out = p_state_stack;
    else
        error('Invalid action code.');
    end
end

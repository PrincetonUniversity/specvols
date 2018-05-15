% RAND_PUSH Push a new state onto the PRNG state stack
%
% Usage
%    rand_push(new_state);
%
% Input
%    new_state: The new state array for the PRNGs. For its definition, see
%       `rand_state` (default zeros(1, 5)).
%
% Description
%    Pushes a new state onto the stack of PRNG states. Combined with
%    `rand_pop`, this allows for initializing the PRNGs, generating some
%    random numbers, and then restoring the previous PRNG states. While
%    this can be done manually using the `rand_state` function, this requires
%    storing the old PRNG states, which is not necessary using the
%    `rand_push` and `rand_pop` functions. The following code illustrates its
%    use:
%
%       rand_push();
%
%       % Do work using 'rand', 'rand', etc.
%
%       rand_pop();
%
%    For more information on the PRNGs affected by these functions, see the
%    documentation for `rand_state`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function rand_push(new_state)
    if nargin < 1 || isempty(new_state)
        new_state = 0;
    end

    rand_stack(0, new_state);
end

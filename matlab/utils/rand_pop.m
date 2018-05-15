% RAND_POP Pops a state off the PRNG state stack
%
% Usage
%    rand_pop();
%
% Description
%    Pops a state off the stack of PRNG states. Combined with `rand_push`,
%    this allows for initializing the PRNGs, generating some random numbers,
%    and then restoring the previous PRNG states. While this can be done
%    manually using the `rand_state` function, this requires storing the old
%    PRNG states, which is not necessary using the `rand_push` and `rand_pop`
%    functions. The following code illustrates its use:
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

function rand_pop()
    rand_stack(1);
end

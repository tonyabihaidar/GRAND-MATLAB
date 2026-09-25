% =========================================================================
% bin_GRAND.m
%
% PURPOSE:
%   Hard-decision GRAND (Guessing Random Additive Noise Decoding) for
%   binary linear codes.
%
% ALGORITHM OVERVIEW:
%   GRAND is a universal decoder that, instead of searching the codebook,
%   searches the noise space.  It guesses noise patterns e in order of
%   decreasing likelihood and tests each guess by checking whether
%   (y XOR e) is a valid codeword — i.e., whether H*(y XOR e)' = 0,
%   or equivalently, whether H*e' = H*y' (mod 2).
%
%   Hard GRAND has no soft information: it uses only the hard-decoded
%   bits y_demod in {0,1}^n and queries noise patterns in increasing
%   Hamming weight order (weight 0, then 1, then 2, ...).  This is the
%   maximum-likelihood ordering for a BSC (Binary Symmetric Channel),
%   and serves as the baseline decoder in Cauchy-noise comparisons.
%
%   For each Hamming weight w, all C(n,w) patterns are enumerated in
%   lexicographic order of their support positions.
%
% COMPLEXITY:
%   Expected number of queries (EG) grows with noise level.  At weight
%   w, there are C(n,w) patterns.  The decoder stops at the FIRST valid
%   codeword found (i.e., the closest codeword in Hamming distance).
%   For a [128,116] RLC with 2^12 = 4096 possible syndromes, EG is
%   typically a few thousand at moderate noise.
%
% INPUTS:
%   H         : (n-k) x n binary parity-check matrix
%   max_query : maximum number of syndrome checks before giving up
%   y_demod   : 1 x n binary vector of hard-decided received bits {0,1}
%
% OUTPUTS:
%   y_decoded      : 1 x n decoded codeword, or -ones(1,n) if abandoned
%   putative_noise : 1 x n estimated noise pattern (zeros if abandoned)
%   n_guesses      : total number of syndrome checks performed
%   abandoned      : 1 if max_query was reached before finding a codeword
% =========================================================================

function [y_decoded, putative_noise, n_guesses, abandoned] = bin_GRAND(H, max_query, y_demod)

    y_demod = double(y_demod(:).');   % ensure a double-precision row vector
    n = size(H, 2);

    if length(y_demod) ~= n
        error('Length of y_demod must equal number of columns of H.');
    end

    n_guesses  = 0;
    abandoned  = 0;
    y_decoded  = -ones(1, n);    % sentinel: overwritten on success
    putative_noise = zeros(1, n);

    % Pre-compute the target syndrome: H*y mod 2.
    % A noise pattern e is valid iff H*e' = Hy (mod 2),
    % since H*(y XOR e)' = H*y' XOR H*e' = 0 requires H*e' = H*y'.
    Hy = mod(H * y_demod.', 2);

    % ------------------------------------------------------------------
    % Query 1: weight-0 noise pattern (no errors assumed).
    % Equivalent to checking whether y_demod itself is a codeword.
    % ------------------------------------------------------------------
    n_guesses = n_guesses + 1;
    if all(Hy == 0)
        % y_demod is already a codeword — no channel errors detected.
        putative_noise = zeros(1, n);
        y_decoded      = y_demod;
        return;
    end

    if n_guesses >= max_query
        abandoned = 1;
        return;
    end

    % ------------------------------------------------------------------
    % Weights w = 1, 2, ..., n: enumerate all C(n,w) patterns per weight.
    % Patterns are represented by their support: a sorted list of w
    % bit positions (1-indexed) to flip.
    % ------------------------------------------------------------------
    for w = 1:n

        combo = 1:w;   % lexicographically first combination: positions {1, 2, ..., w}

        while true

            n_guesses = n_guesses + 1;

            % Syndrome of the candidate noise pattern e with 1s at `combo`.
            % H*e' = sum of columns of H at the flipped positions, mod 2.
            e_syn = mod(sum(H(:, combo), 2), 2);

            if isequal(e_syn, Hy)
                % Valid noise pattern found: (y XOR e) is a codeword.
                putative_noise           = zeros(1, n);
                putative_noise(combo)    = 1;
                y_decoded                = mod(y_demod - putative_noise, 2);
                abandoned                = 0;
                return;
            end

            if n_guesses >= max_query
                abandoned = 1;
                return;
            end

            % Advance to the next combination in lexicographic order.
            combo = next_combination(combo, n);
            if isempty(combo)
                break;   % all C(n,w) patterns at this weight exhausted
            end

        end % while (patterns at weight w)

    end % for w

    % All patterns up to weight n exhausted without a match.
    abandoned = 1;

end


% -------------------------------------------------------------------------
% next_combination(combo, n)
%
% Returns the next k-combination of {1,...,n} after `combo` in
% lexicographic order, or [] if `combo` is the last combination.
%
% Algorithm: find the rightmost position i that can be incremented
% (i.e., combo(i) < n - w + i), then increment it and reset all
% positions to the right to the smallest valid values.
% -------------------------------------------------------------------------
function combo = next_combination(combo, n)

    w = length(combo);
    i = w;

    % Walk left until we find a position that has room to increment.
    while i >= 1 && combo(i) == n - w + i
        i = i - 1;
    end

    if i == 0
        % All positions are at their maximum values: this was the last combo.
        combo = [];
        return;
    end

    % Increment position i and set positions i+1 ... w to consecutive values.
    combo(i) = combo(i) + 1;
    for j = i+1 : w
        combo(j) = combo(j-1) + 1;
    end

end
% =========================================================================
% bin_ORBGRAND.m
%
% PURPOSE:
%   Soft-decision ORBGRAND (Ordered Reliability Bits GRAND) decoder for
%   binary linear codes.
%
% ALGORITHM OVERVIEW:
%   ORBGRAND extends hard GRAND by exploiting soft channel information.
%   Instead of querying noise patterns in Hamming weight order, it uses
%   bit reliability scores (|LLR| values) to define a "logistic weight"
%   ordering that queries the most likely noise patterns first.
%
%   SORTING AND RANKING:
%     Each received bit is assigned a reliability score = |LLR|.
%     Bits are sorted in ASCENDING reliability order (least reliable first).
%     The resulting rank i in {1,...,n} is assigned to each bit, with
%     rank 1 = least reliable (most likely to be in error).
%
%   LOGISTIC WEIGHT:
%     A noise pattern that flips bit ranks {i_1, ..., i_w} has logistic
%     weight  W = sum(i_j).  Patterns with small W flip low-ranked
%     (unreliable) bits, which are the most probable errors.
%     W ranges from w*(w+1)/2 (flip ranks 1..w) to w*n - w*(w-1)/2
%     (flip ranks n-w+1..n).
%
%   QUERY ORDER:
%     Patterns are iterated over increasing W = 1, 2, 3, ...
%     For each W, all (w, rank-set) pairs consistent with that W are
%     enumerated via the Landslide algorithm.
%     This produces a near-ML ordering under any channel model for
%     which reliability ranks correlate with error probability.
%
%   SYNDROME CHECK:
%     To avoid re-permuting at each step, H is pre-reordered as
%     test_H = H(:, ind_order), so that syndrome checks on sorted-domain
%     error vectors are equivalent to syndrome checks in the original domain.
%
% INPUTS:
%   H         : (n-k) x n binary parity-check matrix
%   max_query : maximum number of syndrome checks before giving up
%   y_soft    : 1 x n signed soft LLR vector.
%               sign(y_soft(i)) gives the hard decision for bit i:
%                 y_soft(i) < 0  =>  hard bit = 1
%                 y_soft(i) > 0  =>  hard bit = 0
%               |y_soft(i)| is the reliability score for bit i.
%
%   NOTE: y_soft may be the exact Cauchy LLR or its 3-piece approximation.
%         In either case, sign(y_soft) must match the hard decision, which
%         is enforced by an assertion in the driver.
%
% OUTPUTS:
%   y_decoded : 1 x n decoded codeword, or -ones(1,n) if abandoned
%   err_vec   : 1 x n estimated noise pattern in original bit order
%   n_guesses : total number of syndrome checks performed
%   abandoned : 1 if max_query reached before finding a valid codeword
% =========================================================================

function [y_decoded, err_vec, n_guesses, abandoned] = bin_ORBGRAND(H, max_query, y_soft)

    y_soft = y_soft(:).';   % ensure row vector
    n      = length(y_soft);

    if size(H, 2) ~= n
        error('Length of y_soft must equal number of columns of H.');
    end

    % Hard decision from the sign of the soft LLR:
    %   LLR < 0  =>  bit = 1  (received sample in the negative half-plane)
    %   LLR > 0  =>  bit = 0
    y_demod = double(y_soft < 0);

    abandoned  = 0;
    y_decoded  = -ones(1, n);   % sentinel: overwritten on success
    err_vec    = zeros(1, n);

    % Pre-compute syndrome H*y mod 2.  A noise pattern e is valid iff
    % H*e' = Hy (mod 2), since H*(y-e)' = 0  <=>  H*e' = H*y'.
    Hy = mod(H * y_demod.', 2);

    % ------------------------------------------------------------------
    % Query 1: zero noise pattern — check if hard decisions form a codeword.
    % ------------------------------------------------------------------
    n_guesses = 1;
    if all(Hy == 0)
        y_decoded = y_demod;
        err_vec   = zeros(1, n);
        return;
    end

    if n_guesses >= max_query
        abandoned = 1;
        return;
    end

    % ------------------------------------------------------------------
    % Sort bits from LEAST to MOST reliable (ascending |LLR|).
    % ind_order(j) = original index of the j-th least reliable bit.
    % ind_order(1) = least reliable (smallest |LLR|, rank 1).
    % ind_order(n) = most  reliable (largest  |LLR|, rank n).
    %
    % This rank assignment means logistic weight W = sum(ranks of flipped
    % bits) is minimized by flipping the least-reliable bits first — the
    % correct heuristic for soft GRAND decoding.
    % ------------------------------------------------------------------
    reliability = abs(y_soft);
    [~, ind_order] = sort(reliability);   % ascending sort

    % Pre-permute H columns to match sorted bit order.
    % This avoids re-indexing at every syndrome check:
    %   test_H * err_ordered' = H(:,ind_order) * err_ordered'
    %                         = H * err_vec'      (since err_vec(ind_order) = err_ordered)
    test_H = H(:, ind_order);

    % ------------------------------------------------------------------
    % Main loop: iterate over logistic weight W = 1, 2, 3, ...
    % For each W, iterate over all Hamming weights w consistent with W,
    % then enumerate all rank-sets of size w summing to W via landslide().
    % ------------------------------------------------------------------
    Wmax = n * (n + 1) / 2;   % maximum possible logistic weight

    for W = 1 : Wmax

        if n_guesses >= max_query
            abandoned = 1;
            err_vec   = zeros(1, n);
            return;
        end

        % Maximum Hamming weight w such that w*(w+1)/2 <= W.
        % (Cannot flip more than this many bits at logistic weight W.)
        w_max_hamming = floor((sqrt(1 + 8*W) - 1) / 2);
        w_max_hamming = min(w_max_hamming, n);

        for w = 1 : w_max_hamming

            % Skip if W exceeds the largest achievable logistic weight
            % for w flips in ranks {1,...,n}: W_max = w*n - w*(w-1)/2.
            W_max_for_w = w*n - w*(w-1)/2;
            if W > W_max_for_w
                continue;
            end

            % Skip if W is below the smallest achievable logistic weight
            % for w flips: W_min = 1+2+...+w = w*(w+1)/2.
            if W < w*(w+1)/2
                continue;
            end

            % Enumerate all size-w subsets of {1,...,n} with rank sum = W.
            % landslide() returns these as rows of a matrix; each row is
            % a sorted list of w rank indices (1-based) in the SORTED domain.
            noise_locations = landslide(W, w, n);

            for jj = 1 : size(noise_locations, 1)

                if n_guesses >= max_query
                    abandoned = 1;
                    err_vec   = zeros(1, n);
                    return;
                end

                n_guesses = n_guesses + 1;

                % Construct the noise pattern in sorted-bit domain:
                % flip the bits at the w rank positions returned by landslide.
                err_ordered = zeros(1, n);
                err_ordered(noise_locations(jj, :)) = 1;

                % Syndrome check in sorted domain (equivalent to original domain
                % because test_H = H(:, ind_order)).
                test_syndrome = mod(test_H * err_ordered.', 2);

                if isequal(test_syndrome, Hy)

                    % Valid codeword found.  Map sorted-domain error back to
                    % original bit positions: err_vec(ind_order) = err_ordered.
                    err_vec            = zeros(1, n);
                    err_vec(ind_order) = err_ordered;

                    % Subtract estimated noise from hard decisions (mod 2).
                    y_decoded = mod(y_demod - err_vec, 2);
                    abandoned = 0;
                    return;

                end

            end % jj (patterns for this W, w)
        end % w (Hamming weight)
    end % W (logistic weight)

    % All logistic weights exhausted without finding a codeword.
    abandoned = 1;
    err_vec   = zeros(1, n);

end
function [y_soft, t1, yc, coeffs] = cauchy_llr_approx(r, gamma)
% =========================================================================
% cauchy_llr_approx.m
%
% PURPOSE:
%   Computes a 3-piece piecewise-linear approximation of the signed
%   Cauchy log-likelihood ratio (LLR) for BPSK signalling over a Cauchy
%   noise channel.  The approximation is designed to preserve the
%   RANK ORDER of bit reliabilities (i.e., |approx LLR| ranks bits in
%   the same order as |exact LLR|), which is the only property required
%   for optimal ORBGRAND decoding.
%
% BACKGROUND — EXACT CAUCHY LLR:
%   Under BPSK (bit 0 -> +1, bit 1 -> -1) with additive Cauchy noise
%   of scale parameter gamma > 0, the received sample is r = x + z,
%   where z ~ Cauchy(0, gamma).  The exact signed LLR is:
%
%     L(r) = log f(r | bit=0) / f(r | bit=1)
%           = log(gamma^2 + (r+1)^2) - log(gamma^2 + (r-1)^2)
%
%   Key properties of the reliability R(y) = |L(y)| for y = |r| >= 0:
%     - R(0) = 0  (received signal at the decision boundary is unreliable)
%     - R rises to its peak R(t1) at  t1 = sqrt(1 + gamma^2)
%     - R decays monotonically back to 0 as y -> infinity
%   The decay past the peak arises because extreme received values are
%   LESS reliable under Cauchy noise — large |r| could be explained by
%   a large noise outlier on either transmitted symbol.
%
% APPROXIMATION STRUCTURE — 3-PIECE PIECEWISE LINEAR:
%   On the positive half-axis y = |r| >= 0, R(y) is approximated as:
%
%     R_approx(y) = c1*y + c2*max(0, y-t1) + c3*max(0, y-yc)
%
%   This hinge-function form encodes three linear segments:
%
%     Segment 1  [0,  t1]:   slope = c1              > 0  (rising)
%     Segment 2  [t1, yc]:   slope = c1 + c2         < 0  (steep descent)
%     Segment 3  [yc, inf):  slope = c1 + c2 + c3    < 0  (gentle descent)
%
% BREAKPOINTS:
%   t1 = sqrt(1 + gamma^2)
%       Analytically exact location of the reliability peak (dR/dy = 0).
%
%   yc = last y > t1 where |dR/dy| = 0.1 * |dR/dy|_{y=0+}
%       The "crossover" point where the tail of R(y) has flattened to
%       10% of its initial rising slope.  This is the handoff between
%       the steep Segment 2 and the gentle Segment 3.
%
%       IMPORTANT — WHY THE *LAST* CROSSING:
%       After the peak, |R'(y)| first rises from near-zero to a local
%       maximum, then falls back to zero.  The equation |R'(y)| = target
%       therefore has TWO solutions past t1.  The FIRST solution is very
%       close to t1 (on the rising flank of |R'|) and would collapse
%       Segment 2 to near-zero width, making the approximation nearly
%       flat — which is wrong.  The LAST solution (on the falling flank)
%       is the physically correct tail breakpoint where the descent
%       genuinely slows.  We always use cross_idx(end).
%
% COEFFICIENTS (analytically computed from breakpoints):
%   c1 = R(t1) / t1              [slope of Segment 1: line through origin]
%   c2 = (R(yc)-R(t1))/(yc-t1) - c1   [hinge at t1, always negative]
%   c3 = R'(1.5*yc) - (c1+c2)   [hinge at yc, reduces descent rate]
%       If |slope3| > |slope2|, c3 is clamped to 0 so the tail cannot
%       steepen unphysically (R must be monotone non-increasing for y>t1).
%
% SIGNED OUTPUT CONVENTION:
%   y_soft = sign(r) * max(0, R_approx(|r|))
%   This ensures sign(y_soft) = sign(r) = sign of hard decision, as
%   required by ORBGRAND and asserted in the driver.
%
% INPUTS:
%   r     : received signal vector (real, any size; can be signed)
%   gamma : Cauchy scale parameter (positive scalar)
%
% OUTPUTS:
%   y_soft : signed approximate LLR, same size as r
%   t1     : location of reliability peak
%   yc     : tail crossover breakpoint
%   coeffs : [c1; c2; c3] hinge coefficients
%
% PERFORMANCE:
%   Uses a persistent cache so breakpoints and coefficients are computed
%   only once per unique gamma value across all function calls.
% =========================================================================

    if gamma <= 0
        error('cauchy_llr_approx: gamma must be positive.');
    end

    r_original_size = size(r);
    r = r(:).';   % work with a row vector internally

    % ------------------------------------------------------------------
    % Persistent cache: breakpoints for each gamma are expensive to
    % compute (requires a fine grid search + fzero).  Cache them so
    % that repeated calls with the same gamma reuse prior results.
    % ------------------------------------------------------------------
    persistent cache_gamma cache_t1 cache_yc cache_coeffs

    cache_hit = false;
    if ~isempty(cache_gamma)
        idx = find(abs(cache_gamma - gamma) < 1e-14, 1);
        if ~isempty(idx)
            t1     = cache_t1(idx);
            yc     = cache_yc(idx);
            coeffs = cache_coeffs(:, idx);
            cache_hit = true;
        end
    end

    if ~cache_hit

        % ---- Derivative of R(y) = |L(y)| for y > 0 ----
        % R'(y) > 0 on (0, t1) and R'(y) < 0 on (t1, inf).
        Rprime = @(y) 2*(y+1)./(gamma^2+(y+1).^2) - 2*(y-1)./(gamma^2+(y-1).^2);

        % ---- Breakpoint 1: exact peak location ----
        % t1 satisfies R'(t1) = 0, solved analytically.
        t1 = sqrt(1 + gamma^2);

        % ---- Breakpoint 2: tail crossover yc ----
        % Find the LAST y > t1 where the descent rate equals 10% of
        % the initial rising slope R'(0+) = 4/(1 + gamma^2).
        Rprime0 = 4.0 / (1.0 + gamma^2);   % R'(y) at y = 0+ (initial slope)
        epsilon  = 0.1;                      % 10% threshold
        target   = epsilon * Rprime0;        % positive target value

        % We search for sign changes of  h(y) = -R'(y) - target  past t1.
        % h > 0 when |R'(y)| > target.  Two sign changes exist:
        %   first  (cross_idx(1)):   on the rising flank of |R'|, near t1
        %   last   (cross_idx(end)): on the falling flank — use this one.
        y_lo     = t1 * (1 + 1e-6);
        y_hi     = max(1000*t1, 1000/gamma);
        y_search = logspace(log10(y_lo), log10(y_hi), 200000);
        h_grid   = -Rprime(y_search) - target;

        cross_idx = find(h_grid(1:end-1) .* h_grid(2:end) <= 0);

        if isempty(cross_idx)
            warning('cauchy_llr_approx: yc crossing not found for gamma=%.4g; using 10*t1.', gamma);
            yc = 10 * t1;
        else
            % Refine the LAST crossing to machine precision with fzero.
            ic = cross_idx(end);
            yc = fzero(@(y) -Rprime(y) - target, [y_search(ic), y_search(ic+1)]);
        end

        % ---- Evaluate R at the two breakpoints ----
        R_exact = @(y) abs(log((gamma^2+(y+1).^2) ./ (gamma^2+(y-1).^2)));
        R_t1    = R_exact(t1);   % peak reliability value
        R_yc    = R_exact(yc);   % reliability at tail crossover

        % ---- Compute hinge coefficients ----
        % c1: slope of rising Segment 1, passing through (0,0) and (t1, R_t1).
        c1 = R_t1 / t1;

        % c2: hinge at t1.  Slope of Segment 2 = c1+c2 < 0.
        % Interpolates the line from (t1, R_t1) to (yc, R_yc).
        slope2 = (R_yc - R_t1) / (yc - t1);   % always negative
        c2     = slope2 - c1;                  % negative hinge correction

        % c3: hinge at yc.  Matches the actual derivative at mid-tail.
        % Constrained so |slope3| <= |slope2|: the tail cannot steepen.
        slope3_raw = Rprime(1.5 * yc);   % negative; should be small magnitude

        if abs(slope3_raw) > abs(slope2)
            % Unphysical steepening in the tail — degenerate case.
            % Extend Segment 2's slope into the tail (c3 = 0).
            c3 = 0;
        else
            c3 = slope3_raw - slope2;   % positive: reduces descent rate at yc
        end

        coeffs = [c1; c2; c3];

        % Store in cache for future calls with the same gamma.
        cache_gamma(end+1)    = gamma;
        cache_t1(end+1)       = t1;
        cache_yc(end+1)       = yc;
        cache_coeffs(:,end+1) = coeffs;

    end % ~cache_hit

    % ------------------------------------------------------------------
    % Evaluate the piecewise-linear approximation on the input vector r.
    % Exploit symmetry: R(r) = R(-r), so work with u = |r|.
    % ------------------------------------------------------------------
    u = abs(r);

    % Hinge-function evaluation: R_approx(u) = c1*u + c2*(u-t1)+ + c3*(u-yc)+
    R_approx = coeffs(1)*u ...
             + coeffs(2)*max(0, u - t1) ...
             + coeffs(3)*max(0, u - yc);

    % Clamp to non-negative: the linear tail may undershoot zero for
    % very large |r| where the true R(y) -> 0 asymptotically.
    R_approx = max(R_approx, 0);

    % Preserve a tiny positive value for all u > 0 so that
    % sign(y_soft) = sign(r), as required by the driver's assertion.
    tiny = 1e-12;
    R_approx(u > 0) = max(R_approx(u > 0), tiny);

    % Reconstruct signed LLR: y_soft has the same sign as r.
    y_soft         = sign(r) .* R_approx;
    y_soft(r == 0) = 0;   % boundary case: decision boundary is LLR = 0

    y_soft = reshape(y_soft, r_original_size);

end
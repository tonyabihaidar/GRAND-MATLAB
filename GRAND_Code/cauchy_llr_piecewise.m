function [R_approx, t1, t2, info] = cauchy_llr_piecewise(y, gamma)
%CAUCHY_LLR_PIECEWISE  3-piece approximation of Cauchy |LLR|
%
%   [R_approx, t1, t2, info] = cauchy_llr_piecewise(y, gamma)
%
% Inputs:
%   y      : scalar/vector of channel outputs
%   gamma  : Cauchy scale parameter, gamma > 0
%
% Outputs:
%   R_approx : piecewise-linear approximation of |LLR|
%   t1       : first threshold = sqrt(1 + gamma^2)
%   t2       : second threshold solving R'(t2) = epsilon * R'(0)
%   info     : struct with exact function/derivative values and line params
%
% Approximation:
%   piece 1 on [0, t1]
%   piece 2 on (t1, t2]
%   piece 3 on (t2, inf)
%
% The approximation is continuous at t1 and t2.
%
% Notes:
%   - Uses only |y| because |LLR| is an even function.
%   - epsilon is fixed to 0.1 as requested.

    if nargin < 2
        error('Usage: [R_approx, t1, t2, info] = cauchy_llr_piecewise(y, gamma)');
    end
    if gamma <= 0
        error('gamma must be positive.');
    end

    epsilon = 0.1;

    % Work with nonnegative values because |LLR| is even
    ya = abs(y);

    % Exact reliability function for y >= 0
    R = @(x) log(gamma^2 + (x + 1).^2) - log(gamma^2 + (x - 1).^2);

    % Derivative of R(x) for x >= 0
    Rp = @(x) 2*(x + 1)./(gamma^2 + (x + 1).^2) ...
            - 2*(x - 1)./(gamma^2 + (x - 1).^2);

    % First threshold
    t1 = sqrt(1 + gamma^2);

    % Reference slope at zero
    Rp0 = Rp(0);
    target = epsilon * Rp0;

    % Solve for t2 > t1 such that Rp(t2) = epsilon*Rp(0)
    % Since Rp(t1)=0 and Rp(x) becomes positive again for large x,
    % we search on (t1, +inf) for the rising branch.
    f = @(x) Rp(x) - target;

    left = t1 + 1e-8;
    right = max(t1 + 1, 5*t1 + 5);

    % Expand right bound until a sign change is found
    max_expand = 50;
    count = 0;
    while f(left) * f(right) > 0 && count < max_expand
        right = 2 * right;
        count = count + 1;
    end

    if f(left) * f(right) > 0
        error('Could not bracket t2. Try checking gamma or increasing search range.');
    end

    t2 = fzero(f, [left, right]);

    % Exact values at thresholds
    R0  = R(0);
    R1  = R(t1);
    R2  = R(t2);
    Rp2 = Rp(t2);

    % Piece 1: line through (0,R(0)) and (t1,R(t1))
    m1 = (R1 - R0) / t1;
    b1 = R0;

    % Piece 2: line through (t1,R(t1)) and (t2,R(t2))
    m2 = (R2 - R1) / (t2 - t1);
    b2 = R1 - m2 * t1;

    % Piece 3: tangent at t2
    m3 = Rp2;
    b3 = R2 - m3 * t2;

    % Evaluate piecewise approximation
    R_approx = zeros(size(ya));

    idx1 = (ya <= t1);
    idx2 = (ya > t1) & (ya <= t2);
    idx3 = (ya > t2);

    R_approx(idx1) = m1 * ya(idx1) + b1;
    R_approx(idx2) = m2 * ya(idx2) + b2;
    R_approx(idx3) = m3 * ya(idx3) + b3;

    % Optional diagnostic info
    info = struct();
    info.epsilon = epsilon;
    info.R       = R;
    info.Rp      = Rp;
    info.R0      = R0;
    info.R1      = R1;
    info.R2      = R2;
    info.Rp0     = Rp0;
    info.Rp2     = Rp2;
    info.target  = target;
    info.m1      = m1;
    info.b1      = b1;
    info.m2      = m2;
    info.b2      = b2;
    info.m3      = m3;
    info.b3      = b3;
end
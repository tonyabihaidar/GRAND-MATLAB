% =========================================================================
% validate_cauchy_approx.m
%
% PURPOSE:
%   Validates the 3-piece piecewise-linear approximation of the Cauchy
%   LLR reliability function implemented in cauchy_llr_approx.m.
%
%   Two aspects are tested:
%
%   1. VISUAL SHAPE VALIDATION (Figures 1 and 2):
%      Plots the exact reliability R(y) = |L(y)| alongside the
%      approximation R_approx(y) for a range of gamma values.
%      Figure 1: full signed shape over [-x_max, +x_max].
%      Figure 2: positive side detail with tail region visible.
%      Kendall's tau rank correlation is reported on each subplot
%      as a global shape-quality score.
%
%   2. MONTE CARLO ORDERING PRESERVATION (printed table):
%      For each gamma, draws N_mc = 10000 noise samples, evaluates
%      exact and approximate reliability on each, then checks what
%      fraction of randomly chosen pairs are ordered consistently.
%      Results are reported per region:
%        Region A [0, t1]:       rising segment — most critical for ORBGRAND
%        Region B [t1, yc]:      descent through the outlier zone
%        Region C [yc, 10*yc]:   far tail where both curves are near zero
%
%   WHY ORDERING MATTERS:
%      ORBGRAND queries noise patterns in order of bit reliability.
%      Incorrect ordering causes the decoder to try unlikely patterns
%      before likely ones, degrading complexity and potentially BLER.
%      The approximation is valid for ORBGRAND if and only if it
%      preserves the rank ordering of |LLR| values across received samples.
%
% DEPENDENCIES:
%   cauchy_llr_approx.m
%
% OUTPUTS:
%   approx_validation_shape.png  — Figure 1 (full shape)
%   approx_validation_detail.png — Figure 2 (positive side, tail visible)
%   Console table with per-region ordering preservation percentages.
% =========================================================================

clear; clc; close all;

% Seed the random number generator for reproducible Monte Carlo results.
rng(42, 'twister');

MIN_PAIRS = 100;   % Minimum valid pairs for a trustworthy ordering estimate.
                   % Regions with fewer pairs are flagged as LOW-N.

% gamma values for plotting (a representative subset).
gamma_plot = [0.05, 0.1, 0.2, 0.3, 0.5];
N_plot     = length(gamma_plot);

% gamma values for the Monte Carlo ordering table (extended range).
gamma_mc = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0];

% =========================================================================
% Figure 1: Full signed shape — |LLR(r)| vs r over [-x_max, +x_max]
%
% Shows both the exact and approximate reliability as a function of
% the raw received value r (signed).  The approximation is symmetric
% in r, so the plot is symmetric about zero.  Vertical markers show
% the two breakpoints t1 (reliability peak) and yc (tail crossover).
% =========================================================================
fig1 = figure('Name','Cauchy |LLR|: Full Shape', 'Position', [50 50 1600 320]);

for gi = 1:N_plot
    gamma = gamma_plot(gi);
    [t1, yc, ~] = get_cauchy_params(gamma);   % retrieve breakpoints for this gamma

    % Set plot range to cover at least 2*yc or 4*t1, capped at 30
    % to avoid uninformative near-zero tails at very large |r|.
    x_max   = min(max(2*yc, 4*t1), 30);
    r_range = linspace(-x_max, x_max, 4000);

    R_exact  = abs(cauchy_llr_exact(r_range, gamma));
    R_approx = abs(cauchy_llr_approx(r_range, gamma));

    % Estimate Kendall's tau rank correlation on random positive-side
    % samples.  tau = 1 means perfect rank agreement; tau = 0 means no
    % rank correlation.  For ORBGRAND, we want tau close to 1.
    y_test    = abs(cauchy_sample(gamma, 2000));
    y_test    = y_test(y_test <= x_max);
    R_ex_test = abs(cauchy_llr_exact(y_test, gamma));
    R_ap_test = abs(cauchy_llr_approx(y_test, gamma));
    tau       = corr(R_ex_test(:), R_ap_test(:), 'type', 'Kendall');

    subplot(1, N_plot, gi);
    plot(r_range, R_exact,  'b-',  'LineWidth', 2.0, 'DisplayName', 'Exact');
    hold on;
    plot(r_range, R_approx, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Approx');

    y_top = max(R_exact) * 1.05;   % top of plot for label placement

    % Mark the reliability peak at +/-t1.
    if t1 <= x_max
        xline( t1, 'k:', 'LineWidth', 1.2);
        xline(-t1, 'k:', 'LineWidth', 1.2);
        text(t1, y_top*0.9, 't_1', 'HorizontalAlignment','center','FontSize',8,'Color','k');
    end

    % Mark the tail crossover at +/-yc.
    if yc <= x_max
        xline( yc, 'm:', 'LineWidth', 1.2);
        xline(-yc, 'm:', 'LineWidth', 1.2);
        text(yc, y_top*0.9, 'y_c', 'HorizontalAlignment','center','FontSize',8,'Color','m');
    else
        % yc is off-screen — note it textually instead.
        text(x_max*0.85, y_top*0.5, sprintf('y_c=%.1f\n(off-screen)',yc), ...
            'FontSize',7,'Color','m','HorizontalAlignment','right');
    end

    xlabel('r'); ylabel('|LLR(r)|');
    title(sprintf('\\gamma=%.2f  \\tau=%.3f', gamma, tau));
    legend('Location','north','FontSize',7);
    xlim([-x_max, x_max]); grid on; box on;
end

sgtitle('Cauchy Reliability |LLR(r)|: Exact vs 3-piece Approx (Full Shape)');
set(fig1,'Color','w');
saveas(fig1,'approx_validation_shape.png');
fprintf('Figure 1 saved: approx_validation_shape.png\n');

% =========================================================================
% Figure 2: Positive-side detail — R(y) vs y = |r| >= 0
%
% Zooms into the positive half-axis to show how well the 3-piece linear
% approximation tracks the exact curve, especially in the tail (y > t1).
% The x-axis covers up to 1.5*yc or 5*t1 so the Segment 2 descent and
% Segment 3 gentle tail are both clearly visible.
% =========================================================================
fig2 = figure('Name','Cauchy |LLR|: Positive Side Detail','Position',[50 420 1600 320]);

for gi = 1:N_plot
    gamma = gamma_plot(gi);
    [t1, yc, ~] = get_cauchy_params(gamma);

    y_max = min(max(1.5*yc, 5*t1), 30);   % ensures tail region is visible
    y_pos = linspace(0, y_max, 3000);

    R_ex_pos = abs(cauchy_llr_exact(y_pos, gamma));
    R_ap_pos = abs(cauchy_llr_approx(y_pos, gamma));
    tau      = corr(R_ex_pos(:), R_ap_pos(:), 'type', 'Kendall');

    subplot(1, N_plot, gi);
    plot(y_pos, R_ex_pos, 'b-',  'LineWidth', 2.0, 'DisplayName', 'Exact');
    hold on;
    plot(y_pos, R_ap_pos, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Approx');

    y_top = max(R_ex_pos) * 1.05;
    xline(t1, 'k:', 'LineWidth', 1.2);
    text(t1, y_top*0.92, 't_1', 'HorizontalAlignment','center','FontSize',8,'Color','k');

    if yc <= y_max
        xline(yc, 'm:', 'LineWidth', 1.2);
        text(yc, y_top*0.92, 'y_c', 'HorizontalAlignment','center','FontSize',8,'Color','m');
    else
        text(y_max*0.85, y_top*0.5, sprintf('y_c=%.1f\n(off-screen)',yc), ...
            'FontSize',7,'Color','m','HorizontalAlignment','right');
    end

    xlabel('y = |r|'); ylabel('R(y) = |LLR(y)|');
    title(sprintf('\\gamma=%.2f  \\tau=%.3f  t_1=%.3f  y_c=%.3f', gamma, tau, t1, yc));
    legend('Location','northeast','FontSize',7);
    xlim([0, y_max]); grid on; box on;
end

sgtitle('Cauchy Reliability (Positive Side): Exact vs Approx — Tail Region Visible');
set(fig2,'Color','w');
saveas(fig2,'approx_validation_detail.png');
fprintf('Figure 2 saved: approx_validation_detail.png\n');

% =========================================================================
% Monte Carlo ordering preservation — stratified by region
%
% For each gamma, draws N_mc Cauchy noise samples, evaluates exact and
% approximate reliability on each, and estimates the fraction of randomly
% chosen pairs (i,j) for which sign(R_exact(i) - R_exact(j)) ==
% sign(R_approx(i) - R_approx(j)).  Pairs where either difference is
% below a tolerance (near-ties) are excluded to avoid noise.
%
% The analysis is stratified into three regions:
%   A [0, t1]:      rising segment — determines which bits ORBGRAND queries first
%   B [t1, yc]:     descent — captures the outlier transition zone
%   C [yc, 10*yc]:  far tail — both curves near zero; ordering matters less
% =========================================================================
fprintf('\n%s\n', repmat('=',1,80));
fprintf('MONTE CARLO ORDERING PRESERVATION  (N=10000, seed=42)\n');
fprintf('%s\n', repmat('=',1,80));
fprintf('%-8s  %-8s  %-8s  %-7s\n','gamma','t1','yc','tau_glob');
fprintf('%s\n', repmat('-',1,80));

N_mc    = 10000;
N_pairs = 10000;
tol     = 1e-9;   % exclude near-tied pairs where ordering is numerically ambiguous

for gi = 1:length(gamma_mc)
    gamma = gamma_mc(gi);
    [t1, yc, ~] = get_cauchy_params(gamma);

    r_mc = cauchy_sample(gamma, N_mc);
    y_mc = abs(r_mc);   % work with positive half only (R is even)

    R_ex = abs(cauchy_llr_exact(y_mc, gamma));
    R_ap = abs(cauchy_llr_approx(y_mc, gamma));

    % Global Kendall's tau across all N_mc samples.
    tau_global = corr(R_ex(:), R_ap(:), 'type', 'Kendall');

    fprintf('%-8.3f  %-8.4f  %-8.4f  %-7.4f\n', gamma, t1, yc, tau_global);

    % Define the three analysis regions.
    regions = {
        [0,   t1  ],  'A [0,   t1]   rising  — most critical';
        [t1,  yc  ],  'B [t1,  yc]   descent — outlier capture';
        [yc, 10*yc],  'C [yc, 10yc]  far tail — approx flat  ';
    };

    for ri = 1:size(regions,1)
        lo    = regions{ri,1}(1);
        hi    = regions{ri,1}(2);
        label = regions{ri,2};

        % Select samples that fall within this region.
        in_r  = y_mc >= lo & y_mc < hi;
        idx_r = find(in_r);
        n_r   = numel(idx_r);

        if n_r < 2
            fprintf('    %s: N/A (n=%d)\n', label, n_r);
            continue;
        end

        % Randomly sample N_pairs pairs from within this region.
        n_p = min(N_pairs, n_r*(n_r-1)/2);
        i1  = randi(n_r, 1, n_p);
        i2  = randi(n_r, 1, n_p);
        g1  = idx_r(i1);
        g2  = idx_r(i2);

        % Compute differences in exact and approximate reliability.
        dex   = R_ex(g1) - R_ex(g2);
        dap   = R_ap(g1) - R_ap(g2);

        % Keep only non-self, non-tied pairs.
        valid = (g1 ~= g2) & (abs(dex) > tol) & (abs(dap) > tol);
        n_val = sum(valid);

        if n_val == 0
            fprintf('    %s: no valid pairs\n', label);
            continue;
        end

        % Fraction of valid pairs where exact and approx agree on ordering.
        pct  = 100 * mean(sign(dex(valid)) == sign(dap(valid)));
        if n_val < MIN_PAIRS
            flag = sprintf('  *** LOW-N (n_valid=%d) ***', n_val);
        else
            flag = sprintf('  (n_valid=%d)', n_val);
        end
        fprintf('    %s: %5.1f%%%s\n', label, pct, flag);
    end
    fprintf('\n');
end

fprintf('%s\n', repmat('-',1,80));
fprintf('Interpretation:\n');
fprintf('  Region A: ordering here determines ORBGRAND query quality most.\n');
fprintf('  Region B: captures reliability decay through the outlier zone.\n');
fprintf('  Region C: both curves near zero — ordering matters less.\n');
fprintf('  LOW-N flag: fewer than %d valid pairs — estimate unreliable.\n', MIN_PAIRS);

% =========================================================================
% Local helper functions
% =========================================================================

function L = cauchy_llr_exact(r, gamma)
% cauchy_llr_exact  Exact signed Cauchy LLR for BPSK (bit 0 -> +1, bit 1 -> -1).
% L > 0 when r > 0 (evidence for bit 0); L < 0 when r < 0 (evidence for bit 1).
    L = log(gamma^2 + (r+1).^2) - log(gamma^2 + (r-1).^2);
end

function r = cauchy_sample(gamma, N)
% cauchy_sample  Draw N i.i.d. samples from Cauchy(0, gamma) using the
% inverse-CDF method: if U ~ Uniform(0,1), then gamma*tan(pi*(U-0.5)) ~ Cauchy(0,gamma).
    r = gamma * tan(pi * (rand(1,N) - 0.5));
end

function [t1, yc, coeffs] = get_cauchy_params(gamma)
% get_cauchy_params  Retrieve breakpoints and coefficients for a given gamma.
% Calls cauchy_llr_approx at a dummy point (r=1) solely to populate the
% persistent cache and return the computed parameters.
    [~, t1, yc, coeffs] = cauchy_llr_approx(1, gamma);
end
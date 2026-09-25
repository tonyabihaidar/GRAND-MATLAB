% =========================================================================
% driver_GRAND_Cauchy_comparison.m
%
% PURPOSE:
%   Monte Carlo simulation comparing three decoders over a BPSK channel
%   with additive Cauchy noise, for a [128, 112] Random Linear Code (RLC).
%
%   The three decoders compared are:
%     1. Hard GRAND   — queries noise patterns by Hamming weight (no soft info)
%     2. ORBGRAND exact   — queries by logistic weight using the exact Cauchy LLR
%     3. ORBGRAND approx  — queries by logistic weight using the 3-piece approx LLR
%
%   This isolates the question: does the approximate LLR in cauchy_llr_approx.m
%   preserve enough reliability ordering information to match exact ORBGRAND?
%
% CHANNEL MODEL:
%   Transmitted symbol: x = 1 - 2*c  (BPSK: bit 0 -> +1, bit 1 -> -1)
%   Received signal:    r = x + z,    z ~ Cauchy(0, gamma)
%   Cauchy noise is sampled as  z = gamma * tan(pi*(U - 0.5)), U ~ Uniform(0,1).
%
%   The Cauchy scale parameter gamma controls noise severity:
%   small gamma = good channel, large gamma = heavy-tailed noise.
%
% METRICS COLLECTED PER DECODER PER gamma:
%   BLER  : Block Error Rate = fraction of decoded blocks with at least one error
%   BER   : Bit Error Rate   = fraction of individual bits decoded incorrectly
%   EG    : Average number of codebook queries per block (complexity metric)
%   AbanRate : Fraction of blocks where decoding was abandoned (max_query reached)
%   DemodRate: Fraction of blocks where the hard decision differed from the true codeword
%              (channel-level error rate, before any decoding)
%
% STOPPING RULE:
%   For each gamma, simulation continues until each decoder accumulates
%   at least err_thresh = 100 block errors, ensuring statistical reliability.
%
% OUTPUTS:
%   Saves results to ../RESULTS/CAUCHY_comparison_ORBGRAND_RLC_128_112_1.mat
%   Generates a 4-panel figure: BLER, BER, EG, Abandonment Rate vs gamma.
%
% DEPENDENCIES:
%   bin_GRAND.m, bin_ORBGRAND.m, cauchy_llr_approx.m, make_RLC.m
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% Simulation parameters
% -------------------------------------------------------------------------
% gamma_list: Cauchy scale parameters to sweep.
% Smaller gamma = less noise (fewer and smaller outliers).
gamma_list = [0.03, 0.05, 0.07, 0.10, 0.15, 0.20];

% Minimum number of block errors to collect per decoder per gamma point.
err_thresh = 100;

% Maximum syndrome queries per block before declaring abandonment.
% 5e6 >> 2^(n-k) = 2^12 = 4096, so abandonment is extremely rare.
% When EG << 5e6, all "errors" are undetected errors (decoder finds a
% wrong codeword at low Hamming weight), not abandoned blocks.
max_query = 5e6;

% -------------------------------------------------------------------------
% Code setup: [n=128, k=112] Random Linear Code over GF(2)
% -------------------------------------------------------------------------
code_class = 'RLC';
n = 128;    % codeword length
k = 112;    % information bits
modulation = 'BPSK';
nmodbits = 1;

if ~strcmpi(modulation, 'BPSK')
    error('This Cauchy implementation is currently valid only for BPSK.');
end

if ~exist('../RESULTS', 'dir')
    mkdir('../RESULTS');
end

% Load or generate the shared RLC.  Using a single fixed code for all
% decoders and all gamma values ensures a fair comparison.
rlc_file = ['../RESULTS/RLC_shared_' num2str(n) '_' num2str(k) '.mat'];

if exist(rlc_file, 'file') == 2
    load(rlc_file, 'G', 'H');
    disp('Loaded existing shared RLC.');
else
    [G, H] = make_RLC(k, n, 0.5);
    save(rlc_file, 'G', 'H');
    disp('Generated and saved new shared RLC.');
end

% Sanity checks on code dimensions.
if size(G, 1) ~= k || size(G, 2) ~= n
    error('G has wrong dimensions.');
end
if size(H, 2) ~= n
    error('H has wrong number of columns.');
end

% -------------------------------------------------------------------------
% Results storage: one struct entry per decoder
% -------------------------------------------------------------------------
decoders = {'GRAND', 'ORBGRAND_exact', 'ORBGRAND_approx'};
N_dec    = length(decoders);
N_gamma  = length(gamma_list);

results = struct();
for d = 1:N_dec
    results(d).name          = decoders{d};
    results(d).num_decoded   = zeros(1, N_gamma);   % total blocks processed
    results(d).num_errs      = zeros(1, N_gamma);   % block errors (incl. abandoned)
    results(d).num_bit_errs  = zeros(1, N_gamma);   % bit errors
    results(d).num_queries   = zeros(1, N_gamma);   % total syndrome checks
    results(d).num_aband     = zeros(1, N_gamma);   % abandoned blocks
    results(d).num_demod_err = zeros(1, N_gamma);   % hard-decision errors (channel)
end

% =========================================================================
% Main simulation loop: sweep over Cauchy scale parameter gamma
% =========================================================================
for ii = 1:N_gamma

    gamma = gamma_list(ii);
    fprintf('\n===== gamma = %.4f =====\n', gamma);

    % Per-gamma accumulators for each decoder (indexed 1..N_dec).
    nd     = zeros(1, N_dec);   % blocks decoded
    ne     = zeros(1, N_dec);   % block errors
    nbe    = zeros(1, N_dec);   % bit errors
    nq     = zeros(1, N_dec);   % query counts
    nabd   = zeros(1, N_dec);   % abandonments
    ndemod = zeros(1, N_dec);   % demodulation errors

    % Run until ALL decoders have accumulated err_thresh block errors.
    while any(ne < err_thresh)

        % -----------------------------------------------------------------
        % Generate one transmitted block and its channel output.
        %
        % BPSK convention: c in {0,1}^n, x = 1 - 2*c
        %   bit 0 -> symbol +1, bit 1 -> symbol -1
        % The same block is passed to all three decoders for a fair comparison.
        % -----------------------------------------------------------------
        u = rand(1, k) > 0.5;          % random information bits
        c = mod(u * G, 2);             % encoded codeword

        if any(mod(H * c.', 2) ~= 0)
            error('Generated codeword fails parity check.');
        end

        x = 1 - 2*c;                                          % BPSK symbols
        z = gamma * tan(pi * (rand(1, n) - 0.5));             % Cauchy noise samples
        r = x + z;                                            % received signal

        % Hard decision: r < 0 maps to bit 1 (negative half-plane),
        % r > 0 maps to bit 0.  This is the GRAND input and also the
        % implicit hard decision embedded in the ORBGRAND soft inputs.
        y_hard = double(r < 0);

        % Exact Cauchy LLR for ORBGRAND-exact.
        % L(r) = log f(r|x=+1) / f(r|x=-1)
        %       = log(gamma^2 + (r+1)^2) - log(gamma^2 + (r-1)^2)
        % Positive when r > 0 (bit 0 more likely), negative when r < 0 (bit 1).
        y_soft_exact = log(gamma^2 + (r+1).^2) - log(gamma^2 + (r-1).^2);

        % Approximate Cauchy LLR for ORBGRAND-approx.
        % Returns sign(r) * R_approx(|r|); see cauchy_llr_approx.m.
        [y_soft_approx, t1, yc] = cauchy_llr_approx(r, gamma);

        % Verify that both soft inputs agree with the hard decision on sign.
        % sign(LLR) determines the hard bit; this must be consistent.
        assert(isequal(double(y_soft_exact < 0), y_hard), ...
            'Exact LLR hard decision mismatch.');
        assert(isequal(double(y_soft_approx < 0), y_hard), ...
            'Approx LLR hard decision mismatch.');

        % -----------------------------------------------------------------
        % Channel-level error indicator: did the hard decision differ from
        % the true codeword on at least one bit?  This is the demodulation
        % error rate — it is a property of the channel, not the decoder.
        % Comparing demod errors vs block errors reveals whether residual
        % errors come from the channel or from decoder mis-identification.
        % -----------------------------------------------------------------
        demod_err = ~isequal(y_hard, c);

        % -----------------------------------------------------------------
        % Decoder 1: Hard GRAND (Hamming-weight ordered, no soft info)
        % -----------------------------------------------------------------
        if ne(1) < err_thresh
            nd(1)  = nd(1) + 1;
            if demod_err, ndemod(1) = ndemod(1) + 1; end
            [y_dec, ~, q, abandoned] = bin_GRAND(H, max_query, y_hard);
            nq(1)  = nq(1) + q;
            [ne, nbe, nabd] = update_counts(1, y_dec, c, y_hard, abandoned, ne, nbe, nabd);
        end

        % -----------------------------------------------------------------
        % Decoder 2: ORBGRAND with exact Cauchy LLR reliability scores
        % -----------------------------------------------------------------
        if ne(2) < err_thresh
            nd(2)  = nd(2) + 1;
            if demod_err, ndemod(2) = ndemod(2) + 1; end
            [y_dec, ~, q, abandoned] = bin_ORBGRAND(H, max_query, y_soft_exact);
            nq(2)  = nq(2) + q;
            [ne, nbe, nabd] = update_counts(2, y_dec, c, y_hard, abandoned, ne, nbe, nabd);
        end

        % -----------------------------------------------------------------
        % Decoder 3: ORBGRAND with 3-piece approximate Cauchy LLR
        % -----------------------------------------------------------------
        if ne(3) < err_thresh
            nd(3)  = nd(3) + 1;
            if demod_err, ndemod(3) = ndemod(3) + 1; end
            [y_dec, ~, q, abandoned] = bin_ORBGRAND(H, max_query, y_soft_approx);
            nq(3)  = nq(3) + q;
            [ne, nbe, nabd] = update_counts(3, y_dec, c, y_hard, abandoned, ne, nbe, nabd);
        end

    end % while (error collection loop)

    % Store accumulated counters into results struct.
    for d = 1:N_dec
        results(d).num_decoded(ii)   = nd(d);
        results(d).num_errs(ii)      = ne(d);
        results(d).num_bit_errs(ii)  = nbe(d);
        results(d).num_queries(ii)   = nq(d);
        results(d).num_aband(ii)     = nabd(d);
        results(d).num_demod_err(ii) = ndemod(d);
    end

    % Print diagnostic breakdown for this gamma value.
    fprintf('\n  gamma=%.4f  Diagnostics:\n', gamma);
    for d = 1:N_dec
        if nd(d) > 0
            fprintf('    [%s]\n', decoders{d});
            fprintf('      Blocks:        %d\n',    nd(d));
            fprintf('      Block errors:  %d  (BLER=%.4g)\n', ne(d), ne(d)/nd(d));
            fprintf('      Abandonments:  %d  (rate=%.4g)\n', nabd(d), nabd(d)/nd(d));
            fprintf('      Demod errors:  %d  (rate=%.4g)\n', ndemod(d), ndemod(d)/nd(d));
            fprintf('      Avg queries:   %.1f\n',  nq(d)/nd(d));
            fprintf('      Undetected:    %d  (BLER - Abandon = %.4g)\n', ...
                ne(d) - nabd(d), (ne(d) - nabd(d))/nd(d));
        end
    end

end % gamma loop

% =========================================================================
% Compute summary statistics
% =========================================================================
for d = 1:N_dec
    results(d).BLER      = results(d).num_errs     ./ results(d).num_decoded;
    results(d).BER       = results(d).num_bit_errs ./ (results(d).num_decoded * n);
    results(d).EG        = results(d).num_queries  ./ results(d).num_decoded;
    results(d).AbanRate  = results(d).num_aband    ./ results(d).num_decoded;
    results(d).DemodRate = results(d).num_demod_err ./ results(d).num_decoded;
end

% Save all results to disk for later plotting with plot_cauchy_results.m.
out_filename = ['../RESULTS/CAUCHY_comparison_ORBGRAND_' ...
    code_class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
save(out_filename, 'results', 'gamma_list', 'n', 'k', 'G', 'H', ...
    'code_class', 'modulation', 'err_thresh', 'max_query');
fprintf('\nResults saved to: %s\n', out_filename);

% Generate inline comparison figure.
plot_cauchy_comparison(results, gamma_list, n, k);


% =========================================================================
% Local helper functions
% =========================================================================

function [ne, nbe, nabd] = update_counts(d, y_dec, c, y_hard, abandoned, ne, nbe, nabd)
% update_counts  Accumulate error and abandonment counters for decoder d.
%
% If the block was abandoned (max_query reached):
%   - Count it as a block error (decoding failed).
%   - Use y_hard for bit-error counting (closest available estimate).
% If the block was decoded but incorrectly:
%   - Count it as a block error with bit errors relative to true codeword.
% If the block was decoded correctly: no counters are incremented.

    if abandoned
        ne(d)   = ne(d) + 1;
        nabd(d) = nabd(d) + 1;
        nbe(d)  = nbe(d) + sum(y_hard ~= c);   % hard bits vs truth
    elseif ~isequal(y_dec, c)
        ne(d)  = ne(d) + 1;
        nbe(d) = nbe(d) + sum(y_dec ~= c);     % decoded bits vs truth
    end
end


function plot_cauchy_comparison(results, gamma_list, n, k)
% plot_cauchy_comparison  Generate and save a 4-panel comparison figure.
%
% Panels:
%   1. BLER vs gamma  (log scale): block-level error rate for each decoder
%   2. BER  vs gamma  (log scale): bit-level error rate
%   3. EG   vs gamma  (linear):    average codebook queries (complexity)
%   4. Abandonment rate vs gamma:  fraction of blocks where decoder gave up
%
% The x-axis is reversed so that low-noise (small gamma, good channel)
% is on the right — consistent with the convention that performance curves
% improve toward the right.

    floor_val = 1e-6;   % floor for semilogy to avoid log(0)
    colors    = {'b', 'r', 'g'};
    markers   = {'o', 's', '^'};
    labels    = {'GRAND hard', ...
                 'ORBGRAND exact Cauchy LLR', ...
                 'ORBGRAND approx Cauchy LLR'};

    figure('Name', 'Cauchy GRAND Comparison', 'Position', [100 100 1800 420]);

    % --- Panel 1: Block Error Rate ---
    subplot(1,4,1); hold on; grid on; box on;
    for d = 1:3
        semilogy(gamma_list, max(results(d).BLER, floor_val), ...
            [colors{d} '-' markers{d}], 'LineWidth', 1.5, 'MarkerSize', 7, ...
            'DisplayName', labels{d});
    end
    xlabel('\gamma'); ylabel('BLER');
    title(sprintf('BLER [n=%d, k=%d]', n, k));
    legend('Location', 'northwest');
    set(gca, 'XDir', 'reverse');

    % --- Panel 2: Bit Error Rate ---
    subplot(1,4,2); hold on; grid on; box on;
    for d = 1:3
        semilogy(gamma_list, max(results(d).BER, floor_val), ...
            [colors{d} '-' markers{d}], 'LineWidth', 1.5, 'MarkerSize', 7, ...
            'DisplayName', labels{d});
    end
    xlabel('\gamma'); ylabel('BER');
    title(sprintf('BER [n=%d, k=%d]', n, k));
    legend('Location', 'northwest');
    set(gca, 'XDir', 'reverse');

    % --- Panel 3: Average queries (decoding complexity) ---
    subplot(1,4,3); hold on; grid on; box on;
    for d = 1:3
        plot(gamma_list, results(d).EG, [colors{d} '-' markers{d}], ...
            'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', labels{d});
    end
    xlabel('\gamma'); ylabel('Avg codebook queries');
    title(sprintf('Complexity [n=%d, k=%d]', n, k));
    legend('Location', 'northwest');
    set(gca, 'XDir', 'reverse');

    % --- Panel 4: Abandonment rate ---
    % Abandonment occurs when max_query is reached before finding a codeword.
    % If this is near zero, EG values are not inflated by abandoned blocks.
    subplot(1,4,4); hold on; grid on; box on;
    for d = 1:3
        semilogy(gamma_list, max(results(d).AbanRate, floor_val), ...
            [colors{d} '-' markers{d}], 'LineWidth', 1.5, 'MarkerSize', 7, ...
            'DisplayName', labels{d});
    end
    xlabel('\gamma'); ylabel('Abandonment rate');
    title('Abandonment Rate');
    legend('Location', 'northwest');
    set(gca, 'XDir', 'reverse');

    sgtitle('GRAND vs ORBGRAND under Cauchy Noise');
    saveas(gcf, '../RESULTS/cauchy_comparison_plot.png');
    disp('Figure saved.');
end
% =========================================================================
% check_cauchy_ranking.m
%
% Checks whether exact Cauchy LLR and 3-piece approximate LLR give similar
% bit reliability rankings.
%
% Requires:
%   cauchy_llr_approx.m
% =========================================================================

clear; clc; close all;

gamma_values = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0];

n = 128;        % blocklength
Nblocks = 500; % number of random received blocks per gamma

fprintf('Ranking comparison: exact LLR vs approximate LLR\n\n');

for gi = 1:length(gamma_values)

    gamma = gamma_values(gi);

    total_same_top10 = 0;
    total_same_top20 = 0;
    total_spearman = 0;
    total_kendall = 0;

    example_saved = false;

    for b = 1:Nblocks

        % Random codeword-like bits just for channel testing
        c = rand(1,n) > 0.5;

        % BPSK: 0 -> +1, 1 -> -1
        x = 1 - 2*c;

        % Cauchy noise
        z = gamma * tan(pi * (rand(1,n) - 0.5));
        r = x + z;

        % Exact Cauchy LLR
        L_exact = log(gamma^2 + (r+1).^2) ...
                - log(gamma^2 + (r-1).^2);

        % Approximate Cauchy LLR
        L_approx = cauchy_llr_approx(r, gamma);

        % Reliability = absolute LLR
        R_exact = abs(L_exact);
        R_approx = abs(L_approx);

        % ORBGRAND cares about least reliable bits first
        [~, order_exact] = sort(R_exact, 'ascend');
        [~, order_approx] = sort(R_approx, 'ascend');

        % Compare top unreliable bits
        top10_exact = order_exact(1:10);
        top10_approx = order_approx(1:10);

        top20_exact = order_exact(1:20);
        top20_approx = order_approx(1:20);

        total_same_top10 = total_same_top10 + ...
            numel(intersect(top10_exact, top10_approx)) / 10;

        total_same_top20 = total_same_top20 + ...
            numel(intersect(top20_exact, top20_approx)) / 20;

        % Rank position vectors
        rank_exact = zeros(1,n);
        rank_approx = zeros(1,n);

        rank_exact(order_exact) = 1:n;
        rank_approx(order_approx) = 1:n;

        % Spearman correlation
        total_spearman = total_spearman + corr(rank_exact.', rank_approx.', ...
            'Type', 'Spearman');

        % Kendall correlation
        total_kendall = total_kendall + corr(rank_exact.', rank_approx.', ...
            'Type', 'Kendall');

        % Save one example for plotting
        if ~example_saved
            example_R_exact = R_exact;
            example_R_approx = R_approx;
            example_order_exact = order_exact;
            example_order_approx = order_approx;
            example_rank_exact = rank_exact;
            example_rank_approx = rank_approx;
            example_saved = true;
        end

    end

    avg_same_top10 = 100 * total_same_top10 / Nblocks;
    avg_same_top20 = 100 * total_same_top20 / Nblocks;
    avg_spearman = total_spearman / Nblocks;
    avg_kendall = total_kendall / Nblocks;

    fprintf('gamma = %.3f\n', gamma);
    fprintf('  Same top-10 unreliable bits = %.2f%%\n', avg_same_top10);
    fprintf('  Same top-20 unreliable bits = %.2f%%\n', avg_same_top20);
    fprintf('  Avg Spearman rank corr      = %.4f\n', avg_spearman);
    fprintf('  Avg Kendall rank corr       = %.4f\n\n', avg_kendall);

    % Plot one example per gamma
    figure('Name', sprintf('Ranking check gamma %.3f', gamma), ...
           'Position', [100 100 1300 420]);

    subplot(1,3,1);
    plot(example_R_exact, 'b-o', 'LineWidth', 1.2);
    hold on;
    plot(example_R_approx, 'r--s', 'LineWidth', 1.2);
    xlabel('Bit index');
    ylabel('|LLR|');
    title(sprintf('Reliability values, \\gamma=%.3f', gamma));
    legend('Exact', 'Approx');
    grid on; box on;

    subplot(1,3,2);
    plot(example_order_exact, 'b-o', 'LineWidth', 1.2);
    hold on;
    plot(example_order_approx, 'r--s', 'LineWidth', 1.2);
    xlabel('Reliability rank');
    ylabel('Bit index');
    title('Least reliable to most reliable');
    legend('Exact order', 'Approx order');
    grid on; box on;

    subplot(1,3,3);
    scatter(example_rank_exact, example_rank_approx, 45, 'filled');
    xlabel('Exact rank position');
    ylabel('Approx rank position');
    title('Rank agreement');
    grid on; box on;
    axis square;

    sgtitle(sprintf('Exact vs Approx Reliability Ranking, \\gamma=%.3f', gamma));

end
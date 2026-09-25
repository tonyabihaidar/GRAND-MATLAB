% =========================================================================
% plot_cauchy_results.m
%
% Load saved comparison results and generate publication-quality figures.
% Run after driver_GRAND_Cauchy_comparison.m has finished.
% =========================================================================

clear; clc;

% Load results
filename = '../RESULTS/CAUCHY_comparison_RLC_128_116_1.mat';
if ~exist(filename, 'file')
    error('Results file not found. Run driver_GRAND_Cauchy_comparison.m first.');
end
load(filename);

colors  = {[0.1 0.1 0.8], [0.8 0.1 0.1], [0.1 0.6 0.1]};
markers = {'o', 's', '^'};
lw = 2; ms = 8;

labels = {'GRAND (hard decision)', ...
          'ORBGRAND – exact Cauchy LLR', ...
          'ORBGRAND – 3-piece approx LLR'};

% -------------------------------------------------------------------------
% Figure 1: BLER and BER
% -------------------------------------------------------------------------
fig1 = figure('Name','BLER and BER vs gamma', 'Position',[100 100 900 400]);

subplot(1,2,1);
hold on; grid on; box on;
for d = 1:3
    semilogy(gamma_list, results(d).BLER, ...
        'Color', colors{d}, 'Marker', markers{d}, ...
        'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', labels{d});
end
xlabel('\gamma (Cauchy scale parameter)', 'FontSize', 12);
ylabel('Block Error Rate (BLER)', 'FontSize', 12);
title('BLER vs \gamma', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
set(gca, 'XDir', 'reverse');   % small gamma = good channel on right
ylim([1e-3 1]); xlim([gamma_list(1) gamma_list(end)]);

subplot(1,2,2);
hold on; grid on; box on;
for d = 1:3
    semilogy(gamma_list, results(d).BER, ...
        'Color', colors{d}, 'Marker', markers{d}, ...
        'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', labels{d});
end
xlabel('\gamma (Cauchy scale parameter)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER vs \gamma', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
set(gca, 'XDir', 'reverse');
ylim([1e-5 1]); xlim([gamma_list(1) gamma_list(end)]);

sgtitle(sprintf('Cauchy Noise Decoding  [n=%d, k=%d, R=%.3f]', n, k, k/n), ...
    'FontSize', 13, 'FontWeight', 'bold');
set(fig1, 'Color', 'w');
saveas(fig1, '../RESULTS/fig_BLER_BER.png');
saveas(fig1, '../RESULTS/fig_BLER_BER.pdf');

% -------------------------------------------------------------------------
% Figure 2: Average queries (complexity)
% -------------------------------------------------------------------------
fig2 = figure('Name','Complexity vs gamma', 'Position',[100 550 520 380]);
hold on; grid on; box on;
for d = 1:3
    plot(gamma_list, results(d).EG, ...
        'Color', colors{d}, 'Marker', markers{d}, ...
        'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', labels{d});
end
xlabel('\gamma (Cauchy scale parameter)', 'FontSize', 12);
ylabel('Avg. codebook queries per packet', 'FontSize', 12);
title('Decoding Complexity vs \gamma', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
set(gca, 'XDir', 'reverse');
set(fig2, 'Color', 'w');
saveas(fig2, '../RESULTS/fig_complexity.png');
saveas(fig2, '../RESULTS/fig_complexity.pdf');

% -------------------------------------------------------------------------
% Figure 3: BLER gap between exact and approx ORBGRAND
% -------------------------------------------------------------------------
fig3 = figure('Name','Approx vs Exact gap', 'Position',[640 550 520 380]);
hold on; grid on; box on;
bler_gap = results(3).BLER ./ results(2).BLER;
plot(gamma_list, bler_gap, 'k-o', 'LineWidth', lw, 'MarkerSize', ms);
yline(1, 'b--', 'LineWidth', 1.2, 'Label', 'Exact');
xlabel('\gamma (Cauchy scale parameter)', 'FontSize', 12);
ylabel('BLER_{approx} / BLER_{exact}', 'FontSize', 12);
title('Approx vs Exact ORBGRAND: BLER Ratio', 'FontSize', 13);
set(gca, 'XDir', 'reverse');
ylim([0.5 3]);
set(fig3, 'Color', 'w');
saveas(fig3, '../RESULTS/fig_approx_gap.png');
saveas(fig3, '../RESULTS/fig_approx_gap.pdf');

disp('All figures saved to ../RESULTS/');

% -------------------------------------------------------------------------
% Print summary table
% -------------------------------------------------------------------------
fprintf('\n%-25s | %-6s | %-8s | %-8s | %-8s\n', ...
    'Decoder', 'gamma', 'BLER', 'BER', 'Avg Q');
fprintf('%s\n', repmat('-',1,65));

for d = 1:3
    for ii = 1:length(gamma_list)
        fprintf('%-25s | %-6.3f | %-8.4f | %-8.5f | %-8.1f\n', ...
            results(d).name, gamma_list(ii), ...
            results(d).BLER(ii), results(d).BER(ii), results(d).EG(ii));
    end
    fprintf('%s\n', repmat('-',1,65));
end

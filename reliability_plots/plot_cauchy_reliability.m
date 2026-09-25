function plot_cauchy_reliability()
% =========================================================================
% plot_cauchy_reliability.m
%
% PURPOSE:
%   Plots the average bit reliability profile for BPSK over a Cauchy
%   (impulsive) noise channel, across multiple SNR values.
%
% BACKGROUND:
%   The Cauchy distribution has heavier tails than the Gaussian: large
%   noise outliers occur with non-negligible probability.  This produces
%   a distinctive reliability profile — most bits are highly reliable,
%   but a small fraction receive extreme noise and become very unreliable.
%   Consequently, the sorted reliability profile rises steeply for the
%   majority of bits and then flattens (or drops) for the outlier-affected
%   bits at the high-rank end.
%
%   This shape is qualitatively different from AWGN, where reliability
%   rises smoothly across all ranks.  Understanding this difference
%   motivates using a Cauchy-specific LLR in ORBGRAND rather than the
%   AWGN LLR.
%
% CHANNEL MODEL:
%   BPSK: bit x in {0,1} -> symbol s = 1 - 2*x in {+1, -1}.
%   Cauchy noise: z = gamma * tan(pi*(U - 0.5)), U ~ Uniform(0,1).
%   Received: r = s + z.
%
%   Cauchy LLR for BPSK:
%     L(r) = log(1 + ((r-1)/gamma)^2) - log(1 + ((r+1)/gamma)^2)
%   (or equivalently: log(gamma^2+(r+1)^2) - log(gamma^2+(r-1)^2) with sign flipped)
%
%   NOTE: The LLR formula used here has the argument order reversed
%   relative to the standard convention in the driver.  Both are valid;
%   this version is positive when r > 0 (evidence for bit 1) vs the
%   driver convention (positive when r > 0 means bit 0).  The reliability
%   |LLR| is the same in both cases.
%
% SNR CONVENTION:
%   gamma = 10^(-SNR_dB / 20), matching the convention in plot_awgn and
%   plot_lwgn for cross-noise-model comparison.
%
% RELIABILITY PROFILE:
%   For each frame, sort |LLR| in ascending order (least to most reliable).
%   Average across Nframes to get E[|LLR|_(i)] for rank i = 1..n.
% =========================================================================

    n       = 128;
    Nframes = 20000;           % frames to average for smooth curves
    SNRdB_list = [-10 0 4 7]; % SNR values in dB (lower = noisier / larger gamma)

    figure; hold on; grid on;

    for snrdb = SNRdB_list

        gamma = 10^(-snrdb/20);   % Cauchy scale parameter (heavier tail = larger gamma)

        % Generate a random binary matrix (n bits x Nframes frames) and modulate.
        x = randi([0 1], n, Nframes);
        s = 1 - 2*x;   % BPSK symbols: 0 -> +1, 1 -> -1

        % Generate Cauchy noise using the inverse-CDF method.
        % If U ~ Uniform(0,1), then gamma*tan(pi*(U-0.5)) ~ Cauchy(0,gamma).
        U = rand(n, Nframes);
        z = gamma * tan(pi*(U - 0.5));
        r = s + z;   % received signal matrix (n x Nframes)

        % Compute the exact Cauchy LLR for each received sample.
        % The reliability of each received bit is |LLR|.
        LLR = log( 1 + ((r - 1)./gamma).^2 ) ...
            - log( 1 + ((r + 1)./gamma).^2 );

        Gamma = abs(LLR);   % reliability = |LLR|, always non-negative

        % Sort each frame's reliability vector ascending and average across frames.
        % Result: Gamma_mean(i) = average reliability of the i-th least reliable bit.
        Gamma_sorted = sort(Gamma, 1, 'ascend');   % sort within each column
        Gamma_mean   = mean(Gamma_sorted, 2);       % average across frames

        plot(1:n, Gamma_mean, 'LineWidth', 1.6, ...
            'DisplayName', sprintf('SNR=%g dB, \\gamma=%.5g', snrdb, gamma));
    end

    xlabel('Sorted bit position (1 = least reliable, n = most reliable)');
    ylabel('Average reliability E[|LLR|]');
    legend('Location','northwest');
    title('Cauchy reliability vs sorted bit position (n=128)');
end
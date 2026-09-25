function plot_lwgn_reliability()
% =========================================================================
% plot_lwgn_reliability.m
%
% PURPOSE:
%   Plots the average bit reliability profile for BPSK over a Laplacian
%   (AWLN — Additive White Laplacian Noise) channel, across multiple SNR
%   values.
%
% BACKGROUND:
%   The Laplacian distribution is intermediate between Gaussian and Cauchy:
%   it has heavier tails than the Gaussian but lighter tails than the Cauchy.
%   Its reliability profile under BPSK has a different shape than either,
%   because the Laplacian LLR is piecewise linear in r (vs quadratic for
%   AWGN and the non-monotone Cauchy case).
%
%   Comparing this reliability profile with the AWGN and Cauchy counterparts
%   provides intuition for why noise-matched LLRs improve ORBGRAND decoding.
%
% CHANNEL MODEL:
%   BPSK: bit x in {0,1} -> symbol s = 1 - 2*x in {+1, -1}.
%
%   Laplacian noise with scale b: p(z) = (1/2b)*exp(-|z|/b).
%   Inverse-CDF sampling: if U ~ Uniform(-0.5, 0.5),
%     z = -b * sign(U) * log(1 - 2|U|)  ~ Laplacian(0, b).
%
%   Received: r = s + z.
%
%   Laplacian LLR for BPSK:
%     L(r) = (|r-1| - |r+1|) / b
%   This is piecewise linear in r with a breakpoint at r = 0:
%     r > 1:   L = -2/b  (constant, strongly favors bit 0)
%     0 < r < 1: L = -2r/b (linearly increasing magnitude)
%     r < 0:   L = +2/b  (constant, strongly favors bit 1)  [sign-flipped region]
%
% SNR CONVENTION:
%   sigma = 10^(-SNR_dB / 20), b = sigma / sqrt(2).
%   This makes sigma the standard deviation of the Laplacian distribution
%   (consistent with the AWGN sigma convention for cross-model comparison).
%
% RELIABILITY PROFILE:
%   For each frame, sort |LLR| ascending and average across Nframes.
%   The resulting profile E[|LLR|_(i)] shows how the expected reliability
%   varies from the least to the most reliable bit under this noise model.
% =========================================================================

    n       = 128;
    Nframes = 20000;           % frames to average for smooth curves
    SNRdB_list = [-10 0 4 7]; % SNR values in dB (lower = noisier)

    figure; hold on; grid on;

    for snrdb = SNRdB_list

        sigma = 10^(-snrdb/20);   % noise spread parameter
        b     = sigma / sqrt(2);  % Laplacian scale: Var = 2b^2, so sigma^2 = 2b^2

        % Generate random binary bits and BPSK symbols.
        x = randi([0 1], n, Nframes);
        s = 1 - 2*x;   % BPSK: 0 -> +1, 1 -> -1

        % Sample Laplacian noise via the inverse-CDF method.
        % U ~ Uniform(0,1); shift to Uniform(-0.5, 0.5), then apply formula.
        U = rand(n, Nframes) - 0.5;
        z = -b .* sign(U) .* log(1 - 2*abs(U));
        r = s + z;   % received signal (n x Nframes)

        % Compute the Laplacian LLR for each received sample.
        % |r-1| - |r+1| measures how much closer r is to symbol +1 vs -1.
        LLR = (abs(r-1) - abs(r+1)) / b;

        % Sort and average the reliability |LLR| across frames.
        Gamma_sorted = sort(abs(LLR), 1, 'ascend');   % ascending within each frame
        Gamma_mean   = mean(Gamma_sorted, 2);          % average across frames

        plot(1:n, Gamma_mean, 'LineWidth', 1.6, ...
            'DisplayName', sprintf('SNR=%g dB, b=%.5g', snrdb, b));
    end

    xlabel('Sorted bit position (1 = least reliable, n = most reliable)');
    ylabel('Average reliability E[|LLR|]');
    title('Laplacian (AWLN) reliability vs sorted bit position (n=128)');
    legend('show');
end
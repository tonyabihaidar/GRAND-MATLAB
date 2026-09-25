function plot_awgn_reliability()
% =========================================================================
% plot_awgn_reliability.m
%
% PURPOSE:
%   Plots the average bit reliability profile for BPSK over an AWGN
%   (Additive White Gaussian Noise) channel, across multiple SNR values.
%
% BACKGROUND:
%   In ORBGRAND, bits are sorted by their LLR magnitude (reliability)
%   before decoding.  The shape of this sorted reliability profile
%   determines how quickly the decoder finds the correct noise pattern:
%   a steep, well-spread profile means most bits are highly reliable
%   (easy to decode), while a flat profile means many bits have similar
%   reliability (harder to distinguish errors).
%
%   This function averages the sorted reliability profile over many
%   random transmitted frames, producing a smooth expected curve for
%   each SNR level.  This is compared against the Cauchy and Laplacian
%   counterparts in plot_cauchy_reliability.m and plot_lwgn_reliability.m.
%
% CHANNEL MODEL:
%   BPSK: bit x in {0,1} is mapped to symbol s = 1 - 2*x in {+1, -1}.
%   Received: r = s + sigma * randn,   where sigma = 10^(-SNR_dB / 20).
%
%   AWGN LLR for BPSK:
%     L(r) = log P(r | s=+1) / P(r | s=-1) = (2/sigma^2) * r
%   This is linear in r, so reliability |L(r)| = (2/sigma^2) * |r|.
%
% RELIABILITY PROFILE:
%   For each frame, sort |L| in ascending order (least to most reliable).
%   Average across Nframes to estimate E[|L|_(i)] for rank i = 1..n.
%   Bit at rank 1 is on average the least reliable (most likely in error),
%   and is thus the first candidate to flip in ORBGRAND.
%
% PARAMETERS:
%   n         : codeword length (128 bits)
%   Nframes   : number of frames to average over (20000)
%   SNRdB_list: SNR values in dB to plot
%
% NOTE ON sigma CONVENTION:
%   sigma = 10^(-SNR_dB / 20) maps SNR_dB = 0 -> sigma = 1,
%   SNR_dB > 0 -> sigma < 1 (less noise).  This matches the convention
%   used in the companion Cauchy and Laplacian reliability scripts.
% =========================================================================

    n      = 128;
    Nframes = 20000;              % frames to average for smooth curves
    SNRdB_list = [-10 0 4 7];    % SNR values in dB (lower = noisier)

    figure; hold on; grid on;

    for snrdb = SNRdB_list

        sigma  = 10^(-snrdb/20);   % noise standard deviation
        sigma2 = sigma^2;

        % Accumulator for the sorted reliability sum across frames.
        acc = zeros(n,1);

        for t = 1:Nframes

            % Generate a random binary codeword and BPSK modulate.
            x = randi([0 1], n, 1);
            s = 1 - 2*x;            % BPSK: 0 -> +1, 1 -> -1

            % Add AWGN.
            r = s + sigma*randn(n,1);

            % Compute the AWGN LLR and its magnitude (reliability).
            % LLR = (2/sigma^2)*r  for BPSK under AWGN.
            LLR   = (2*r) / sigma2;
            Gamma = abs(LLR);

            % Sort reliabilities ascending: Gamma_sorted(1) = least reliable bit.
            Gamma_sorted = sort(Gamma, 'ascend');

            % Accumulate sorted reliabilities for averaging.
            acc = acc + Gamma_sorted;
        end

        % Average over frames to get the expected sorted reliability profile.
        Gamma_mean = acc / Nframes;

        plot(1:n, Gamma_mean, 'LineWidth', 1.6, ...
            'DisplayName', sprintf('SNR=%ddB, σ=%0.5g', snrdb, sigma));
    end

    xlabel('Sorted bit position (1 = least reliable, n = most reliable)');
    ylabel('Average reliability E[|LLR|]');
    title('AWGN reliability vs sorted bit position (n=128)');
    legend show
end
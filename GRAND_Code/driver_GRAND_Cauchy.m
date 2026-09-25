% =========================================================================
% driver_GRAND_Cauchy.m
%
% Single-decoder simulation under BPSK + Cauchy noise.
%
% Decoder choices:
%   DECODER = 'GRAND'
%   DECODER = 'ORBGRAND'
%
% Reliability choices for ORBGRAND:
%   RELIABILITY = 'exact'
%   RELIABILITY = 'approx'
%
% No ORBGRAND1 is used.
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% Choose experiment
% -------------------------------------------------------------------------
DECODER = 'ORBGRAND';       % 'GRAND' or 'ORBGRAND'
RELIABILITY = 'approx';     % 'exact' or 'approx'; ignored if DECODER='GRAND'

gamma_list = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0];
err_thresh = 100;
max_query = 5e6;

% -------------------------------------------------------------------------
% Code setup
% -------------------------------------------------------------------------
code_class = 'RLC';
n = 128;
k = 116;
modulation = 'BPSK';
nmodbits = 1;

if ~strcmpi(modulation, 'BPSK')
    error('This Cauchy implementation is currently valid only for BPSK.');
end

if ~exist('../RESULTS', 'dir')
    mkdir('../RESULTS');
end

if strcmpi(code_class, 'RLC')

    rlc_file = ['../RESULTS/RLC_shared_' num2str(n) '_' num2str(k) '.mat'];

    if exist(rlc_file, 'file') == 2
        load(rlc_file, 'G', 'H');
        disp('Loaded existing shared RLC.');
    else
        [G, H] = make_RLC(k, n, 0.5);
        save(rlc_file, 'G', 'H');
        disp('Generated and saved new shared RLC.');
    end

else
    error('This simplified driver currently includes only RLC. Add other code classes if needed.');
end

if size(G, 1) ~= k || size(G, 2) ~= n
    error('G has wrong dimensions.');
end

if size(H, 2) ~= n
    error('H has wrong number of columns.');
end

% -------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------
snr_db = -20 * log10(gamma_list);

num_decoded = zeros(1, length(gamma_list));
num_demod_errs = zeros(1, length(gamma_list));
num_demod_bit_errs = zeros(1, length(gamma_list));
num_errs = zeros(1, length(gamma_list));
num_bit_errs = zeros(1, length(gamma_list));
num_aband = zeros(1, length(gamma_list));
num_queries = zeros(1, length(gamma_list));

% -------------------------------------------------------------------------
% Simulation
% -------------------------------------------------------------------------
for ii = 1:length(gamma_list)

    gamma = gamma_list(ii);

    fprintf('\n===== gamma = %.4f =====\n', gamma);

    while num_errs(ii) < err_thresh

        num_decoded(ii) = num_decoded(ii) + 1;

        % Encode
        u = rand(1, k) > 0.5;
        c = mod(u * G, 2);

        if any(mod(H * c.', 2) ~= 0)
            error('Generated codeword does not satisfy H*c''=0.');
        end

        % BPSK: bit 0 -> +1, bit 1 -> -1
        x = 1 - 2*c;

        % Add Cauchy noise
        z = gamma * tan(pi * (rand(1, n) - 0.5));
        r = x + z;

        % Hard decision
        y_demod = r < 0;

        % Exact Cauchy LLR
        y_soft_exact = log(gamma^2 + (r+1).^2) ...
                     - log(gamma^2 + (r-1).^2);

        % Approx Cauchy LLR
        [y_soft_approx, t1, yc] = cauchy_llr_approx(r, gamma);

        % Demodulation error count
        if ~isequal(y_demod, c)
            num_demod_errs(ii) = num_demod_errs(ii) + 1;
            num_demod_bit_errs(ii) = num_demod_bit_errs(ii) + sum(abs(y_demod - c));
        end

        % Decode
        if strcmpi(DECODER, 'GRAND')

            [y_decoded, ~, n_guess, abandoned] = bin_GRAND(H, max_query, y_demod);

        elseif strcmpi(DECODER, 'ORBGRAND')

            if strcmpi(RELIABILITY, 'exact')
                soft_in = y_soft_exact;
            elseif strcmpi(RELIABILITY, 'approx')
                soft_in = y_soft_approx;
            else
                error('RELIABILITY must be exact or approx.');
            end

            if ~isequal((soft_in < 0), y_demod)
                error('Soft-input sign does not match hard demodulation.');
            end

            [y_decoded, ~, n_guess, abandoned] = bin_ORBGRAND(H, max_query, soft_in);

        else
            error('DECODER must be GRAND or ORBGRAND.');
        end

        num_queries(ii) = num_queries(ii) + n_guess;

        % Decoding error count
        if ~isequal(y_decoded, c)

            num_errs(ii) = num_errs(ii) + 1;

            if abandoned
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_demod - c));
                num_aband(ii) = num_aband(ii) + 1;
            else
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_decoded - c));
            end

            fprintf('n=%d k=%d R=%.2f gamma=%.4f decoder=%s reliability=%s Dec=%d Err=%d BLER=%.4g BER=%.4g EG=%.0f t1=%.3f yc=%.3f\n', ...
                n, k, k/n, gamma, DECODER, RELIABILITY, ...
                num_decoded(ii), num_errs(ii), ...
                num_errs(ii)/num_decoded(ii), ...
                num_bit_errs(ii)/(num_decoded(ii)*n), ...
                num_queries(ii)/num_decoded(ii), t1, yc);
        end
    end
end

% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
if strcmpi(DECODER, 'GRAND')
    rel_label = 'hard';
else
    rel_label = RELIABILITY;
end

filename = ['../RESULTS/CAUCHY_' DECODER '_' rel_label '_' ...
    code_class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];

code.class = code_class;
code.n = n;
code.k = k;
code.R = k/n;
code.modulation = modulation;
code.nmodbits = nmodbits;
code.G = G;
code.H = H;
code.gamma = gamma_list;
code.snr = snr_db;
code.num_decoded = num_decoded;
code.num_demod_errs = num_demod_errs;
code.num_demod_bit_errs = num_demod_bit_errs;
code.num_errs = num_errs;
code.num_bit_errs = num_bit_errs;
code.num_aband = num_aband;
code.num_queries = num_queries;
code.BLERdemod = num_demod_errs ./ num_decoded;
code.BERdemod = num_demod_bit_errs ./ (num_decoded * n);
code.BLER = num_errs ./ num_decoded;
code.BER = num_bit_errs ./ (num_decoded * n);
code.EG = num_queries ./ num_decoded;
code.max_query = max_query;
code.decoder = DECODER;
code.reliability = rel_label;

save(filename, 'code');

fprintf('\nDone. Results saved to %s\n', filename);
function [y_decoded,err_vec,n_guesses,abandoned] = bin_ORBGRAND1(H,max_query,y_soft)

    % Force row vector
    y_soft = y_soft(:).';

    % Hard demodulate
    y_demod = (y_soft < 0);
    n = length(y_demod);

    n_guesses = 1;
    err_vec = zeros(1,n);
    abandoned = 0;

    y_decoded = -1*ones(size(y_demod));

    % First query: demodulated word
    Hy = mod(H*y_demod',2);

    if isequal(Hy, zeros(size(Hy)))
        y_decoded = y_demod;
        return;
    end

    % Reliability ordering: least reliable first
    reliability = abs(y_soft);
    [L, ind_order] = sort(reliability,'ascend');

    % Inverse permutation
    inv_perm = zeros(1,n);
    for ii = 1:n
        inv_perm(ind_order(ii)) = ii;
    end

    % Robust slope estimate
    mid = round(n/2);
    beta = (L(mid) - L(1)) / max(mid-1,1);

    % Protect against zero/tiny beta due to ties or approximation clipping
    if beta <= 1e-12 || isnan(beta) || isinf(beta)
        c = 0;
    else
        c = max(round(L(1)/beta - 1), 0);
    end

    % Reordered parity-check matrix
    test_H = H(:,ind_order);

    % Total logistic weight
    wt = c + 1;

    while n_guesses < max_query && wt <= c*n + n*(n+1)/2

        w = max(1, ceil((1 + 2*(n+c) - sqrt((1 + 2*(n+c))^2 - 8*wt))/2));

        while w <= n

            W = wt - c*w;

            if W < w*(w+1)/2
                break;
            end

            noise_locations = landslide(W,w,n);

            for jj = 1:size(noise_locations,1)

                n_guesses = n_guesses + 1;

                err_ordered = zeros(1,n);
                err_ordered(noise_locations(jj,:)) = 1;

                if isequal(Hy, mod(test_H*err_ordered',2))

                    err_vec = err_ordered(inv_perm);
                    y_decoded = mod(y_demod - err_vec, 2);
                    abandoned = 0;
                    return;

                end
            end

            w = w + 1;
        end

        wt = wt + 1;
    end

    abandoned = 1;
    err_vec = zeros(size(y_demod));

end
function [rate,denominator] = compute_rate(A, v, p, alpha, sigma2_dB)
    sigma2 = 10^(sigma2_dB/10);
    [num_user, num_irs_elements, ~] = size(A);
    cond = abs(sum(abs(v))-num_irs_elements)<1e-3;
    assert(cond,'phase shifts v not valid')
    rate = 0;
    for kk=1:num_user
        a_kk =A(kk,:,kk).';
        p_k = p(kk);
        numerator = abs(a_kk.'*v)^2*p_k;
        denominator = sigma2;
        for jj=1:num_user
            if jj ~=kk
                a_kj = A(jj,:,kk).';
                p_j = p(jj);
                snr_jk = abs(a_kj.'*v)^2*p_j;
                denominator = denominator+snr_jk;
            end
        end
        sinr_k = numerator/denominator;
        rate_k = alpha(kk)*log2(1+sinr_k);
        rate = rate+rate_k;
        denominator = denominator-sigma2;
    end
end
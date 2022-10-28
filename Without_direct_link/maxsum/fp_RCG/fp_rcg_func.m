function v = fp_rcg_func(A,p,alpha,sigma2_dB,v0)
[num_user, num_elements_irs, ~] = size(A);
sigma2 = 10^(sigma2_dB/10);
y = zeros(num_user,1);
gamma = zeros(num_user,1);
max_iter = 1000;
obj = nan(max_iter,1);
v = v0;
for ii = 1:max_iter
    %% update gamma
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
        gamma(kk) = numerator/denominator;
    end
    
    %% update y
    for kk=1:num_user
        a_kk =A(kk,:,kk).';
        p_k = p(kk);
        numerator = sqrt(p_k*alpha(kk)*(1+gamma(kk)))*(v'*conj(a_kk));
        denominator = sigma2;
        for jj=1:num_user
            a_kj = A(jj,:,kk).';
            p_j = p(jj);
            snr_jk = abs(a_kj.'*v)^2*p_j;
            denominator = denominator+snr_jk;
        end
        y(kk) = numerator/denominator;
    end

    %% update v
    v = update_v(v,A,y,gamma,1000);

    %% update y
    for kk=1:num_user
        a_kk =A(kk,:,kk).';
        p_k = p(kk);
        numerator = sqrt(p_k*alpha(kk)*(1+gamma(kk)))*(v'*conj(a_kk));
        denominator = sigma2;
        for jj=1:num_user
            a_kj = A(jj,:,kk).';
            p_j = p(jj);
            snr_jk = abs(a_kj.'*v)^2*p_j;
            denominator = denominator+snr_jk;
        end
        y(kk) = numerator/denominator;
    end
   obj(ii) = compute_rate(A, v, p, alpha, sigma2_dB);
   if ii>1 && obj(ii)-obj(ii-1)<1e-3
%        plot(obj);
       break;
   end
end

end
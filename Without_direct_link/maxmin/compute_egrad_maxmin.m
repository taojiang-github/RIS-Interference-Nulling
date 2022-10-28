function egrad_v = compute_egrad_maxmin(R1,R2,kk, x, alpha, sigma2_dB)
    sigma2 = 10^(sigma2_dB/10);
    tmp = squeeze(R1(kk,:,:));
    A_k = tmp*x/(x'*tmp*x+sigma2);
    tmp = squeeze(R2(kk,:,:));
    A_k = A_k-tmp*x/(x'*tmp*x+sigma2);
    egrad_v = 2*alpha(kk)*A_k;  
end
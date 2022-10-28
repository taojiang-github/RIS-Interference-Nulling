function egrad_v = compute_egrad(R1,R2, A, x, alpha,p, sigma2_dB)
    [ ~, tmp] = compute_minrate(A, x, p, alpha, sigma2_dB);
    kk = tmp.k;
    sigma2 = 10^(sigma2_dB/10);
    tmp = squeeze(R1(kk,:,:));
    A_k = tmp*x/(x'*tmp*x+sigma2);
    tmp = squeeze(R2(kk,:,:));
    A_k = A_k-tmp*x/(x'*tmp*x+sigma2);
    egrad_v = 2*alpha(kk)*A_k;  
end
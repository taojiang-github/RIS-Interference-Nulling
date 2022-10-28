function egrad_v = compute_egrad(R1,R2,c1,c2,d1,d2, num_user, x, alpha, sigma2_dB)
    sigma2 = 10^(sigma2_dB/10);
    egrad_v = 0;
    for kk =1:num_user
        tmp = squeeze(R1(kk,:,:));
        tmp2 = squeeze(c1(kk,:,:)).'; 
        A_k = (tmp*x+tmp2)/(x'*tmp*x+2*real(x'*tmp2)+d1(kk)+sigma2);
        tmp = squeeze(R2(kk,:,:));
        tmp2 = squeeze(c2(kk,:,:)).'; 
        A_k = A_k-(tmp*x+tmp2)/(x'*tmp*x+2*real(x'*tmp2)+d2(kk)+sigma2);
        egrad_v = 2*alpha(kk)*A_k+egrad_v;
    end    
end
function [v,obj_f] = subgradient_maxmin(A,sigma2_dB,v)
[num_user,num_irs_elements,~] = size(A);
p = ones(num_user,1);
alpha=ones(num_user,1);
R1 = zeros(num_user,num_irs_elements,num_irs_elements);
R2 = zeros(num_user,num_irs_elements,num_irs_elements);
for kk =1:num_user
    tmp=0;
    for jj=1:num_user
        a_kj = A(jj,:,kk)';
        tmp=tmp+a_kj*a_kj'*p(jj);
    end
    R1(kk,:,:) = tmp;
    a_kk = A(kk,:,kk)';
    tmp = tmp-a_kk*a_kk'*p(kk);
    R2(kk,:,:) = tmp;
end
iter_max = 10000;
rate_tract = nan(1,iter_max);
v_tract = nan(num_irs_elements,iter_max);
step_size = 0.5;  
no_increase = 0;
rate_best = 0;
for i_iter = 1:iter_max
    v_tract(:,i_iter) = v;
    [ rate_old, tmp] = compute_minrate(A, v, p, alpha, sigma2_dB);
    rate_tract(1,i_iter) = rate_old;
    if  rate_old > rate_best
        no_increase = 0;
        rate_best = rate_old;
    else
        no_increase = no_increase+1;
    end
    if no_increase>1000
        break
    end
    grad_v = compute_egrad_maxmin(R1,R2,tmp.k, v, alpha, sigma2_dB);
    v=v+step_size*grad_v/norm(grad_v);
    v = v./abs(v);
end

% plot(1:iter_max, rate_tract);
[obj_f, idx] = max(rate_tract);
v = v_tract(:,idx);
end
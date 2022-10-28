function [v,obj_f] = sca_func(A,sigma2_dB,v0)
[num_user,num_irs_elements,~] = size(A);
sigma2 = 10^(sigma2_dB/10);
p = ones(num_user,1);
alpha=ones(num_user,1);

%% initialization
iter_max_ini = 1000;
v_old = v0;
obj_v1_list = nan(1,iter_max_ini+1);
obj_v1_list(1) = obj_v1(A, v_old, sigma2_dB);
for i_iter = 1:iter_max_ini
    cvx_clear
    cvx_begin quiet
        variable v1(num_irs_elements,1) complex
        expression obj_list(num_user,1)
        for kk=1:num_user
            y_k = sigma2;
            tmp = 0;
            for jj=1:num_user
                if jj~=kk
                    a_kj = A(jj,:,kk).';
                    y_k = y_k+abs(a_kj.'*v_old)^2;
                    tmp = tmp+pow_abs(a_kj.'*v1,2);
                end
            end
            a_kk = A(kk,:,kk)';
            b_k = (a_kk'*v_old)/y_k;
            c_k = abs(v_old'*a_kk)^2/(y_k^2);
            rho_k = 2*real(b_k*v1'*a_kk)-c_k*tmp- sigma2*c_k;
            obj_list(kk) = rho_k;
        end
        obj = min(obj_list);
        maximize (obj)
        subject to
            for nn=1:num_irs_elements
               pow_abs(v1(nn),2)<=1;
            end
    cvx_end
    v_old = v1;
    obj_v1_list(1,i_iter+1) = obj_v1(A, v_old, sigma2_dB);
    if obj_v1_list(1,i_iter+1)-obj_v1_list(1,i_iter)<1e-3
       break; 
    end
end
% plot(obj_v1_list);

%% SCA
mu = 10*obj_v1_list(i_iter+1)/(1/norm(v_old,2)^2-1/num_irs_elements);
% mu = min(mu,1e8);
if mu>1e9
   scaling = mu/1e5;
else
   scaling = 1;
end
iter_max = 1000;
obj_v2_list = nan(1,iter_max+1);
obj_v2_list(1) = obj_v2(A, v_old, sigma2_dB, mu);
rate_tract = nan(1,1+iter_max);
rate_tract(1) = compute_minrate(A, v_old./abs(v_old), p, alpha, sigma2_dB);
for i_iter = 1:iter_max
    cvx_clear
    cvx_begin quiet
        variable v2(num_irs_elements,1) complex
        expression obj_list(num_user,1)
        expression l_theta
        for kk=1:num_user
            y_k = sigma2;
            tmp = 0;
            for jj=1:num_user
                if jj~=kk
                    a_kj = A(jj,:,kk).';
                    y_k = y_k+abs(a_kj.'*v_old)^2;
                    tmp = tmp+pow_abs(a_kj.'*v2,2);
                end
            end
            a_kk = A(kk,:,kk)';
            b_k = (a_kk'*v_old)/y_k;
            c_k = abs(v_old'*a_kk)^2/(y_k^2);
            rho_k = 2*real(b_k*v2'*a_kk)-c_k*tmp- sigma2*c_k;
            obj_list(kk) = rho_k;
        end
        l_theta= inv_pos(2*real(v_old'*v2)-norm(v_old,2)^2);
        obj = min(obj_list)+mu*(1/num_irs_elements-l_theta);
        maximize (obj/scaling)
        subject to
            (2*real(v_old'*v2)-norm(v_old,2)^2)>=0;
            for nn=1:num_irs_elements
               pow_abs(v2(nn),2)<=1;
            end
    cvx_end
    v_old = v2;
    rate_tract(1,i_iter+1) = compute_minrate(A, v_old./abs(v_old), p, alpha, sigma2_dB);
    obj_v2_list(1,i_iter+1) = obj_v2(A, v_old, sigma2_dB, mu);
    if obj_v2_list(1,i_iter+1)-obj_v2_list(1,i_iter)<1e-3
       break; 
    end
end
v = v_old./abs(v_old);
obj_f = compute_minrate(A, v, p, alpha, sigma2_dB);
end

function rate = obj_v1(A, v, sigma2_dB)
    sigma2 = 10^(sigma2_dB/10);
    [num_user, num_irs_elements, ~] = size(A);
    rate = Inf;
    for kk=1:num_user
        a_kk =A(kk,:,kk).';
        numerator = abs(a_kk.'*v)^2;
        denominator = sigma2;
        for jj=1:num_user
            if jj ~=kk
                a_kj = A(jj,:,kk).';
                snr_jk = abs(a_kj.'*v)^2;
                denominator = denominator+snr_jk;
            end
        end
        sinr_k = numerator/denominator;
        rate_k = log2(1+sinr_k);
        if rate_k<rate
           rate = rate_k;
        end
    end
end

function rate = obj_v2(A, v, sigma2_dB,mu)
    sigma2 = 10^(sigma2_dB/10);
    [num_user, num_irs_elements, ~] = size(A);
    rate = Inf;
    for kk=1:num_user
        a_kk =A(kk,:,kk).';
        numerator = abs(a_kk.'*v)^2;
        denominator = sigma2;
        for jj=1:num_user
            if jj ~=kk
                a_kj = A(jj,:,kk).';
                snr_jk = abs(a_kj.'*v)^2;
                denominator = denominator+snr_jk;
            end
        end
        sinr_k = numerator/denominator;
        rate_k = log2(1+sinr_k);
        if rate_k<rate
           rate = rate_k;
        end
    end
    rate = rate+mu*(1/num_irs_elements-1/norm(v,2)^2);
end
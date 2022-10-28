function [v,obj_f] = sdr_maxmin(A,sigma2_dB)
    sigma2 = 10^(sigma2_dB/10);
    [num_user,num_irs_elements,~] = size(A);
    R1 = {};
    R2 = {};
    for kk=1:num_user
       R1{kk} =  conj(A(kk,:,kk).')*A(kk,:,kk);
       tmp = 0;
       for jj=1:num_user
          if jj~=kk
             tmp = tmp+ conj(A(jj,:,kk).')*A(jj,:,kk);
          end
       end
       R2{kk} = tmp;
    end
    
    %%
    t_min=0;
    t_max=10000;
    feasible = 0;
    while true
        if t_max-t_min<=1e-3 && feasible
           break; 
        end
        t = (t_min+t_max)/2;
        %cvx_begin
        cvx_begin sdp quiet
            variable V(num_irs_elements,num_irs_elements)  hermitian semidefinite
            minimize 1
            subject to         
                for kk=1:num_user
                    real(t*trace(R2{kk}*V)+t*sigma2) <= real(trace(R1{kk}*V));
                end
                real(diag(V)) == 1;
        cvx_end
        %cvx_end
        if strcmp(cvx_status,'Solved')
            t_min = t;
            feasible = 1;
        else
            t_max = t;
            feasible = 0;
        end
    end
    obj_f.V = t;
    
    %% Obtain v
    if rank(V,1e-3) == 1
        [u,s,v_tmp] = svd(V);
        v = u(:,1)*sqrt(s(1,1));
        v = v./abs(v);
        obj_f.v = check_feasible_t(A,v,sigma2,R1,R2);
    else
        t_best = -1;
        for ii = 1:1000
            zi = chol(V,'lower');
            xi = (randn(num_irs_elements,1)+1i*randn(num_irs_elements,1))/sqrt(2);
            xi = zi*xi;
            xi = xi./abs(xi);
            t = check_feasible_t(A,xi,sigma2,R1,R2);
            if t > t_best
                v = xi;
                t_best = t;
            end
        end
        obj_f.v = t_best;
    end
end
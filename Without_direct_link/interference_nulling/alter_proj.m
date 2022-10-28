function [v,obj_f_all] = alter_proj(A)
    [num_user,num_irs_elements,~] = size(A);
    num_iter = 10000;
    obj_f = nan(num_iter+1,1);
    obj_f2 = nan(num_iter+1,1);
    A0=A;
    A = permute(A,[2,1,3]);
    A = reshape(A, [num_irs_elements,num_user*(num_user)]);
    indx = 1:num_user+1:num_user*num_user;
    A(:,indx) = [];

    v = randn(num_irs_elements,1)+1j*randn(num_irs_elements,1);
    v = v./abs(v);
   
    obj_f(1) = norm(A.'*v)^2;
    [signal_power,interfence_power, sir] = compute_SIR(A0, v);
    obj_f2(1) = max(1./sir);
    B = A.';
    A_pre = B'/( B* B')*B;
    for ii=1:num_iter
        v = v-A_pre*v;
        v(abs(v)<1e-10) = 1;
        v = v./abs(v);
        obj_f(ii) = norm(A.'*v)^2;
        [signal_power,interfence_power, sir] = compute_SIR(A0, v);
        obj_f2(ii) = max(1./sir);
        if ii>1 && ((obj_f(ii) < 1e-8) || (abs((obj_f(ii)-obj_f(ii-1))/obj_f(ii-1))<1e-3))
            break;
        end
%         if mod(ii,100)==1
%            fprintf('ii=%3d,obj=%4.3f\n',ii,obj_f(ii))
%         end
    end
%     semilogy(obj_f)
   obj_f_all.interference = obj_f;
   obj_f_all.interference2 = obj_f2;
end
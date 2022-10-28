function [v,obj_f] = alter_proj_2(A,C)
    [num_user,num_irs_elements,~] = size(A);
    num_iter = 10000;
    obj_f = nan(num_iter+1,1);
    A = permute(A,[2,1,3]);
    A = reshape(A, [num_irs_elements,num_user*(num_user)]);
    indx = 1:num_user+1:num_user*num_user;
    A_direct = A(:,indx);
    A(:,indx) = [];
    c = C(:);
    indx = 1:num_user+1:num_user*num_user;
    c(indx) = [];
    
    R = zeros(num_irs_elements,num_irs_elements);
    for kk=1:num_user
        R = R+conj(A_direct(:,kk))*A_direct(:,kk).';
    end
    [V,D] = eig(R);
    
    v = V(:,1);
    v = v./abs(v);
   
    obj_f(1) = norm(A.'*v+c)^2;
%     [signal_power,interfence_power, sir] = compute_SIR(A0,C0,v);
%     obj_f2(1) = max(1./sir);
    B = A.';
    A_pre = B'/( B* B')*B;
    c_pre = B'/( B* B')*c;
    for ii=1:num_iter
        v = v-(A_pre*v+c_pre);
        v(abs(v)<1e-10) = 1;
        v = v./abs(v);
        obj_f(ii) = norm(A.'*v+c)^2;
%         [signal_power,interfence_power, sir] = compute_SIR(A0, C0, v);
%         obj_f2(ii) = max(1./sir);
        if ii>1 && ((obj_f(ii) < 1e-8) || (abs((obj_f(ii)-obj_f(ii-1))/obj_f(ii-1))<1e-3))
            break;
        end
%         if mod(ii,100)==1
%            fprintf('ii=%3d,obj=%4.3f\n',ii,obj_f(ii))
%         end
    end
%     semilogy(obj_f)
%    obj_f_all.interference = obj_f;
%    obj_f_all.interference2 = obj_f2;
end
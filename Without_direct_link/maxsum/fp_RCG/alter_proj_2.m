function [v,obj_f] = alter_proj_2(A)
    [num_user,num_irs_elements,~] = size(A);
    num_iter = 10000;
    obj_f = nan(num_iter+1,1);
    A = permute(A,[2,1,3]);
    A = reshape(A, [num_irs_elements,num_user*(num_user)]);
    indx = 1:num_user+1:num_user*num_user;
    A_direct = A(:,indx);
    A(:,indx) = [];
    
    R = zeros(num_irs_elements,num_irs_elements);
    for kk=1:num_user
        R = R+conj(A_direct(:,kk))*A_direct(:,kk).';
    end
    [V,D] = eig(R);
    
    v = V(:,1);
    v = v./abs(v);
   
    obj_f(1) = norm(A.'*v)^2;
    B = A.';
    A_pre = B'/( B* B')*B;
%     A_pre = B'*B;
    for ii=1:num_iter
        v = v-A_pre*v;
        v(abs(v)<1e-10) = 1;
        v = v./abs(v);
        obj_f(ii) = norm(A.'*v)^2;
        if ii>1 && ((obj_f(ii) < 1e-8) || (abs((obj_f(ii)-obj_f(ii-1))/obj_f(ii-1))<1e-3))
            break;
        end
%         if mod(ii,100)==1
%            fprintf('ii=%3d,obj=%4.3f\n',ii,obj_f(ii))
%         end
    end
%     semilogy(obj_f)
end
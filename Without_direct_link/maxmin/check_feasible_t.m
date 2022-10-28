function t_best = check_feasible_t(A,v,sigma2,R1,R2)
[num_user,num_irs_elements,~] = size(A);
if nargin==3
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
end
V = v*v';
t_best = Inf;
for kk=1:num_user
    t = real(trace(R1{kk}*V))/real(trace(R2{kk}*V)+sigma2);
    if t<t_best
       t_best = t; 
    end
end

end
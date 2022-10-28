function [channel_irs_user,location_user_set]=generate_channel_Rician(num_irs_elements, num_user, irs_Nh, ...
    num_samples, Rician_factor, scale_factor,is_Tx)
location_irs = [0,0,0];
channel_irs_user = nan(num_samples,num_irs_elements,num_user);
location_user_set = nan(num_samples,num_user,3);
% tmp = 0;
for ii =1:num_samples
    location_user = generate_location(num_user,is_Tx);
    location_user_set(ii,:,:) = location_user;
    [pathloss, aoa_irs_y,aoa_irs_z] = generate_pathloss_angles(location_user,location_irs);
    pathloss = pathloss-scale_factor/2;
    pathloss = sqrt(10.^(-pathloss/10));
%     tmp = tmp+pathloss;
    i1 = mod(0:num_irs_elements-1,irs_Nh);
    i2 = floor((0:num_irs_elements-1)/irs_Nh);
    h_NLOS =(randn(num_irs_elements,num_user)+1j*randn(num_irs_elements,num_user))/sqrt(2);
    for kk=1:num_user
       a_irs_user = exp(1j*pi*(i1.*aoa_irs_y(kk)+i2.*aoa_irs_z(kk))).';
       channel_irs_user(ii,:,kk) =  sqrt(Rician_factor/(1+Rician_factor))*a_irs_user...
           +sqrt(1/(1+Rician_factor))*h_NLOS(:,kk);
       channel_irs_user(ii,:,kk) = channel_irs_user(ii,:,kk)*pathloss(kk);
    end
end
% tmp/num_samples
end

%% 
function location_user=generate_location(num_user,is_Tx)
location_user = nan(num_user,3);
if is_Tx == 1
    x_lim = [5,45];
    y_lim = [-45,-5];
    for kk=1:num_user
        x = x_lim(1)+(x_lim(2)-x_lim(1)) * rand;
        y = y_lim(1)+(y_lim(2)-y_lim(1)) * rand;
        location_user(kk,:) = [x,y,-20];
    end
elseif is_Tx == 2
    x_lim = [5,45];
    y_lim = [5,45];
    for kk=1:num_user
        x = x_lim(1)+(x_lim(2)-x_lim(1)) * rand;
        y = y_lim(1)+(y_lim(2)-y_lim(1)) * rand;
        location_user(kk,:) = [x,y,-20];
    end    
end
end

%%
function [pathloss, aoa_irs_y, aoa_irs_z] = generate_pathloss_angles(location_user,location_irs)
[num_user,~] = size(location_user);
aoa_irs_y = nan(num_user,1);
aoa_irs_z = nan(num_user,1);
pathloss = nan(num_user,1);
for kk=1:num_user
   d_k = norm(location_user(kk,:) - location_irs);
   pathloss(kk) = 30+22*log10(d_k);
   aoa_irs_y(kk) = (location_user(kk,2)-location_irs(2))/d_k;
   aoa_irs_z(kk) = (location_user(kk,3)-location_irs(3))/d_k;   
end
end
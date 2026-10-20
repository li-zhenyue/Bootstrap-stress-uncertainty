function sfm = SyntheticFocalMechanism( sigma, plane )
% 根据应力张量计算理论震源机制
%   sigma 应力张量，挤压为正
%   plane = [ strike, dip ] 为断层面参数

plane = plane * pi/180;
num = size( plane, 1 );

rake = zeros( num, 1 );
for i = 1 : num
    
    strike = plane( i, 1 ); dip = plane( i, 2 );
    n = [ - sin( strike ) * sin( dip ), cos( strike ) * sin( dip ), - cos( dip ) ]; % 断层面法向
    p = sigma * n';
    sigma_n = n * p; % 正应力大小
    shear = - ( p - sigma_n * n' ); % 剪应力矢量，负号是因为符号规定挤压为正
    sigma_s = norm( shear ); % 剪应力大小
    
    str_vec = [ cos( strike ), sin( strike ), 0 ];
    dip_vec = [ cos( dip ) * sin( strike ), - cos( dip ) * cos( strike ), - sin( dip ) ];
    rake( i ) = atan2( dip_vec * shear / sigma_s, str_vec * shear / sigma_s );
    
end

sfm = [ plane, rake ] * 180/pi;


end



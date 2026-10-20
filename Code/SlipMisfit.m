function misfit = SlipMisfit( stress_tensor, plane )
% 计算震源机制的滑动误差
% 输入的应力张量符号规定为：拉张为正

plane = plane * pi/180;
nbr = size( plane, 1 );
angle = nan( nbr, 1 );

for i = 1 : nbr
    
    strike = plane( i, 1 );
    dip = plane( i, 2 );
    rake = plane( i, 3 );
    
    % 剪应力方向
    n = [ -sin(strike) * sin(dip); cos(strike) * sin(dip); -cos(dip) ]; % 断层面法向
    p = stress_tensor * n;
    ns = n' * p; % normal stress
    ss = p - ns * n; % shear stress
    syn_slip = ss / norm(ss);
    
    % 震源机制给出的观测滑动矢量
    obs_slip = zeros( 1, 3 );
    obs_slip( 1 ) = cos(rake) * cos(strike) + cos(dip) * sin(rake) * sin(strike);
    obs_slip( 2 ) = cos(rake) * sin(strike) - cos(dip) * sin(rake) * cos(strike);
    obs_slip( 3 ) =-sin(rake) * sin(dip);
    
    angle( i ) = acos( obs_slip * syn_slip );
    
end

misfit = angle * 180/pi;


end



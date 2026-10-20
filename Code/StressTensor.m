function sigma = StressTensor( phi, delta, lamuda, R )
% 地理（北东下）坐标系下的应力张量，挤压为正
% phi: sigma1 轴走向，delta: sigma1 轴倾伏角
% lamuda: sigma3 轴在以 sigma1 轴为法向的平面内的滑动角
% R = ( s1 - s2 ) / ( s1 - s3 )

phi = phi * pi/180;
delta = delta * pi/180;
lamuda = lamuda * pi/180;

s1_vector = [ cos( delta ) * cos( phi ), cos( delta ) * sin( phi ), sin( delta ) ]; % sigma1 轴方向
s2_vector = [   sin( phi ) * sin( lamuda ) + cos( phi ) * sin( delta ) * cos( lamuda ),...
              - cos( phi ) * sin( lamuda ) + sin( phi ) * sin( delta ) * cos( lamuda ),...
              - cos( delta ) * cos( lamuda ) ]; % sigma2 轴方向
s3_vector = [ - sin( phi ) * cos( lamuda ) + cos( phi ) * sin( delta ) * sin( lamuda ),...
                cos( phi ) * cos( lamuda ) + sin( phi ) * sin( delta ) * sin( lamuda ),...
              - cos( delta ) * sin( lamuda ) ]; % sigma3 轴方向

s1 = 1; s2 = 1 - 2 * R; s3 = -1; % 挤压为正，拉张为负
M = [ s1_vector; s2_vector; s3_vector ];
sigma = M' * diag( [ s1, s2, s3 ] ) * M; 


end
function stress_tensor = linear_stress_inversion( mechanism )
% Michael 线性化反演应力张量的方法，输出的应力张量 [ 挤压为负 ]

strike = mechanism( :, 1 ) * pi/180;
dip = mechanism( :, 2 ) * pi/180;
rake = mechanism( :, 3 ) * pi/180;

% 法向和滑动方向
n1 = - sin( dip ) .* sin( strike );
n2 =   sin( dip ) .* cos( strike );
n3 = - cos( dip );

v1 =   cos( rake ) .* cos( strike ) + cos( dip ) .* sin( rake ) .* sin( strike );
v2 =   cos( rake ) .* sin( strike ) - cos( dip ) .* sin( rake ) .* cos( strike );
v3 = - sin( rake ) .* sin( dip );

% 构建系数矩阵 A
A11 = n1 .* ( n2.^2 + 2 * n3.^2 );
A12 = n2 .* ( 1 - 2 * n1.^2 );
A13 = n3 .* ( 1 - 2 * n1.^2 );
A14 = n1 .* ( - n2.^2 + n3.^2 );
A15 = - 2 * n1 .* n2 .* n3;

A21 = n2 .* ( - n1.^2 + n3.^2 );
A22 = n1 .* ( 1 - 2 * n2.^2 );
A23 = - 2 * n1 .* n2 .* n3;
A24 = n2 .* ( n1.^2 + 2 * n3.^2 );
A25 = n3 .* ( 1 - 2 * n2.^2 );

A31 = n3 .* ( - 2 * n1.^2 - n2.^2 );
A32 = - 2 * n1 .* n2 .* n3;
A33 = n1 .* ( 1 - 2 * n3.^2 );
A34 = n3 .* ( - n1.^2 - 2 * n2.^2 );
A35 = n2 .* ( 1 - 2 * n3.^2 );

N = length( strike ) * 3;
A = zeros( N, 5 );
A( 1 : 3 : N, 1 ) = A11; A( 1 : 3 : N, 2 ) = A12; A( 1 : 3 : N, 3 ) = A13; A( 1 : 3 : N, 4 ) = A14; A( 1 : 3 : N, 5 ) = A15;
A( 2 : 3 : N, 1 ) = A21; A( 2 : 3 : N, 2 ) = A22; A( 2 : 3 : N, 3 ) = A23; A( 2 : 3 : N, 4 ) = A24; A( 2 : 3 : N, 5 ) = A25;
A( 3 : 3 : N, 1 ) = A31; A( 3 : 3 : N, 2 ) = A32; A( 3 : 3 : N, 3 ) = A33; A( 3 : 3 : N, 4 ) = A34; A( 3 : 3 : N, 5 ) = A35;

% 滑动矢量矩阵 s
v = [ v1, v2, v3 ]';
s = reshape( v, N, 1 );

% 求解应力张量
stress_vector = real( pinv( A ) * s );
stress_vector( 6 ) = - stress_vector( 1 ) - stress_vector( 4 );
stress_tensor = [ stress_vector( 1 ) stress_vector( 2 ) stress_vector( 3 );
                  stress_vector( 2 ) stress_vector( 4 ) stress_vector( 5 );
                  stress_vector( 3 ) stress_vector( 5 ) stress_vector( 6 ) ];

stress_tensor = stress_tensor / norm( stress_tensor );


end



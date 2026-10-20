function [ ns, ss ] = faults_in_mohr_circle( stress_tensor, faults )
% 计算断层在莫尔圆中的位置，输入的应力张量挤压为正
% s1 和 S3 轴内侧的断层，剪应力为正，外侧的为负

% 断层法向
strike = faults( :, 1 ) * pi/180;
dip    = faults( :, 2 ) * pi/180;
n = [ - sin( dip ) .* sin( strike ), sin( dip ) .* cos( strike ), - cos( dip ) ];

[ V, D ] = eig( stress_tensor );
[ ~, ind ] = sort( diag( D ), 'descend' );

s1_vector = V( :, ind( 1 ) );
if s1_vector( 3 ) < 0
    s1_vector = -s1_vector;
end

s3_vector = V( :, ind( 3 ) );
if s3_vector( 3 ) < 0
    s3_vector = -s3_vector;
end

% 法向在 s1 和 s3 同侧的断层， 剪应力为负
symbol = ones( length( strike ), 1 );
ind = ( n * s1_vector ) .* ( n * s3_vector ) > 0;
symbol( ind ) = -1;

[ ns, ss ] = normal_shear_stress( stress_tensor, faults );
ss = ss .* symbol;


end


function [ ns, ss ] = normal_shear_stress( stress_tensor, faults )

strike = faults( :, 1 ) * pi/180;
dip    = faults( :, 2 ) * pi/180;
N = length( strike );

% 法向
n1 = - sin( dip ) .* sin( strike );
n2 =   sin( dip ) .* cos( strike );
n3 = - cos( dip );

ns = nan( N, 1 ); % normal stress
ss = nan( N, 1 ); % shear stress
for i = 1 : N
    
    normal = [ n1( i ), n2( i ), n3( i ) ];
    p = stress_tensor * normal';
    ns( i ) = normal * p;
    ss_vector = p - ns( i ) * normal';
    ss( i ) = norm( ss_vector );
    
end


end
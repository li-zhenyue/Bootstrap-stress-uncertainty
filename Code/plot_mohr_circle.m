function plot_mohr_circle( stress_tensor, faults, varargin )
% 绘制莫尔圆
% 至少需要输入主应力值或应力张量之一

p = inputParser;
p.KeepUnmatched = false;
p.CaseSensitive = false;
addParameter( p, 'CircleWidth', 1.2 )
addParameter( p, 'CircleColor', 'w' )
addParameter( p, 'IsoangleLine', { '-', 'LineWidth', 0.6, 'color', [0.86, 0.86, 0.86] } )
addParameter( p, 'CompressionSign', -1 ) % 默认输入应力张量的符号为【挤压为负】
addParameter( p, 'marker', { 'b+', 'LineWidth', 1, 'MarkerSize', 8 } )
parse( p, varargin{ : } )

stress_tensor = p.Results.CompressionSign * stress_tensor;

d = sort( eig( stress_tensor ), 'descend' );
s1 = d( 1 );
s2 = d( 2 );
s3 = d( 3 );

figure( 'color', 'w' )
axis equal
axis off
hold on

base_map( s1, s2, s3, p.Results.CircleWidth, p.Results.CircleColor, p.Results.IsoangleLine ) % 绘制三个莫尔圆

if ~isempty( faults )
    [ ns, ss ] = faults_in_mohr_circle( stress_tensor, faults );
    plot( ns, ss, p.Results.marker{ : } );
end


end


function base_map( s1, s2, s3, circle_width, color, IsoangleLine )

theta = linspace( 0, 2*pi, 100 );

% circle 13
r13 = ( s1 - s3 ) / 2; % 半径
o13 = ( s1 + s3 ) / 2; % 圆心
Xcir13 = r13 * cos( theta ) + o13; 
Ycir13 = r13 * sin( theta );
fill( Xcir13, Ycir13, color, 'EdgeColor', 'none' )

% circle 23
r23 = ( s2 - s3 ) / 2; 
o23 = ( s2 + s3 ) / 2;
Xcir23 = r23 * cos( theta ) + o23; 
Ycir23 = r23 * sin( theta );
fill( Xcir23, Ycir23, 'w', 'EdgeColor', 'none' )

% circle 12
r12 = ( s1 - s2 ) / 2; 
o12 = ( s1 + s2 ) / 2;
Xcir12 = r12 * cos( theta ) + o12; 
Ycir12 = r12 * sin( theta );
fill( Xcir12, Ycir12, 'w', 'EdgeColor', 'none' )

equal_angle_line( s1, s2, s3, IsoangleLine ) % 绘制等角度线
plot( Xcir13, Ycir13, 'k-', 'LineWidth', circle_width )
plot( Xcir23, Ycir23, 'k-', 'LineWidth', circle_width )
plot( Xcir12, Ycir12, 'k-', 'LineWidth', circle_width )
line( [ s1, s3 ], [ 0, 0 ], 'Color', 'k', 'LineWidth', circle_width )

% 标记3个主应力
fs = 15;
text( s1, 0, '\sigma_1', 'FontSize', fs )
text( s2, 0, '\sigma_2', 'FontSize', fs )
text( s3, 0, '\sigma_3', 'FontSize', fs )


end

function equal_angle_line( s1, s2, s3, LineDetail )
% 分别绘制与 s1 和 s3 轴夹角相同的断层面在莫尔圆中的投影
% LineDetail = { '-', 'LineWidth', 0.6, 'Color', [0.86, 0.86, 0.86] };

theta = ( 15 : 15 : 75 ) * pi/180;
n_theta = length( theta );
np = 100; % 弧线上的点数

% 与 s3 轴夹角相同的断层面
for i = 1 : n_theta
    n3_square = cos( theta( i ) )^2;
    n1_square = linspace( 0, 1 - n3_square, np );
    n2_square = 1 - n1_square - n3_square;
    ns = s1 * n1_square + s2 * n2_square + s3 * n3_square;
    ss = sqrt( s1^2 * n1_square + s2^2 * n2_square + s3^2 * n3_square - ns.^2 );
    plot( ns,  ss, LineDetail{ : } )
    plot( ns, -ss, LineDetail{ : } )
end

% 与 s1 轴夹角相同的断层面
for i = 1 : n_theta
    n1_square = cos( theta( i ) )^2;
    n3_square = linspace( 0, 1 - n1_square, np );
    n2_square = 1 - n1_square - n3_square;
    ns = s1 * n1_square + s2 * n2_square + s3 * n3_square;
    ss = sqrt( s1^2 * n1_square + s2^2 * n2_square + s3^2 * n3_square - ns.^2 );
    plot( ns,  ss, LineDetail{ : } )
    plot( ns, -ss, LineDetail{ : } )
end


end

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

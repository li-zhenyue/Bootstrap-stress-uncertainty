function faults_and_stress_axes_mapping( varargin )

p = inputParser;
p.KeepUnmatched = false;
p.CaseSensitive = false;
addParameter( p, 'radius', 1 )
addParameter( p, 'faults', nan )
addParameter( p, 'axes', nan )
addParameter( p, 'slip', 'off' )
addParameter( p, 'FaultsDetail', { '-', 'color',[0 0.4470 0.7410],'LineWidth', 0.8 } )
addParameter( p, 'PDetail', { 'ro', 'MarkerSize', 9 } )
addParameter( p, 'BDetail', { 'gd', 'MarkerSize', 9 } )
addParameter( p, 'TDetail', { 'b+', 'MarkerSize', 9 } )
addParameter( p, 'SlipDetail', { 'LineWidth', 1, 'color', 'r' , 'AutoScaleFactor', 0.2, 'MaxHeadSize', 1 } )
addParameter( p, 'CreatFigure', 'on' )
parse( p, varargin{ : } )


if strcmp( p.Results.CreatFigure, 'on' )
    % 绘制投影的赤道面
    figure( 'Color', 'w' )
    hold on
    theta = linspace( 0, 2 * pi, 200 );
    cirx = p.Results.radius * cos( theta );
    ciry = p.Results.radius * sin( theta );
    plot( cirx, ciry, 'k-', 'LineWidth', 1.2 )
%     plot( 0, 0, 'k+', 'LineWidth', 1.5, 'MarkerSize', 10 );
    marking_nswe( p.Results.radius )
end
hold on
axis image
axis off

if ~isnan( p.Results.faults )
    % 绘制断层面投影
    coord = faults_Schmidt_projection( p.Results.faults, p.Results.radius );
    for i = 1 : size( coord, 1 )
        plot( coord{ i, 1 }, coord{ i, 2 },  p.Results.FaultsDetail{ : } )
    end
    
    % 绘制滑动方向
    if strcmp( p.Results.slip, 'on' )
        [ x, y, vector_sign ] = slip_direction_projection( p.Results.faults,  p.Results.radius );
        dx = x ./ sqrt( x.^2 + y.^2 ) .* vector_sign;
        dy = y ./ sqrt( x.^2 + y.^2 ) .* vector_sign;
        quiver( x, y, dx, dy, p.Results.SlipDetail{ : } )
    end
end

% 绘制应力轴
if ~isnan( p.Results.axes )
    axes_column = size( p.Results.axes, 2 );
    switch axes_column
        case 2
            P = p.Results.axes;
            [ px, py ] = axes_Schmidt_projection( P, p.Results.radius );
            plot( px, py, p.Results.PDetail{ : } )
        case 4
            P = p.Results.axes( :, [ 1, 2 ] );
            T = p.Results.axes( :, [ 3, 4 ] );
            [ px, py ] = axes_Schmidt_projection( P, p.Results.radius );
            [ tx, ty ] = axes_Schmidt_projection( T, p.Results.radius );
            plot( px, py, p.Results.PDetail{ : } )
            plot( tx, ty, p.Results.TDetail{ : } )
        case 6
            P = p.Results.axes( :, [ 1, 2 ] );
            B = p.Results.axes( :, [ 3, 4 ] );
            T = p.Results.axes( :, [ 5, 6 ] );
            [ px, py ] = axes_Schmidt_projection( P, p.Results.radius );
            [ bx, by ] = axes_Schmidt_projection( B, p.Results.radius );
            [ tx, ty ] = axes_Schmidt_projection( T, p.Results.radius );
            plot( px, py, p.Results.PDetail{ : } )
            plot( bx, by, p.Results.BDetail{ : } )
            plot( tx, ty, p.Results.TDetail{ : } )
        otherwise
            error( 'axes azimuth and plunge input error!' )
    end
end


end


function coord = faults_Schmidt_projection( faults, radius )
% 断层下半球等面积投影的坐标

n_faults = size( faults, 1 );
n_rake = 100;

coord = cell( n_faults, 2 );
for i = 1 : n_faults
    strike = ones( n_rake, 1 ) * faults( i, 1 );
    dip    = ones( n_rake, 1 ) * faults( i, 2 );
    rake   = linspace( -180, 0, n_rake )';
    [ coord{ i, 1 }, coord{ i, 2 }, ~ ] = slip_direction_projection( [ strike, dip, rake ], radius );
end


end

function [ x, y, vector_sign ] = slip_direction_projection( sdr, radius )

sdr = sdr * pi/180;
strike = sdr( :, 1 );
dip    = sdr( :, 2 );
rake   = sdr( :, 3 );

% 断层走向为 x 轴，走向顺时针旋转 90°为 y 轴，竖直向下为 z 轴建立坐标系
slip_z = - sin( rake ) .* sin( dip );
slip_x =   cos( rake );
slip_y = - sin( rake ) .* cos( dip );

ih = acos( abs( slip_z ) ); % 离源角
strk = atan2( slip_y, slip_x ); % 断层坐标系下滑动方向的走向

ind = strk < 0;
strk( ind ) = strk( ind ) + pi;
slip_strike = strike + strk;

x = radius * sqrt( 2 ) * sin( ih / 2 ) .* sin( slip_strike );
y = radius * sqrt( 2 ) * sin( ih / 2 ) .* cos( slip_strike );

vector_sign = ones( length( strike ), 1 );
vector_sign( ind ) = -1;


end

function [ tx, ty ] = axes_Schmidt_projection( T, radius )
% project axis

Tstr = T( :, 1 ) * pi/180;
Tpl  = T( :, 2 ) * pi/180;

% coordinate of T axis
tx = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* sin( Tstr );
ty = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* cos( Tstr );


end

function marking_nswe( r )

line( [ 0, 0 ], [  r,  1.1 * r ], 'color', 'k', 'LineWidth', 1.2 ); text( 0,  1.16 * r, 'N', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' ) % N
line( [ 0, 0 ], [ -r, -1.1 * r ], 'color', 'k', 'LineWidth', 1.2 ); text( 0, -1.16 * r, 'S', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' ) % S
line( [ -1.1 * r, -r ], [ 0, 0 ], 'color', 'k', 'LineWidth', 1.2 ); text( -1.18*r,   0, 'W', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' ) % W
line( [ r,   1.1 * r ], [ 0, 0 ], 'color', 'k', 'LineWidth', 1.2 ); text(  1.16*r,   0, 'E', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' ) % E

end



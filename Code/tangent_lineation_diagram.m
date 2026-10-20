
function tangent_lineation_diagram( n, v, varargin )
% n 和 v 分别是断层面法向和滑动方向

[ mPlane, faults, quiverDetail, nswe, faultDetail, quiverLen,...
    createFigure, net ] = input_parser( varargin );

radius = 1;
nbr = size( n, 1 );

%%
if createFigure

    % figure( 'Color', 'w', 'Position', [730, 550, 580, 450] )
    hold on
    axis equal
    axis off

    % 投影赤道面
    circle_width = 0.8;
    color = [0.9, 0.9, 0.9];
    
    phi = linspace( 0, 2 * pi, 200 );
    circle = @( radius, phi, color, width ) plot( radius * cos( phi ), radius * sin( phi ),...
        'color', color, 'LineWidth', width );

    if net

        for r = 0.2 : 0.2 : 1
            if r < 1
                circle( r, phi, color, 0.7 * circle_width )
            end
        end

        cline = @( theta, radius, color, width ) line( sind(theta) * [-1, 1] * radius,...
            cosd(theta) * [-1, 1] * radius, 'LineStyle', '-', 'color', color, 'LineWidth', width );
        for theta = 0 : 30 : 150
            cline( theta, radius, color, 0.7 * circle_width );
        end

    end
    
    circle( radius, phi, 'k', circle_width )
    % scatter( 0, 0, 70, 'k', 'Marker', '+', 'LineWidth', 1.2 )

    % marking_nswe( radius, nswe )

end

%% 投影点与滑动方向
x  = nan( nbr, 1 );
y  = nan( nbr, 1 );
dx = nan( nbr, 1 );
dy = nan( nbr, 1 );

for i = 1 : nbr
    [ y( i ), x( i ), dy( i ), dx( i ) ] = stereoproj( n( i, : ), v( i, : ), quiverLen );
end

%% M-面、其他断层

if mPlane
    for i = 1 : nbr
        [ m_strike, m_dip ] = m_plane_parameters( n( i, : ), v( i, : ) );
        coord = faults_Schmidt_projection( [ m_strike, m_dip ], radius );
        plot( coord{ 1 }, coord{ 2 }, '-', 'color',[0 0.4470 0.7410],...
            'LineWidth', 0.8, faultDetail{:} )
    end
end

if ~isempty( faults )
    coord = faults_Schmidt_projection( faults, radius );
    for i = 1 : size( coord, 1 )
        plot( coord{ i, 1 }, coord{ i, 2 }, '-', 'color',[0 0.4470 0.7410],...
            'LineWidth', 0.8, faultDetail{:} )
    end
end

%% 绘制滑动方向箭头
quiver( x, y, dx, dy, 'AutoScale', 'off', ...
    'Alignment', 'center',...
    'Color', 'r', 'LineWidth', 0.8,...
    'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k',...
    'MaxHeadSize', 1, quiverDetail{:} )


end



function [ strike, dip ] = m_plane_parameters( n, v )

b = cross( n, v );
if b(3) > 0
    b = -b;
end

% 倾角
dip = acos( - b( 3 ) );

% 走向
strike = atan2( b( 2 ), b( 1 ) ) - pi/2;
if strike < 0
    strike = strike + 2 * pi; 
end

strike = strike * 180/pi;
dip = dip * 180/pi;

end


function [ mPlane, faults, quiverDetail, nswe, faultDetail, quiverLen,...
    createFigure, net ] = input_parser( varargin )

% 是否绘制 M-Plane
idx = find( strcmp( varargin{ : }, 'mPlane' ) );
if ~isempty( idx )
    mPlane = varargin{ : }{ idx + 1 };
else
    mPlane = false;
end

idx = find( strcmp( varargin{ : }, 'faults' ) );
if ~isempty( idx )
    faults = varargin{ : }{ idx + 1 };
else
    faults = [];
end

idx = find( strcmp( varargin{ : }, 'quiverDetail' ) );
if ~isempty( idx )
    quiverDetail = varargin{ : }{ idx + 1 };
else
    quiverDetail = {};
end

idx = find( strcmp( varargin{ : }, 'nswe' ) );
if ~isempty( idx )
    nswe = varargin{ : }{ idx + 1 };
else
    nswe = true;
end

idx = find( strcmp( varargin{ : }, 'faultDetail' ) );
if ~isempty( idx )
    faultDetail = varargin{ : }{ idx + 1 };
else
    faultDetail = {};
end

idx = find( strcmp( varargin{ : }, 'quiverLen' ) );
if ~isempty( idx )
    quiverLen = varargin{ : }{ idx + 1 };
else
    quiverLen = 0.15;
end

idx = find( strcmp( varargin{ : }, 'createFigure' ) );
if ~isempty( idx )
    createFigure = varargin{ : }{ idx + 1 };
else
    createFigure = true;
end

idx = find( strcmp( varargin{ : }, 'net' ) );
if ~isempty( idx )
    net = varargin{ : }{ idx + 1 };
else
    net = true;
end


end


function [ x, y, dx, dy ] = stereoproj( n, s, len )

n = -n;
s = -s;

k = 1 / sqrt( 1 + n(3) );

x = k * n( 1 );
y = k * n( 2 );

dx = k * s( 1 ) - 1/2 * x * k^3 * s( 3 ); 
dy = k * s( 2 ) - 1/2 * y * k^3 * s( 3 );

L = sqrt( dx^2 + dy^2 );

dx = dx / L * len;
dy = dy / L * len;


end


function marking_nswe( r, nswe )

line( [ 0, 0 ], [  r,  1.08 * r ], 'color', 'k', 'LineWidth', 1.2 );
line( [ 0, 0 ], [ -r, -1.08 * r ], 'color', 'k', 'LineWidth', 1.2 );
line( [ -1.08 * r, -r ], [ 0, 0 ], 'color', 'k', 'LineWidth', 1.2 );
line( [ r,   1.08 * r ], [ 0, 0 ], 'color', 'k', 'LineWidth', 1.2 );

if nswe
    text( 0,  1.16 * r, 'N', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' );
    text( 0, -1.16 * r, 'S', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' );
    text( -1.18*r,   0, 'W', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' );
    text(  1.16*r,   0, 'E', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'color', 'k', 'HorizontalAlignment', 'center' );
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

clc, clear

%%
critical_instability = 0.8;
friction = 0.6;
nr = 11;
R = linspace( 0, 1, nr );

%% 应力方向
phi = 30;
delta = 40;
lamuda = 40;

optimal = [ 30,	40,	157.45,	35.93, -88.34, 29.50 ];

%%
figure( 'color', 'w','Position', [ 400, 100, 1000, 870]  )
tiledlayout( 3, 4 );

for i = 1 : nr

    nexttile

    stress = StressTensor( phi, delta, lamuda, R( i ) );

    % 寻找合成震源机制的断层空间
    tangent_lineation_diagram( nan(1,3), nan(1,3) )

    fault_space( stress, friction, critical_instability )

    text( 0, 1.5, ['R = ',num2str( R( i ) )], 'FontSize', 15, 'FontName','TimesNewRoman',...
        'HorizontalAlignment', 'center' )

    % 合成震源机制
    faults = synthetic_faults( stress );
    plot_optimal_axes( optimal )

end



%%
function faults = synthetic_faults( stress )

optimal = [ 30,	40,	157.45,	35.93, -88.34, 29.50 ];
s1_azm = optimal( 1 );
s1_plg = optimal( 2 );
s3_azm = optimal( 5 );
s3_plg = optimal( 6 );

friction = 0.6;
critical_instability = 0.8;

nd = 8;
dip = linspace( 0, 90, nd );
ns = round( linspace(6, 36, nd ) );
faults = nan( sum(ns), 2 );
for i = 1 : nd

    strike = linspace( 0, 360, ns(i) )';
    
    k2 = sum( ns( 1 : i ) );
    k1 = k2 - ns( i ) + 1;
    faults( k1 : k2, : ) = [ strike, dip(i) * ones( ns(i), 1 ) ];
    
end


instability = FaultInstability( stress, faults, friction, 'CompressionSign', 'positive' );
idx = instability >= critical_instability;
faults = faults( idx, : );

sfm = SyntheticFocalMechanism( stress, faults );
[ ~, ~, outer_ind ] = clusteringFaults( s1_azm, s1_plg, s3_azm, s3_plg, sfm );

[ n, v ] = sdr2nv( sfm );
tangent_lineation_diagram( n( outer_ind, :), v( outer_ind, :), 'quiverLen', 0.12,...
    'quiverDetail', {'color', 'b', 'Marker', 'none', 'MarkerSize', 1.5,...
    'LineWidth', 0.8 }, 'createFigure', false )
tangent_lineation_diagram( n(~outer_ind, :), v(~outer_ind, :), 'quiverLen', 0.12,...
    'quiverDetail', {'color', 'r', 'Marker', 'none', 'MarkerSize', 1.5,...
    'LineWidth', 0.8 }, 'createFigure', false )


end


function plot_optimal_axes( optimal )

radius = 1;

s1_ap = optimal( 1 : 2 );
s2_ap = optimal( 3 : 4 );
s3_ap = optimal( 5 : 6 );

[ s1_x, s1_y ] = axes_Schmidt_projection( s1_ap, radius );
[ s2_x, s2_y ] = axes_Schmidt_projection( s2_ap, radius );
[ s3_x, s3_y ] = axes_Schmidt_projection( s3_ap, radius );

plot( s1_x, s1_y, 'ko', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 0.75 );
plot( s2_x, s2_y, 'kd', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 0.75 );
plot( s3_x, s3_y, 'k^', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 0.75 );

% dev_x = 0.25;
% dev_y = 0.1;
% text( s1_x + dev_x, s1_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
% text( s1_x + dev_x + 0.15, s1_y + dev_y - 0.12, '1', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
% 
% text( s2_x + dev_x, s2_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
% text( s2_x + dev_x + 0.15, s2_y + dev_y - 0.12, '2', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
% 
% text( s3_x + dev_x, s3_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
% text( s3_x + dev_x + 0.15, s3_y + dev_y - 0.12, '3', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )


end


function [ tx, ty ] = axes_Schmidt_projection( T, radius )
% project axis

Tstr = T( :, 1 ) * pi/180;
Tpl  = T( :, 2 ) * pi/180;

% coordinate of T axis
tx = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* sin( Tstr );
ty = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* cos( Tstr );


end


function fault_space( stress, friction, critical_instability )

strike = 0 : 1 : 360;
dip = 0 : 1 : 90;
ns = length( strike );
nd = length( dip );

[strike, dip] = meshgrid( strike, dip );
ins = FaultInstability( stress, [strike(:), dip(:)], friction, 'CompressionSign', 'positive' );

np = ns * nd;
n = nan( np, 3 );
n( :, 1 ) = -sind( dip(:) ) .* sind( strike(:) );
n( :, 2 ) =  sind( dip(:) ) .* cosd( strike(:) );
n( :, 3 ) = -cosd( dip(:) );
n = -n;

x = nan( np, 1 );
y = nan( np, 1 );
for j = 1 : np
    k = 1 / sqrt( 1 + n( j, 3 ) );
    y( j ) = k * n( j, 1 );
    x( j ) = k * n( j, 2 );
end

x = reshape( x, nd, ns );
y = reshape( y, nd, ns );
ins = reshape( ins, nd, ns );
% contour(x, y, ins, critical_instability * ones(1,2), 'LineWidth', 0.8, 'EdgeColor', 'g' )
contourf( x, y, ins, critical_instability * ones(1,2), 'FaceColor', [0.85, 0.93, 0.99], 'LineWidth', 0.2 )

end
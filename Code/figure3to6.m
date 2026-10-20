
% Set the last argument of "CreateSeismicFaults" on line 51 of the code below to 'both',...
% this reproduces Figures 3 and 4 of the manuscript. ...
% Setting it instead to 'inner' or 'outer' reproduces Figures 5 and 6.

clc; clear

%%
nfms = 50;
critical_instability = 0.8;
friction = 0.6;
realizations = 1000;
noise = 15;
nr = 11;
R = linspace( 0, 1, nr );
auxiliary_ratio = 0.5;

%% stress orientation
phi = 30;
delta = 40;
lamuda = 40;

%%
nbr = round( nfms * auxiliary_ratio );
nbins = 10;
R_range = nan( nr, 2 );
R_optimum = nan( nr, 1 );
sd_statistics = cell( nr, 1 );
R_statistics = cell( nr, 1 );

for j = 1 : nr
    
    stress = StressTensor( phi, delta, lamuda, R( j ) );
    
    % synthetic focal mechanisms
    % [ faults, instability ] = CreateSeismicFaults( stress, nfms, critical_instability, friction, 'both' );
    % noised_faults = nan( nfms, 3 );
    % for k = 1 : nfms
    %     noised_faults( k, : ) = noisy_mechanism( faults( k, : ), noise );
    % end
    % auxiliary = PlaneConvert( noised_faults );
    
    inverted_R = nan( realizations, 1 );
    inverted_s1 = nan( realizations, 2 );
    inverted_s2 = nan( realizations, 2 );
    inverted_s3 = nan( realizations, 2 );

    for i = 1 : realizations
    
        % synthetic focal mechanisms
        [ faults, instability ] = CreateSeismicFaults( stress, nfms, critical_instability, friction, 'both' );
        noised_faults = nan( nfms, 3 );
        for k = 1 : nfms
            noised_faults( k, : ) = noisy_mechanism( faults( k, : ), noise );
        end
        auxiliary = PlaneConvert( noised_faults );

        idx = false( nfms, 1 );
        idx( randperm( nfms, nbr ) ) = true;
        
        sdr = nan( nfms, 3 );
        sdr( idx, : ) = auxiliary( idx, : );
        sdr(~idx, : ) = noised_faults(~idx, : );
        
        inverted_stress = linear_stress_inversion( sdr ); % 挤压为负
        [ s1_vector, s2_vector, s3_vector, inverted_R( i ) ] = stress_tensor_decomposition( inverted_stress, 'CompressionSign', -1 );
        [ inverted_s1( i, 1 ), inverted_s1( i, 2 ) ] = vector_to_azimuth_and_plunge( s1_vector );
        [ inverted_s2( i, 1 ), inverted_s2( i, 2 ) ] = vector_to_azimuth_and_plunge( s2_vector );
        [ inverted_s3( i, 1 ), inverted_s3( i, 2 ) ] = vector_to_azimuth_and_plunge( s3_vector );
        
    end
    
    R_range( j, : ) = [ min( inverted_R ), max( inverted_R )];
    [ N, edges ] = histcounts( inverted_R, nbins );
    idx = find( N == max( N ) );
    R_optimum( j ) = mean( edges( [idx, idx+1] ) );
    sd_statistics{ j } = [ inverted_s1, inverted_s2, inverted_s3 ];
    R_statistics{ j } = inverted_R;
    
end

%%
figure( 'color', 'w','Position', [ 600, 300, 600, 500]  )
hold on
box on
grid on
scatter( R, R, 70, 's', 'MarkerFaceColor', [255  0   0]/255, 'MarkerEdgeColor', 'k' )

neg = R_optimum - R_range( :, 1 );
pos = R_range( :, 2 ) - R_optimum;
errorbar( R, R_optimum, neg, pos, 's', 'MarkerSize', 8, 'MarkerFaceColor',...
    'b', 'Color', 'b', 'MarkerEdgeColor', 'none', 'LineWidth', 1 )

axis( [ -0.1, 1.1, -0.1, 1.1 ] )
xlabel( 'Actual R value' )
ylabel( 'Inverted R value')
set( gca, 'position', [ 0.15, 0.18, 0.7, 0.7 ], 'FontSize', 14,...
    'FontName','TimesNewRoman','LineWidth', 1 )

%%    
figure( 'color', 'w','Position', [ 400, 100, 800, 800]  )
optimal = [ 30,	40,	157.45,	35.93, -88.34, 29.50 ];

tiledlayout( 3, 4 );
for j = 1 : nr
    nexttile
    
    StressStereonet( nan(1, 6), sd_statistics{ j },  'radius', 1, 'position', [ 0, 0 ],...
        'CircleWidth', 0.8, 'MarkerSize', 20, 'OptimalWidth', 2, 'CreatFigure', false )
    
    text( 0, 1.5, ['R = ',num2str( R( j ) )], 'FontSize', 15, 'FontName','TimesNewRoman',...
        'HorizontalAlignment', 'center' )

    plot_optimal_axes( optimal )

    axis off
end





%%
function plot_optimal_axes( optimal )

radius = 1;

s1_ap = optimal( 1 : 2 );
s2_ap = optimal( 3 : 4 );
s3_ap = optimal( 5 : 6 );

[ s1_x, s1_y ] = axes_Schmidt_projection( s1_ap, radius );
[ s2_x, s2_y ] = axes_Schmidt_projection( s2_ap, radius );
[ s3_x, s3_y ] = axes_Schmidt_projection( s3_ap, radius );

plot( s1_x, s1_y, 'ko', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 1 );
plot( s2_x, s2_y, 'kd', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 1 );
plot( s3_x, s3_y, 'k^', 'MarkerSize', 6.5, 'MarkerFaceColor', 'y', 'LineWidth', 1 );

dev_x = 0.25;
dev_y = 0.1;
text( s1_x + dev_x, s1_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
text( s1_x + dev_x + 0.15, s1_y + dev_y - 0.12, '1', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )

text( s2_x + dev_x, s2_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
text( s2_x + dev_x + 0.15, s2_y + dev_y - 0.12, '2', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )

text( s3_x + dev_x, s3_y + dev_y, '\sigma', 'FontName', 'TimesNewRoman', 'FontSize', 17,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
text( s3_x + dev_x + 0.15, s3_y + dev_y - 0.12, '3', 'FontName', 'TimesNewRoman', 'FontSize', 10,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )


end


function [ tx, ty ] = axes_Schmidt_projection( T, radius )
% project axis

Tstr = T( :, 1 ) * pi/180;
Tpl  = T( :, 2 ) * pi/180;

% coordinate of T axis
tx = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* sin( Tstr );
ty = radius * sqrt( 2 ) * sin( ( pi/2 - Tpl ) / 2 ) .* cos( Tstr );


end

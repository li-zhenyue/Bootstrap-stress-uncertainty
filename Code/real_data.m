
clc; clear

slickensides = load( 'TYM_sdr.dat' );
nbr = size( slickensides, 1 );
realizations = 2000;

auxiliary = PlaneConvert( slickensides );
statistics = nan( realizations, 7 );

for i = 1 : realizations
    
    idx = false( nbr, 1 );
    idx( randperm( nbr, round( nbr/2 ) ) ) = true;
    
    sdr = nan( nbr, 3 );
    sdr( idx, : ) = auxiliary( idx, : );
    sdr(~idx, : ) = slickensides(~idx, : );
    
    inverted_stress = linear_stress_inversion( sdr ); % 挤压为负
    statistics( i, : ) = tensor2apr( inverted_stress );
    
end

opt_ST = linear_stress_inversion( slickensides );
optimum = tensor2apr( opt_ST );

%%
figure( 'color', 'w' )
axis off
sd_statistics = statistics( :, 1 : 6 );
StressStereonet( optimum( 1 : 6 ), sd_statistics, 'radius', 1, 'position', [ 0, 0 ],...
        'CircleWidth', 1.2, 'MarkerSize', 60, 'OptimalWidth', 2,...
        'orientation', 1, 'CreatFigure', false )

%% R uncertity
figure( 'color', 'w' )
hold on
box on
R_statistics = statistics( :, 7 );
histogram( R_statistics, 'FaceColor', [30 144 255]/255, ...
    'BinWidth', 0.05, 'Normalization', 'probability' )
xlim( [0, 1] )
set( gca, 'xtick', 0 : 0.2 : 1, 'position', [0.15, 0.15, 0.7, 0.7],...
    'FontSize', 14, 'FontName', 'TimesNewRoman', 'LineWidth', 1,...
    'TickLength', [0.013, 0.006])
xlabel( 'R' )
ylabel( 'probability' )
plot( optimum(7) * ones( 1, 2 ), ylim, 'r-', 'LineWidth', 1.2 )
text( optimum(7) + 0.03, 0.25, 'Optimum', 'FontSize', 13, 'Color', 'r' )




%%
function apr = tensor2apr( ST )

[ s1_vector, s2_vector, s3_vector, shape_ratio ] = stress_tensor_decomposition( ST, 'CompressionSign', -1 );
[ s1_azm, s1_plg ] = vector_to_azimuth_and_plunge( s1_vector );
[ s2_azm, s2_plg ] = vector_to_azimuth_and_plunge( s2_vector );
[ s3_azm, s3_plg ] = vector_to_azimuth_and_plunge( s3_vector );

apr = [s1_azm, s1_plg, s2_azm, s2_plg, s3_azm, s3_plg, shape_ratio ];


end
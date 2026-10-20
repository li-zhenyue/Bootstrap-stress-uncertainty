
% 对擦痕数据执行应力反演
clc; clear

sdr = load( 'TYM_sdr.dat' );
realizations = 2000;
noise = 15;

m2st = @( m ) [m(1), m(2), m(3); m(2), m(4), m(5); m(3), m(5), m(6)];

[ optm, uncm ] = invert_slickensides( sdr, realizations, noise );
optimum = tensor2apr( m2st( optm ) );

statistics = nan( realizations, 7 );
for i = 1 : realizations
    statistics( i, : ) = tensor2apr( m2st( uncm(i, :) ) );
end


%% 绘图

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
histogram( R_statistics, 'FaceColor', [1, 0.41, 0.16], ...
    'BinWidth', 0.05, 'Normalization', 'probability' )
xlim( [0, 1] )
set( gca, 'xtick', 0 : 0.2 : 1, 'position', [0.15, 0.15, 0.7, 0.7],...
    'FontSize', 14, 'FontName', 'TimesNewRoman', 'LineWidth', 1,...
    'TickLength', [0.013, 0.006])
xlabel( 'R' )
ylabel( 'probability' )
plot( optimum( 7 ) * ones( 1, 2 ), ylim, '-', 'LineWidth', 1.2,'Color', [1, 0, 1] )
text( 0.68, 0.25, ['Optimum = ', num2str( round( optimum( 7 ), 2 ) )], 'FontSize', 13, 'Color', [1, 0, 1] )




%%

function [ optimum, uncertainty ] = invert_slickensides( sdr, realizations, noise )
% 通过随机扰动擦痕的法向和滑动方向，反演应力张量

nbr = size( sdr, 1 );
st2m = @(T) [T(1,1), T(1,2), T(1,3), T(2,2), T(2,3), T(3,3)];

ST = linear_stress_inversion( sdr );
optimum = st2m( ST );

uncertainty = nan( realizations, 6 );
for i = 1 : realizations
    
    noisy_sdr = nan( nbr, 3 );
    for j = 1 : nbr
        noisy_sdr( j, : )= noisy_mechanism( sdr( j, : ), noise );
    end

    ST = linear_stress_inversion( noisy_sdr );
    uncertainty( i, : ) = st2m( ST );

end


end


function apr = tensor2apr( ST )

[ s1_vector, s2_vector, s3_vector, shape_ratio ] = stress_tensor_decomposition( ST, 'CompressionSign', -1 );
[ s1_azm, s1_plg ] = vector_to_azimuth_and_plunge( s1_vector );
[ s2_azm, s2_plg ] = vector_to_azimuth_and_plunge( s2_vector );
[ s3_azm, s3_plg ] = vector_to_azimuth_and_plunge( s3_vector );

apr = [s1_azm, s1_plg, s2_azm, s2_plg, s3_azm, s3_plg, shape_ratio ];


end
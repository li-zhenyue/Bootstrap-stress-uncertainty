
% Set the data loaded on lines 6 and 7 to AVB_results.mat and AVB_sdr.dat, ...
% respectively, then run the program to obtain the results in Figure 7. ...
% Setting them instead to TYM_results.mat and TYM_sdr.dat gives the ...
% results in Figure 8.

clc; clear

load TYM_result.mat
slickensides = load( 'TYM_sdr.dat' );

%% R uncertity
figure( 'color', 'w' )
hold on
box on

R_ambugious = statistics_ambiguous_plane( :, 7 );
R_noisy = statistics_noisy_realizations( :, 7 );

histogram( R_ambugious, 'FaceColor', [30 144 255]/255, ...
    'BinWidth', 0.05, 'Normalization', 'probability' )
histogram( R_noisy, 'FaceColor', [1, 0.41, 0.16], ...
    'BinWidth', 0.05, 'Normalization', 'probability' )

xlim( [0, 1] )
% ylim( [0, 0.38] )
set( gca, 'xtick', 0 : 0.2 : 1, 'position', [0.15, 0.15, 0.7, 0.7],...
    'FontSize', 14, 'FontName', 'TimesNewRoman', 'LineWidth', 1,...
    'TickLength', [0.013, 0.006])
xlabel( 'R' )
ylabel( 'probability' )

ylim( [0, 0.7] )

plot( optimum( 7 ) * ones( 1, 2 ), ylim, '-', 'LineWidth', 1.2,'Color', [1, 0, 1] )
text( 0.67, 0.26, ['Optimum = ', num2str( round( optimum( 7 ), 2 ) )], 'FontSize', 13, 'Color', [1, 0, 1] )


%% data
faults_and_stress_axes_mapping( 'faults', slickensides, 'slip', 'on',...
    'SlipDetail', { 'LineWidth', 1, 'color', 'r' , 'AutoScaleFactor', 0.8, 'MaxHeadSize', 1 })

%% stress direction uncertainty
figure( 'color', 'w', 'Position', [300, 300, 650, 450] )
sd_ambugious = statistics_ambiguous_plane( :, 1 : 6 );
sd_noisy = statistics_noisy_realizations( :, 1 : 6 );

tiledlayout( 1, 2, 'Padding', 'compact' )

nexttile
axis off
StressStereonet( optimum( 1 : 6 ), sd_ambugious, 'radius', 1, 'position', [ 0, 0 ],...
        'CircleWidth', 1, 'MarkerSize', 60, 'OptimalWidth', 2,...
        'orientation', 1, 'CreatFigure', false )
text(0, -1.5, 'invert focal mechanisms', 'FontSize', 14, 'HorizontalAlignment', 'center' )

nexttile
axis off
StressStereonet( optimum( 1 : 6 ), sd_noisy, 'radius', 1, 'position', [ 0, 0 ],...
        'CircleWidth', 1, 'MarkerSize', 60, 'OptimalWidth', 2,...
        'orientation', 1, 'CreatFigure', false )
text(0, -1.5, 'invert slickensides', 'FontSize', 14, 'HorizontalAlignment', 'center' )


%% slip misfit

opt_ST_noisy = construct_stress_tensor( optimum(1), optimum(2), optimum(5), optimum(6), optimum(7) );
opt_ST_ambugious = construct_stress_tensor( optimum(1), optimum(2), optimum(5), optimum(6), 0.72 );
noisy_misfit = SlipMisfit( opt_ST_noisy, slickensides );
ambugious_misfit = SlipMisfit( opt_ST_ambugious, slickensides );

figure( 'color', 'w' )
hold on
box on

histogram( ambugious_misfit, 'FaceColor', [30 144 255]/255, ...
    'BinWidth', 5, 'Normalization', 'count' )
histogram( noisy_misfit, 'FaceColor', [1, 0.41, 0.16], ...
    'BinWidth', 5, 'Normalization', 'count' )

% xlim( [0, 1] )
set( gca, 'xtick', 0 : 10 : 140, 'position', [0.15, 0.15, 0.7, 0.7],...
    'FontSize', 14, 'FontName', 'TimesNewRoman', 'LineWidth', 1,...
    'TickLength', [0.013, 0.006])
xlabel( 'Slip misfit' )
ylabel( 'Count' )

legend( 'Bootstrap', 'Slickensides' )

ylim( [0, 14] )

%%
function ST = construct_stress_tensor( s1_azm, s1_plg, s3_azm, s3_plg, R )

s1_vector = [ cosd( s1_plg ) * cosd( s1_azm ), cosd( s1_plg ) * sind( s1_azm ), sind( s1_plg ) ]; % sigma1
s3_vector = [ cosd( s3_plg ) * cosd( s3_azm ), cosd( s3_plg ) * sind( s3_azm ), sind( s3_plg ) ]; % sigma3
s2_vector = cross( s1_vector, s3_vector );

s1 = -1; s2 = -1 + 2 * R; s3 = 1; % 拉张为正

M = [ s1_vector; s2_vector; s3_vector ];
ST = M' * diag( [ s1, s2, s3 ] ) * M; 

end
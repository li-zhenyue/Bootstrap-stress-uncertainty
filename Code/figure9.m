
clc; clear

% AVB
% slickensides = load( 'AVB_sdr.dat' );
% s1_azm =  76.6650;
% s1_plg =  74.0956;
% s3_azm = -27.0149;
% s3_plg =  3.8552;
% optmum_idx = 27;

% TYM
slickensides = load( 'TYM_sdr.dat' );
s1_azm = -132.2141;
s1_plg =  83.7797;
s3_azm = -24.9401;
s3_plg =  1.8537;
optmum_idx = 37;

R = 0 : 0.025 : 1;
nr = length( R );

sm = nan( 1, nr );
unc = nan( 1, nr );
for i = 1 : nr
    ST = construct_stress_tensor( s1_azm, s1_plg, s3_azm, s3_plg, R( i ) );
    misfit = SlipMisfit( -ST, slickensides );
    sm( i ) = mean( misfit );
    unc( i ) = std( misfit, 1 );
end

%%
figure( 'Color', 'w' )
hold on
box on
grid on
fill( [ R, flip(R)], [ sm + unc, flip( sm - unc ) ], [0.8, 0.8, 0.8], 'EdgeColor','none', 'FaceAlpha', 0.5 )
plot( R, sm, 'k-', 'LineWidth', 1.2 )
set( gca, 'position', [0.15, 0.15, 0.7, 0.7], 'FontSize', 14, 'LineWidth', 1, 'TickLength', [0.013, 0.006] )
scatter( R( optmum_idx ), sm( optmum_idx ), 100, 'Marker', '^', 'MarkerFaceColor', [1, 0.41, 0.16], 'MarkerEdgeColor', 'none' )

scatter( R( 30 ), sm( 30 ), 130, 'Marker', 'square', 'MarkerFaceColor', [30 144 255]/255, 'MarkerEdgeColor', 'none' )


%%
function sigma = construct_stress_tensor( s1_azm, s1_plg, s3_azm, s3_plg, R )

stressDirection = @( phi, delta ) [ cosd( delta ) * cosd( phi ), cosd( delta ) * sind( phi ), sind( delta ) ];
s1_vector = stressDirection( s1_azm, s1_plg );
s3_vector = stressDirection( s3_azm, s3_plg );
s2_vector = cross( s1_vector, s3_vector );
s2_vector = s2_vector / norm( s2_vector );

s1 = 1; s2 = 1 - 2 * R; s3 = -1; % 挤压为正，拉张为负
M = [ s1_vector; s2_vector; s3_vector ];
sigma = M' * diag( [ s1, s2, s3 ] ) * M; 


end
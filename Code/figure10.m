
% Comparison of slip misfit on the auxiliary fault planes for different R values

clc; clear

%%
nfms = 6000;
critical_instability = 0.8;
friction = 0.6;
noise = 0;
R = 0 : 0.01 : 1 ;
nr = length( R );

% stress orientation
phi = 30;
delta = 40;
lamuda = 40;

%%
slip_misfit = nan( nr, 2 );

for i = 1 : nr
     
    stress = StressTensor( phi, delta, lamuda, R( i ) );
    faults = CreateSeismicFaults( stress, nfms, critical_instability, friction, 'both' );
    noisy_faults = nan( nfms, 3 );
    for k = 1 : nfms
        noisy_faults( k, : ) = noisy_mechanism( faults( k, : ), noise );
    end
    auxiliary = PlaneConvert( noisy_faults );

    angle = SlipMisfit( -stress, auxiliary );
    slip_misfit( i, 1 ) = mean( real( angle ) );
    slip_misfit( i, 2 ) = std( real( angle ), 1 );

end

%%
sm = smooth( slip_misfit( :, 1 ), 5 );
unc = smooth( slip_misfit( :, 2 ), 5 );

figure( 'Color', 'w' )
hold on
box on
grid on
fill( [ R, flip( R ) ], [ (sm + unc)', flip( (sm - unc)' ) ], [0.8, 0.8, 0.8], 'EdgeColor','none', 'FaceAlpha', 0.5 )
plot( R, sm, 'k-', 'LineWidth', 1.2 )
set( gca, 'position', [0.15, 0.15, 0.7, 0.7], 'FontSize', 14, 'LineWidth', 1, 'TickLength', [0.013, 0.006] )

ylim( [-5, 50] )


clc; clear; close all

% Stress parameters
phi = 30; 
delta = 40; 
lamuda = 40; 
R = 0.6;
ST = StressTensor( phi, delta, lamuda, R ); % Stress tensor

% the azimuth and plunge of the principal stress axes,...
% as well as the shape ratio
apr = tensor2apr( -ST );

nbr = 50; % The number of synthetic faults
critical_instability = 0.93;
friction = 0.6;

% Synthetic faults on either side of the σ₁ axis
inner = CreateSeismicFaults( ST, nbr, critical_instability, friction, 'inner' );
outer = CreateSeismicFaults( ST, nbr, critical_instability, friction, 'outer' );

% Obtain the parameters of the auxiliary fault planes and of the P, B and T axes
[ inner_aux, inner_P, inner_B, inner_T ] = PlaneConvert( inner );
[ outer_aux, outer_P, outer_B, outer_T ] = PlaneConvert( outer );


%% Mohr circle

plot_mohr_circle( ST, [], 'CompressionSign', 1 );
[ inner_ns, inner_ss ] = faults_in_mohr_circle( ST, inner );
[ outer_ns, outer_ss ] = faults_in_mohr_circle( ST, outer );
plot( inner_ns, inner_ss, 'r+', 'MarkerSize', 7, 'LineWidth', 1 )
plot( outer_ns, outer_ss, 'b+', 'MarkerSize', 7, 'LineWidth', 1 )
plot_Mohr_Coulomb_failure_line( 1, -1, critical_instability, friction, -1.1, 0.25,...
    'LineDetail', { 'color', 'k', 'LineWidth', 1 } )


%% Plot the fault nodal planes and the principal stress axes 

faults_and_stress_axes_mapping( 'faults', inner_aux,...
    'FaultsDetail', { '-', 'color', 0.75 * ones(1,3), 'LineWidth', 0.8 } )

faults_and_stress_axes_mapping( 'faults', outer_aux,...
    'FaultsDetail', { '-', 'color', 0.75 * ones(1,3), 'LineWidth', 0.8 }, ...
    'CreatFigure', 'off')

faults_and_stress_axes_mapping( 'faults', inner, ...
    'FaultsDetail', { '-', 'color', 'r', 'LineWidth', 0.8 }, ...
    'CreatFigure', 'off')

faults_and_stress_axes_mapping( 'faults', outer,...
    'FaultsDetail', { '-', 'color', 'b', 'LineWidth', 0.8 },...
    'CreatFigure', 'off' )

sz = 13;
faults_and_stress_axes_mapping( 'axes', apr(1 : 6), ...
    'PDetail', { 'ko', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'BDetail', { 'kd', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'TDetail', { 'k^', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'CreatFigure', 'off'  )


%% P, B, and T axes

faults_and_stress_axes_mapping( 'axes', [inner_P, inner_B, inner_T], ...
    'PDetail', { 'ro', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'BDetail', { 'rd', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'TDetail', { 'r+', 'MarkerSize', 8, 'LineWidth', 1 } )

faults_and_stress_axes_mapping( 'axes', [outer_P, outer_B, outer_T], ...
    'PDetail', { 'bo', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'BDetail', { 'bd', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'TDetail', { 'b+', 'MarkerSize', 8, 'LineWidth', 1 }, ...
    'CreatFigure', 'off' )

faults_and_stress_axes_mapping( 'axes', apr(1 : 6), ...
    'PDetail', { 'ko', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'BDetail', { 'kd', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'TDetail', { 'k^', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'CreatFigure', 'off'  )

%% The distribution of the P, B and T axes after adding noise

noise = 15;
noisy_inner = nan( nbr, 3 );
noisy_outer = nan( nbr, 3 );
for i = 1 : nbr
    noisy_inner( i, : ) = noisy_mechanism( inner( i, : ), noise );
    noisy_outer( i, : ) = noisy_mechanism( outer( i, : ), noise );
end

[ ~, noisy_inner_P, noisy_inner_B, noisy_inner_T ] = PlaneConvert( noisy_inner );
[ ~, noisy_outer_P, noisy_outer_B, noisy_outer_T ] = PlaneConvert( noisy_outer );

faults_and_stress_axes_mapping( 'axes', [noisy_inner_P, noisy_inner_B, noisy_inner_T], ...
    'PDetail', { 'ro', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'BDetail', { 'rd', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'TDetail', { 'r+', 'MarkerSize', 8, 'LineWidth', 1 } )

faults_and_stress_axes_mapping( 'axes', [noisy_outer_P, noisy_outer_B, noisy_outer_T], ...
    'PDetail', { 'bo', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'BDetail', { 'bd', 'MarkerSize', 7, 'LineWidth', 1 }, ...
    'TDetail', { 'b+', 'MarkerSize', 8, 'LineWidth', 1 }, ...
    'CreatFigure', 'off' )

faults_and_stress_axes_mapping( 'axes', apr(1 : 6), ...
    'PDetail', { 'ko', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'BDetail', { 'kd', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'TDetail', { 'k^', 'MarkerSize', sz, 'MarkerFaceColor', 'y', 'LineWidth', 1.2 }, ...
    'CreatFigure', 'off'  )



%%
function apr = tensor2apr( ST )

[ s1_vector, s2_vector, s3_vector, shape_ratio ] = stress_tensor_decomposition( ST, 'CompressionSign', -1 );
[ s1_azm, s1_plg ] = vector_to_azimuth_and_plunge( s1_vector );
[ s2_azm, s2_plg ] = vector_to_azimuth_and_plunge( s2_vector );
[ s3_azm, s3_plg ] = vector_to_azimuth_and_plunge( s3_vector );

apr = [s1_azm, s1_plg, s2_azm, s2_plg, s3_azm, s3_plg, shape_ratio ];


end
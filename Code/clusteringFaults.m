
function [ inner, outer, ind ] = clusteringFaults( s1_azm, s1_plg, s3_azm, s3_plg, faults )


strike = faults( :, 1 );
dip = faults( :, 2 );

fnormal = @(phi, delta) [ -sind(delta) .* sind(phi), sind(delta) .* cosd(phi), -cosd(delta) ]; % fault normal
sdirection = @(azm, plg) [ cosd(plg) * cosd(azm), cosd(plg) * sind(azm), sind(plg) ]; % principal stress direction

s1_direction = sdirection( s1_azm, s1_plg );
s3_direction = sdirection( s3_azm, s3_plg );
normal = fnormal( strike, dip );

ind = ( normal * s1_direction' ) .* ( normal * s3_direction' ) > 0;
outer = faults( ind, : );
inner = faults(~ind, : );


end
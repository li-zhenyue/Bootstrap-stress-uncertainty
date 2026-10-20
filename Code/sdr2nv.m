
function [ n, v ] = sdr2nv( sdr )

sdr = sdr * pi/180;

strike = sdr( :, 1 );
dip = sdr( :, 2 );
rake = sdr( :, 3 );

nbr = size( sdr, 1 );

n = nan( nbr, 3 );
n( :, 1 ) = -sin( dip ) .* sin( strike );
n( :, 2 ) =  sin( dip ) .* cos( strike );
n( :, 3 ) = -cos( dip );

v = nan( nbr, 3 );
v( :, 1 ) =  cos( rake ) .* cos( strike ) + cos( dip ) .* sin( rake ) .* sin( strike );
v( :, 2 ) =  cos( rake ) .* sin( strike ) - cos( dip ) .* sin( rake ) .* cos( strike );
v( :, 3 ) = -sin( rake ) .* sin( dip );


end
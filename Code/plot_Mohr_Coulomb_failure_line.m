function plot_Mohr_Coulomb_failure_line( sigma1, sigma3, instability, friction, ns1, ns2, varargin )

p = inputParser;
p.KeepUnmatched = false;
p.CaseSensitive = false;
addParameter( p, 'LineDetail', { 'color', 'k', 'LineWidth', 1.2 } )
addParameter( p, 'branch', 'bilateral' ) % or unilateral
parse( p, varargin{ : } )

r = ( sigma1 - sigma3 ) / 2;
x0 = sigma1 - r * ( 1 + sqrt( 1 + friction^2 ) / friction ) * instability;

mohr_coulomb_line = @( ns ) friction * ( ns - x0 );
ss1 = mohr_coulomb_line( ns1 );
ss2 = mohr_coulomb_line( ns2 );

if strcmp( p.Results.branch, 'unilateral' )
    line( [ ns1, ns2 ],   [ ss1, ss2 ], p.Results.LineDetail{ : } )
elseif strcmp( p.Results.branch, 'bilateral' )
    line( [ ns1, ns2 ],   [ ss1, ss2 ], p.Results.LineDetail{ : } )
    line( [ ns1, ns2 ], - [ ss1, ss2 ], p.Results.LineDetail{ : } )
else
    error( 'parameter input error !')
end


end



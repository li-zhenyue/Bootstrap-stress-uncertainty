function [ faults, instability ] = CreateSeismicFaults( stress_tensor, number, critical_instability, friction, clusters )
% 输入的应力张量挤压为正

[ V, D ] = eig( stress_tensor );
[ ~, idx ] = sort( diag( D ), 'descend' );

s1_vector = V( :, idx( 1 ) );
if s1_vector( 3 ) < 0
    s1_vector = - s1_vector;
end

s3_vector = V( :, idx( 3 ) );
if s3_vector( 3 ) < 0
    s3_vector = - s3_vector;
end

normal = @(phi, delta) [ -sind(delta) .* sind(phi), sind(delta) .* cosd(phi), -cosd(delta) ];

% 生成一定数量满足失稳条件的断层
num = 500;
while true
    
    strike = rand( num, 1 ) * 360; 
    dip = rand( num, 1 ) * 90;
    instability = FaultInstability( stress_tensor, [ strike, dip ], friction, 'CompressionSign', 'positive' );
    
    % failure criterion 
    indf = instability >= critical_instability;

    % clusters
    if strcmp( clusters, 'outer' )
        indc = ( normal( strike, dip ) * s1_vector ) .* ( normal( strike, dip ) * s3_vector ) > 0;
    elseif strcmp( clusters, 'inner' )
        indc = ( normal( strike, dip ) * s1_vector ) .* ( normal( strike, dip ) * s3_vector ) < 0;
    elseif strcmp( clusters, 'both' )
        indc = true( num, 1 );
    else
        error( 'clusters code error !' )
    end

    % 符合标准的断层
    ind = indf & indc;
    strike = strike( ind ); 
    dip = dip( ind ); 
    instability = instability( ind );
    
    if length( strike ) < number
        num = num * 2;
    else
        r = randperm( length( strike ), number );
        strike = strike( r );
        dip = dip( r );
        instability = instability( r );
        break
    end
    
end

faults = SyntheticFocalMechanism( stress_tensor, [ strike, dip ] ); % 应力挤压为正


end



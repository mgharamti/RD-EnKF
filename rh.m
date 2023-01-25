function Histo = rh(dXf, refX)

    r = zeros( size( dXf,1 ),1 );
    for i = 1:size( dXf,1 )
        % sorting the forecasts and obs 
        s= sort( [ dXf(i,:),refX(i) ] );      
        % finding the position of the obs
        r(i)= find( s==refX(i) );              
    end
    ranks= 1:size( dXf,2 )+1;
    Histo= histc( r,ranks )/size( dXf,1 ); 
    
end

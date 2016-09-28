function Rxx_new = adjustMapSize(Rxx,Nr)

% Size of map, new size must be Nr
N = size(Rxx,1);

if N == Nr 
    Rxx_new = Rxx;
    return
end

Rxx_new = nan(Nr);

if N > Nr
    diff = (N - Nr) / 2;
    
    diff= ceil(diff);
    
    for ii = 1:Nr
        
        for jj = 1:Nr
            
            ii=ceil(ii);
            jj=ceil(jj);
            
            Rxx_new(ii,jj) = Rxx(ii+diff,jj+diff);
           
        end
    end
    
end

if N < Nr
    diff = (Nr - N) / 2;
    for ii = 1:N
        for jj = 1:N
            Rxx_new(ii+diff,jj+diff) = Rxx(ii,jj);
        end
    end
end
disp('')
function M= turnNansIntoLocalMeans(M)

M_orig= M;

%pad with edges
[size_M]= size(M);
M= zeros(size_M+2); 
M(2:end-1, 2:end-1)= M_orig;

[i j]= find(isnan(M));

for len= 1: length(i)
    M(i(len),j(len))= nanmean([M(i(len)+1,j(len)) M(i(len)-1,j(len)) ...
        M(i(len),j(len)+1) M(i(len),j(len)-1)]);
end

%remove padding
M_new= M;
M= nan(size_M);
M= M_new(2:end-1, 2:end-1);
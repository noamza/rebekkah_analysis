function gaussian_mat= createGaussianMat(mat, circle_size, centers, centervalue)

gsize= size(mat);

gaussian_mat= zeros(gsize);

for h= 1:length(centers)
gaussian_mat(centers(h,1), centers(h,2))= centervalue*25;
end

gaussian_mat_orig= gaussian_mat;
gaussian_mat= zeros(gsize + 20); % pad with edges
new_size= size(gaussian_mat);
gaussian_mat(11:new_size(1)-10, 11:new_size(2)-10)= gaussian_mat_orig; % gaussian mat with padded edges 


%kernalsize= sigma;

% for r=:kernalsize
%     for c=1:kernalsize
%         for m=1:length(centers)
%         gaussian_mat(r,c) = gaussC(r,c, sigma, centers(m,:));
%         end
%     end
% end

parms.sigma= 2;

% if mod(parms.sigma,2)~= 0
%     parms.sigma= parms.sigma+1;
% end

kernal_size= round(circle_size*2) +1;
gaussian_mat= SmoothGaussian(gaussian_mat, parms, kernal_size);
gaussian_mat= gaussian_mat(11:new_size(1)-10, 11:new_size(2)-10);

disp('');
 
% function val = gaussC(x, y, sigma, center)
% xc = center(1);
% yc = center(2);
% exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
% val       = (exp(-exponent));
%  end


end

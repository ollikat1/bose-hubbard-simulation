function [ coeff_mat ] = coefficients(N,p)
% Creates the coefficient matrix for N bosons, p lattices, describing the
% possibilities of distributing the bosons on lattices.
%
% Example for N=2, p=3:
%
% coeff_mat = [0, 0, 2
%              0, 1, 1
%              0, 2, 0
%              1, 0, 1
%              1, 1, 0
%              2, 0, 0]
%
% For fixed p, one could do this (example p=5):
%
% i=1;
% for n_5=0:N,
%     for n_4=0:N-n_5,
%         for n_3=0:N-n_5-n_4,
%             for n_2=0:N-n_5-n_4-n_3,
%                 for n_1=0:N-n_5-n_4-n_3-n_2,
%                     if n_1+n_2+n_3+n_4+n_5==N,
%                         coeff_vec(i,1) = n_1;
%                         coeff_vec(i,2) = n_2;
%                         coeff_vec(i,3) = n_3;
%                         coeff_vec(i,4) = n_4;
%                         coeff_vec(i,5) = n_5;
%                         i=i+1;
%                     end
%                 end
%             end
%         end
%     end
% end
%
% However, we want to pass p, i.e., the number of for-loops, also as 
% variable argument. This could be done with a recursive function, for example.
%
% Instead, we follow a much faster trick given in
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/143037

c = nchoosek(1:N+p-1,p-1);
m = size(c,1);
t = ones(m,N+p-1);
t(repmat((1:m).',1,p-1)+(c-1)*m) = 0;
u = [zeros(1,m);t.';zeros(1,m)];
v = cumsum(u,1);
coeff_mat = diff(reshape(v(u==0),p+1,m),1).';

end
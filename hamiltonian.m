function [ H_diag, H_offdiag ] = hamiltonian(coeffs, type )
% Returns the Hamiltonian matrix for the system
% type=1: 1D lattice without periodic boundary conditions
% type=2: 1D lattice with periodic boundary conditions
basis_size=length(coeffs);
lattice_size=length(coeffs(1,:));
H_diag=zeros(basis_size,basis_size);
H_offdiag=zeros(basis_size,basis_size);

if type==1 || type==2, 
    % Diagonal elements
    for i=1:basis_size,
        for j=1:lattice_size,
            if coeffs(i,j)>1, % Only nonzero if occupancy 2 or greater.
                H_diag(i,i) = H_diag(i,i) + coeffs(i,j);
            end
        end
    end
    % Off-diagonal site hopping, nearest neighbours with adjacent indices in 1D lattice
    for i=1:basis_size,
        for j=i+1:basis_size,
            for k=1:lattice_size-1,
                if coeffs(i,k)==coeffs(j,k)+1 && coeffs(i,k+1)==coeffs(j,k+1)-1,
                    H_offdiag(i,j) = H_offdiag(i,j) + sqrt(coeffs(i,k)*coeffs(j,k+1));
                    H_offdiag(j,i) = H_offdiag(i,j);
                end
                if coeffs(i,k)==coeffs(j,k)-1 && coeffs(i,k+1)==coeffs(j,k+1)+1,
                    H_offdiag(i,j) = H_offdiag(i,j) + sqrt(coeffs(i,k+1)*coeffs(j,k));
                    H_offdiag(j,i) = H_offdiag(i,j);
                end
            end
            % Periodic boundaries: connect site 1 and "end"
            if type==2,
                if coeffs(i,1)==coeffs(j,1)+1 && coeffs(i,end)==coeffs(j,end)-1,
                    H_offdiag(i,j) = H_offdiag(i,j) + sqrt(coeffs(i,1)*coeffs(j,end));
                    H_offdiag(j,i) = H_offdiag(i,j);
                end
                if coeffs(i,1)==coeffs(j,1)-1 && coeffs(i,end)==coeffs(j,end)+1,
                    H_offdiag(i,j) = H_offdiag(i,j) + sqrt(coeffs(i,end)*coeffs(j,1));
                    H_offdiag(j,i) = H_offdiag(i,j);
                end
            end
        end
    end
end

% Make the matrices sparse

H_diag = sparse(H_diag);
H_offdiag = sparse(H_offdiag);

end
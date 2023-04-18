clear all
close all

%  1D lattice with/without periodic boundary conditions

J=1; % Site hopping strength
N=3; % Number of bosons
p=5; % Number of lattice points
ratio_max=100; % Max ratio of U/J

colors = {'r','g','b','m','k'}; % Color cell used in plotting

coeffs = coefficients(N,p); % Coefficient matrix, where the rows correspond
                            % to all different combinations of N bosons on
                            % p lattices.
basis_size = length(coeffs); % Size of the state vector, where each element
                             % correspond to the Fock state described by
                             % the rows of the coeffs matrix.

[H_diag,H_offdiag] = hamiltonian(coeffs,1); % Non-periodic Hamiltonian
[H_diag_periodic,H_offdiag_periodic] = hamiltonian(coeffs,2); % Periodic Hamiltonian

% Vectors to store the energy values
energy = zeros(ratio_max/0.1,1);
energy_periodic = zeros(ratio_max/0.1,1);
% Vectors to store the occupancy expectation value for different sites
occupancy = zeros(ratio_max/0.1,p);
occupancy_periodic = zeros(ratio_max/0.1,p);

index=1; % index counter

for ratio=0.1:0.1:ratio_max, % Loop through different values of U/J
    U = ratio*J;
    % Diagonalize the Hamiltonians
    [v,e] = eigs(U*H_diag - J*H_offdiag,1,'sa');
    [v_periodic,e_periodic] = eigs(U*H_diag_periodic - J*H_offdiag_periodic,1,'sa');
    % Store energies
    energy(index) = e;
    energy_periodic(index) = e_periodic;
    % Store occupancies
    for i=1:p,
        for j=1:basis_size,
            occupancy(index,i) = occupancy(index,i) + abs(v(j))^2*coeffs(j,i);
            occupancy_periodic(index,i) = occupancy_periodic(index,i) + abs(v_periodic(j))^2*coeffs(j,i);
        end
    end
    index = index + 1;
end

% Plotting

figure
hold on
plot(0.1:0.1:ratio_max,energy,'-',0.1:0.1:ratio_max,energy_periodic,'-')
xlabel('U/J')
ylabel('E')
legend('1D lattice','periodic boundary conditions')
hold off

figure
hold on
for i=1:p,
    plot(0.1:0.1:ratio_max,occupancy(:,i),'color',colors{mod(i+4,5)+1},'DisplayName',['n_' num2str(i)])
end
xlabel('U/J')
ylabel('<n_i>')
legend(gca,'show')
hold off
figure
hold on
for i=1:p,
    plot(0.1:0.1:ratio_max,occupancy_periodic(:,i),'color',colors{mod(i+4,5)+1},'DisplayName',['n^{per}_' num2str(i)])
end
xlabel('U/J')
ylabel('<n_i>')
legend(gca,'show')
hold off
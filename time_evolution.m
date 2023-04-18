clear all
close all

% Time evolution 

J=1; % Site hopping strength
U=J; % On-site repulsion
N=3; % Number of bosons
p=3; % Number of lattice points

dt = 1e-3; % Timestep
Tmax = 20; % Total time

colors = {'r','g','b','m','k'};

coeffs = coefficients(N,p); % Coefficient matrix
basis_size = length(coeffs); % Size of the state vector

[H_diag,H_offdiag] = hamiltonian(coeffs,1); % Non-periodic Hamiltonian
[H_diag_periodic,H_offdiag_periodic] = hamiltonian(coeffs,2); % Periodic Hamiltonian

% Vectors for plotting
occupancy = zeros(Tmax/dt+1,p);
occupancy_periodic = zeros(Tmax/dt+1,p);

psi = zeros(length(coeffs),1); % Initial state vector
psi(end) = 1; % Initialize such that all probability amplitude is at final basis vector 
			  % (for the example N=3, p=3 this is the state where all bosons are at site 1)
psi_periodic = zeros(length(coeffs),1); % Initial state vector
psi_periodic(end) = 1;

index = 1;

for t=0:dt:Tmax, % Loop through times
    psi = expm(-1i*dt*(U*H_diag - J*H_offdiag))*psi; % The state after the time-step
    psi_periodic = expm(-1i*dt*(U*H_diag_periodic - J*H_offdiag_periodic))*psi_periodic; 
    for i=1:p,
        for j=1:basis_size,
            occupancy(index,i) = occupancy(index,i) + abs(psi(j))^2*coeffs(j,i);
            occupancy_periodic(index,i) = occupancy_periodic(index,i) + abs(psi_periodic(j))^2*coeffs(j,i);
        end
    end
    index = index + 1;
end

figure
hold on
for i=1:p,
    plot(0:dt:Tmax,occupancy(:,i),'color',colors{mod(i+4,5)+1},'DisplayName',['n_' num2str(i)])
end
xlabel('time')
ylabel('<n_i>')
legend(gca,'show')
hold off
figure
hold on
for i=1:p,
    plot(0:dt:Tmax,occupancy_periodic(:,i),'color',colors{mod(i+4,5)+1},'DisplayName',['n^{per}_' num2str(i)])
end
xlabel('time')
ylabel('<n_i>')
legend(gca,'show')
hold off
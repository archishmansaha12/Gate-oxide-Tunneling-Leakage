clc;
clear all

%% --------parameters--------

nm = 1e-9;  %define 1 nm length
q = 1.6e-19; %unit charge (1e)
m0 = 9.1e-31; %% mass of electron
h = 6.626e-34; 
hbar = h/2/pi; 
m_eff = 0.4*m0;

%Vg = 0.03;  %Gate voltage (INPUT)

Ea = 4.1*q; %Electron affinity of Si substrate = 4.1eV
Eg = 1.1*q; %Bandgap of Si = 1.1eV
X_c = 3.15*q; %Ec offset between Si and SiO2 = 3.25eV
phi_bi = X_c/q;
Eg_ox = 9*q;  %Bandgap of SiO2 = 9eV

t_oxide = [1.5*nm 3*nm];  %oxide thickness (INPUT)
N_sub = 5e18*1e6; %enter +ve for NMOS (INPUT)
rho_fix = 0;%2*1e12*1e4; %fixed oxide charge density (INPUT)
D_it = 0;%3*1e12*1e4; %interface trap density (INPUT)

eps_0 = 8.85e-12; %permittivity of free space
eps_sub = 11.6*eps_0; %permittivity of Si substrate
eps_ox = 3.9*eps_0; %permittivity of SiO2

ni = 1.5e10*1e6; %intrinsic carrier concentration
p_sub = ni;
n_sub = ni;

kbT = 25.6e-3*q;

Vg = 0.01:0.01:3;

%% ------------------Solve Poisson equation----------------------

N = length(Vg);

J_g1 = zeros(length(t_oxide),N);
J_g2 = zeros(length(t_oxide),N);

for k=1:length(t_oxide)
    
t_ox = t_oxide(k);

F_ox_arr_p = zeros(1,N);

%psi_s_arr = zeros(1,N);
%psi_s_poisson_arr = zeros(1,N);

for i=1:N

Vg1 = Vg(i)
%Vg1 = 1.3
[psi_s, F_ox] = poisson(m0,q,eps_sub,eps_ox,t_ox,N_sub,Eg,ni,Vg1,Ea);
%[Vg, F_ox] = coupled_poisson_quantum(psi_s,eps_sub,eps_ox,t_ox,25.6e-3,N_sub);

%Vg_arr(k,i) = Vg;

F_ox_arr_p(i) = F_ox;
%psi_s_arr(i) = psi_s;
%psi_s_poisson_arr(i) = psi_s_poisson;

end

% figure;
% plot(Vg_arr,psi_s_arr,'k--','linewidth',2);

%% ----------------modelling reflection term in tunneling------------------

m_ox = 0.61*m0;  %oxide effective mass
m_Si_per = 0.98*m0; %perpendicular effective mass in Silicon
m_Si_par = 0.19*m0; %parallel effective mass in Silicon
eta = 2;
J_g_wo_ref = zeros(1,length(F_ox_arr_p));
J_g = zeros(1,length(F_ox_arr_p));
f = zeros(1,length(F_ox_arr_p));

for i1 = 1:length(F_ox_arr_p)
    
    F_ox = F_ox_arr_p(i1);
    field = F_ox/1e8
    
%     V_ox = F_ox*t_ox;
%     
%     C = Vg_arr(k,i)/t_ox*N_inv*exp(20/(phi_bi)*(abs(V_ox)/phi_bi)^0.6*(1-abs(V_ox)/phi_bi));
%     J_g(i1) = q^2/(8*pi*h*eps_sub*phi_bi)*C*exp(-8*pi*sqrt(2*m_eff)*(q*phi_bi)^1.5/(3*h*q*F_ox)*(1-(1-abs(V_ox)/phi_bi)^1.5));
%     
    E_Si_per = 0.6*((3*pi*hbar*q*m_Si_per)^(2/3))/(2*m_Si_per)*(eps_ox*F_ox/eps_sub)^(2/3);
    f(i1) = 0.6*2*q/((3*pi*hbar*q*m_Si_per)^(1/3))*(eps_ox*F_ox/eps_sub)^(2/3);

    v_Si_per = sqrt(2*E_Si_per/m_Si_per); 

    E_Si_par_max = pi*hbar^2/(q*m_Si_par*eta)*eps_ox*F_ox;

    E_Si_par_arr = linspace(0,E_Si_par_max,10000); %range of E_Si_par
    dE = E_Si_par_arr(2) - E_Si_par_arr(1);
    T_WKB = zeros(1,length(E_Si_par_arr)); %WKB tunneling probability
    T_R = zeros(1,length(E_Si_par_arr)); %Reflection probability
    
    %%%%%% ------------------ without integral --------------------------
    
    E_Si_par = 0.5*E_Si_par_max;
    qphi_cat = X_c - (E_Si_per + E_Si_par);
    qphi_an = X_c - (E_Si_per + E_Si_par) - q*F_ox*t_ox;

    E_ox = zeros(1,2);
    E_ox(1) = qphi_cat;
    E_ox(2) = qphi_an;

    gamma_prime = (1-2*E_ox/Eg_ox);
    gamma = E_ox.*(1-E_ox/Eg_ox);
    v_ox = sqrt(2*gamma/m_ox)./gamma_prime;

    exponent = Eg_ox*(sqrt(2*m_ox)/(4*hbar*q*F_ox))*(2*gamma_prime.*sqrt(gamma)+sqrt(Eg_ox).*asin(gamma_prime));
    T_WKB_dir = exp(exponent(1)- exponent(2));
    v_Si_per1 = sqrt(2*(E_Si_per+q*F_ox*t_ox)/m_Si_per);
    T_R_dir = (4*v_Si_per*v_ox(1)/(v_Si_per^2 + v_ox(1)^2))*(4*v_Si_per1*v_ox(2)/(v_Si_per1^2 + v_ox(2)^2));
    
    J_g_wo_ref(i1) = eps_ox*F_ox*f(i1)*T_WKB_dir;
    J_g(i1) = eps_ox*F_ox*f(i1)*T_WKB_dir*T_R_dir;
    
%     %%%%%%% -------------- with integral -----------------------
% 
%     for i=1:length(E_Si_par_arr)
% 
%         E_Si_par = E_Si_par_arr(i);
%         qphi_cat = X_c - (E_Si_per + E_Si_par);
%         qphi_an = X_c - (E_Si_per + E_Si_par) - q*F_ox*t_ox;
% 
%         E_ox = zeros(1,2);
%         E_ox(1) = qphi_cat;
%         E_ox(2) = qphi_an;
% 
%         gamma_prime = (1-2*E_ox/Eg_ox);
%         gamma = E_ox.*(1-E_ox/Eg_ox);
%         v_ox = sqrt(2*gamma/m_ox)./gamma_prime;
% 
%         exponent = Eg_ox*(sqrt(2*m_ox)/(4*hbar*q*F_ox))*(2*gamma_prime.*sqrt(gamma)+sqrt(Eg_ox).*asin(gamma_prime));
%         T_WKB(i) = exp(exponent(1)- exponent(2));
%         v_Si_per1 = sqrt(2*(E_Si_per+q*F_ox*t_ox)/m_Si_per);
%         T_R(i) = (4*v_Si_per*v_ox(1)/(v_Si_per^2 + v_ox(1)^2))*(4*v_Si_per1*v_ox(2)/(v_Si_per1^2 + v_ox(2)^2));
% 
%     end
% 
%     J_wo_ref = sum(T_WKB)*dE;
%     J_w_ref = sum(T_WKB.*T_R)*dE;
%     J_g_wo_ref(i1) = eta*q*m_Si_par*f(i1)/(pi*hbar^2)*J_wo_ref;
%     J_g(i1) = eta*q*m_Si_par*f(i1)/(pi*hbar^2)*J_w_ref;

end

J_g1(k,:) = J_g;
J_g2(k,:) = J_g_wo_ref;

end
%% ---------------------------
figure;
semilogy(Vg,abs(J_g2(1,:))/1e4,'k:','linewidth',1.3);
hold on;
semilogy(Vg,abs(J_g1(1,:))/1e4,'k','linewidth',1.1);
hold on;
semilogy(Vg,abs(J_g2(2,:))/1e4,'k--','linewidth',2);
hold on;
semilogy(Vg,abs(J_g1(2,:))/1e4,'k','linewidth',2);
grid on;
xlim([0 3]);
xlabel('Gate Voltage V_{G} (in V)');
ylabel('Gate Current J_G (in A-cm^{-2})');
legend('t_{ox}=1.5nm without reflection term','t_{ox}=1.5nm with reflection term','t_{ox}=3.0nm without reflection term','t_{ox}=3.0nm with reflection term','Location','best');
title('Gate Current vs Gate voltage using Poisson solver');


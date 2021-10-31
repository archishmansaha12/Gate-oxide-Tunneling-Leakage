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
X_c_il = 3.5*q; %Ec offset between Si and SiO2 = 3.5eV
X_c_hk = 1.5*q; %Ec offset between Si and HfO2 = 1.5eV

phi_bi = X_c_il/q;

Eg_il = 9*q;   %Bandgap of SiO2 = 9eV
Eg_hk = 5.6*q; %Bandgap of HfO2 = 5.6eV

t_il = [0.54*nm 0.47*nm 0.41*nm 0.32*nm];  %IL thickness (INPUT)
t_hk = 2.2*nm;    % HK layer thickness
N_sub = 5e18*1e6; %enter +ve for NMOS (INPUT)
rho_fix = 0;%2*1e12*1e4; %fixed oxide charge density (INPUT)
D_it = 0;%3*1e12*1e4; %interface trap density (INPUT)

eps_0 = 8.85e-12; %permittivity of free space
eps_sub = 11.6*eps_0; %permittivity of Si substrate
eps_il = 3.9*eps_0; %permittivity of SiO2
eps_hk = 22*eps_0; %permittivity of HfO2

ni = 1.5e10*1e6; %intrinsic carrier concentration
p_sub = ni;
n_sub = ni;

kbT = 25.6e-3*q;

%Vg = 0.1:0.1:10;

%% ------------------Solve Schrodinger-Poisson equation----------------------


psi_s_arr = 0.01:0.01:2;
N = length(psi_s_arr);

J_g1 = zeros(length(t_il),N);
%J_g2 = zeros(length(t_il),N);
Vg_arr = zeros(length(t_il),N);

for k=1:length(t_il)
    
t_ox = t_il(k);

F_SiO2_arr = zeros(1,N);
F_HfO2_arr = zeros(1,N);

N_inv_arr = zeros(1,N);
%psi_s_poisson_arr = zeros(1,N);

for i=1:N

psi_s = psi_s_arr(i)
%Vg1 = 1.3
[Vg, N_inv, F_il, F_hk] = schrodinger_poisson_stack(m0,q,eps_sub,eps_il,eps_hk,t_ox,t_hk,N_sub,Eg,ni,psi_s,Ea);
%[Vg, F_ox] = coupled_poisson_quantum(psi_s,eps_sub,eps_il,t_ox,25.6e-3,N_sub);

Vg_arr(k,i) = Vg;
N_inv_arr(i) = N_inv;

F_SiO2_arr(i) = F_il;
F_HfO2_arr(i) = F_hk;
%psi_s_arr(i) = psi_s;
%psi_s_poisson_arr(i) = psi_s_poisson;

end

% figure;
% plot(Vg_arr,psi_s_arr,'k--','linewidth',2);

%% ----------------modelling reflection term in tunneling------------------

m_il = 0.61*m0;  % il effective mass
m_hk = 0.25*m0;  % high-k effective mass

m_Si_per = 0.98*m0; %perpendicular effective mass in Silicon
m_Si_par = 0.19*m0; %parallel effective mass in Silicon

eta = 2;
%J_g_wo_ref = zeros(1,length(F_SiO2_arr));
J_g = zeros(1,length(F_SiO2_arr));
f1 = zeros(1,length(F_SiO2_arr));
f2 = zeros(1,length(F_SiO2_arr));

for i1 = 1:length(F_SiO2_arr)
    
    F_SiO2 = F_SiO2_arr(i1);
    F_HfO2 = F_HfO2_arr(i1);
    field = F_SiO2/1e8;
    

    %%%%%% ------------------- for il layer------------------------------

    E_Si_per1 = 0.6*((3*pi*hbar*q*m_Si_per)^(2/3))/(2*m_Si_per)*(eps_il*F_SiO2/eps_sub)^(2/3);
    f1(i1) = 0.6*2*q/((3*pi*hbar*q*m_Si_per)^(1/3))*(eps_il*F_SiO2/eps_sub)^(2/3)
    
    v_Si_per = sqrt(2*E_Si_per1/m_Si_per);
    E_Si_par_max1 = pi*hbar^2/(q*m_Si_par*eta)*eps_il*F_SiO2;

    E_Si_par1 = 0.5*E_Si_par_max1
    qphi_cat1 = X_c_il - (E_Si_per1 + E_Si_par1);
    qphi_an1 = X_c_il - (E_Si_per1 + E_Si_par1) - q*F_SiO2*t_ox;

    E_ox1 = zeros(1,2);
    E_ox1(1) = qphi_cat1;
    E_ox1(2) = qphi_an1;

    gamma_prime1 = (1-2*E_ox1/Eg_il);
    gamma1 = E_ox1.*(1-E_ox1/Eg_il);
    v_ox1 = sqrt(2*gamma1/m_il)./gamma_prime1;

    exponent1 = Eg_il*(sqrt(2*m_il)/(4*hbar*q*F_SiO2))*(2*gamma_prime1.*sqrt(gamma1)+sqrt(Eg_il).*asin(gamma_prime1));
    T_WKB_dir1 = exp(exponent1(1)- exponent1(2));
    v_Si_per1 = sqrt(2*(E_Si_per1+q*F_SiO2*t_ox)/m_Si_per);
    T_R_dir1 = (4*v_Si_per*v_ox1(1)/(v_Si_per^2 + v_ox1(1)^2))*(4*v_Si_per1*v_ox1(2)/(v_Si_per1^2 + v_ox1(2)^2));
    
    
    %%%%%% ------------------- for hk layer------------------------------

    E_Si_per2 = 0.6*((3*pi*hbar*q*m_Si_per)^(2/3))/(2*m_Si_per)*(eps_hk*F_HfO2/eps_sub)^(2/3);
    f2(i1) = 0.6*2*q/((3*pi*hbar*q*m_Si_per)^(1/3))*(eps_hk*F_HfO2/eps_sub)^(2/3)

    v_Si_per = sqrt(2*E_Si_per2/m_Si_per);
    E_Si_par_max2 = pi*hbar^2/(q*m_Si_par*eta)*eps_hk*F_HfO2;

    E_Si_par2 = 0.5*E_Si_par_max2
    qphi_cat2 = X_c_hk - (E_Si_per2 + E_Si_par2);
    qphi_an2 = X_c_hk - (E_Si_per2 + E_Si_par2) - q*F_HfO2*t_hk;

    E_ox2 = zeros(1,2);
    E_ox2(1) = qphi_cat2;
    E_ox2(2) = qphi_an2;

    gamma_prime2 = (1-2*E_ox2/Eg_hk);
    gamma2 = E_ox2.*(1-E_ox2/Eg_hk);
    v_ox2 = sqrt(2*gamma2/m_hk)./gamma_prime2;

    exponent2 = Eg_hk*(sqrt(2*m_hk)/(4*hbar*q*F_HfO2))*(2*gamma_prime2.*sqrt(gamma2)+sqrt(Eg_hk).*asin(gamma_prime2));
    T_WKB_dir2 = exp(exponent2(1)- exponent2(2));
    v_Si_per2 = sqrt(2*(E_Si_per2+q*F_HfO2*t_hk)/m_Si_per);
    T_R_dir2 = (4*v_Si_per*v_ox2(1)/(v_Si_per^2 + v_ox2(1)^2))*(4*v_Si_per2*v_ox2(2)/(v_Si_per2^2 + v_ox2(2)^2));
    
    
    %J_g_wo_ref(i1) = eps_ox*F_ox*f(i1)*T_WKB_dir1;
    J_g(i1) = eps_il*F_SiO2*f1(i1)*T_WKB_dir1*T_R_dir1*T_WKB_dir2*T_R_dir2;
    
    
%     for i=1:length(E_Si_par_arr1)
% 
%         E_Si_par = E_Si_par_arr1(i);
%         qphi_cat = X_c - (E_Si_per1 + E_Si_par);
%         qphi_an = X_c - (E_Si_per1 + E_Si_par) - q*F_ox*t_ox;
% 
%         E_ox = zeros(1,2);
%         E_ox(1) = qphi_cat;
%         E_ox(2) = qphi_an;
% 
%         gamma_prime = (1-2*E_ox/Eg_ox);
%         gamma = E_ox.*(1-E_ox/Eg_ox);
%         v_ox = sqrt(2*gamma/m_il)./gamma_prime;
% 
%         exponent = Eg_ox*(sqrt(2*m_il)/(4*hbar*q*F_ox))*(2*gamma_prime.*sqrt(gamma)+sqrt(Eg_ox).*asin(gamma_prime));
%         T_WKB1(i) = exp(exponent(1)- exponent(2));
%         v_Si_per1 = sqrt(2*(E_Si_per1+q*F_ox*t_ox)/m_Si_per);
%         T_R1(i) = (4*v_Si_per1*v_ox(1)/(v_Si_per1^2 + v_ox(1)^2))*(4*v_Si_per1*v_ox(2)/(v_Si_per1^2 + v_ox(2)^2));
% 
%     end
% 
%     J_wo_ref = sum(T_WKB1)*dE1;
%     J_w_ref = sum(T_WKB1.*T_R1)*dE1;
%     J_g_wo_ref(i1) = eta*q*m_Si_par*f1(i1)/(pi*hbar^2)*J_wo_ref;
%     J_g(i1) = eta*q*m_Si_par*f1(i1)/(pi*hbar^2)*J_w_ref;

end

J_g1(k,:) = J_g;
%J_g2(k,:) = J_g_wo_ref;

end
%% ---------------------------
figure;
semilogy(Vg_arr(1,:),abs(J_g1(1,:))/1e4,'k','linewidth',1.1);
hold on;
semilogy(Vg_arr(2,:),abs(J_g1(2,:))/1e4,'k:','linewidth',2);
grid on;
semilogy(Vg_arr(3,:),abs(J_g1(3,:))/1e4,'k--','linewidth',1.5);
hold on;
semilogy(Vg_arr(4,:),abs(J_g1(4,:))/1e4,'k','linewidth',2);
grid on;
xlim([0 2]);
xlabel('Gate Voltage V_{G} (in V)');
ylabel('Gate Current J_G (in A-cm^{-2})');
legend('t_{eff}=0.97nm','t_{eff}=0.90nm','t_{eff}=0.84nm','t_{eff}=0.75nm','Location','best');
title('Gate Current vs Gate Voltage for HKMG stack using Schrodinger-Poisson solver (HfO_2-SiO_2 bilayer)');


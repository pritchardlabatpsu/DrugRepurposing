%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         FBS Model        %%%
%%% Scott Leighow - 04/13/21 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Objective: Compare steady-state drug-target complex values using
%%% numerical solutions of binding kinetic equations

%% Set up

clear
close all

%% Initialize

% Species concentrations (total) [uM]
% using numbers for imatinib (see SerumBindingAnalytical.m for rationale)

% In vivo drug concentration (ASH TKI 2015 poster)
D_vivo = 3.4;   % uM

% In vivo target concentration
BCRABL_actin_ratio = 0.33;  % median value reported in Barnes et al Oncogene 2005
actin_conc = 500;    % [uM] intracellular actin conc (Molecular Cell Biology. 4th edition Section 18.1)
BCRABL_conc = BCRABL_actin_ratio*actin_conc;

blood_vol = 10;     % L (5 L peripheral blood + 5 L marrow)
WBC_count = 98e9;   % [/L] median WBC count at presentation of CML (Qin et al Medicine 2016)
WBC_totcount = WBC_count*blood_vol;   
WBC_vol = 130e-15;  % [L] (from um^3)
WBC_totvol = WBC_totcount*WBC_vol;  % [L]

BCRABL_moles = WBC_totvol*BCRABL_conc;  % micro-moles of BCR-ABL molecules in blood
T_vivo = BCRABL_moles/blood_vol; % uM

% In vivo HSA concentration
HSA_gmL = 6.9/300; % [g/mL]
HSA_uM = HSA_gmL*1000/66500*1e6; % uM
H_vivo = HSA_uM;   % in vivo human serum albumin (HSA) concentration

% In vitro target concentration
BaF3_vol = WBC_vol;
BaF3_seedingcount = 3e3;
BaF3_totvol = BaF3_vol*BaF3_seedingcount;
BaF3_BCRABLconc = 1000; % assume high BCRABL expression from CMV promoter
BCRABL_vitromoles = BaF3_totvol*BaF3_BCRABLconc; % uM

BaF3_seedingvol = 100/1e6; % [L]
T_vitro = BCRABL_vitromoles/BaF3_seedingvol; % [uM]

H_vitro = H_vivo;   % in vitro HSA concentration

B_mgmL = 45;        % BSA in vitro conc [mg/mL]
B_uM = B_mgmL/66000*1e6; % [uM]
B_vitro = B_uM;

% Affinities [uM]
K_T_vec = 0.001; %logspace(-3,-1,3); % drug-target binding affinity [uM]
K_B_vec = .1; %logspace(-3,-1,3); % drug-BSA binding affinity [uM]
K_H_vec = .0001; % logspace(-3,-1,3); % drug-HSA binding affinity [uM]

% Model parameters
eps = 1e-3;     % Run simulation until within eps of steady state

% Storage variables
fracbounds = NaN(length(K_T_vec),length(K_B_vec),length(K_H_vec));
fracbound_vivo = fracbounds;
fracbound_vitro_freeCave = fracbounds;
fracbound_vitro_effCave = fracbounds;

% Initialize
j = 1;
k = 1;
m = 1;

%% Simulations

for K_T = K_T_vec
    for K_B = K_B_vec
        for K_H = K_H_vec
            
            % Define rate constants based on affinities
            % Assume typical on-rate 1e5 /M/s (discussed in Cell Signaling
            % Lim et al 2015)
            kon_T = 1e-1;       % drug-target on rate [/uM-s]
            koff_T = K_T*kon_T; % drug-target off rate [/s]
            
            kon_B = 1e-1;       % drug-BSA on rate [/uM-s]
            koff_B = K_B*kon_B; % drug-BSA off rate [/s]
            
            kon_H = 1e-1;       % drug-HSA on rate [/uM-s]
            koff_H = K_H*kon_H; % drug-HSA on rate [/s]
            
            %% In vivo (target and HSA present)
            
            % State: x = [D T H DT DH]
            dxdt_vivo = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(4)-kon_H*x(1)*x(3)+koff_H*x(5)
                                -kon_T*x(1)*x(2)+koff_T*x(4)
                                -kon_H*x(1)*x(3)+koff_H*x(5)
                                kon_T*x(1)*x(2)-koff_T*x(4)
                                kon_H*x(1)*x(3)-koff_H*x(5)];
            
            % Initial conditions
            x0 = [D_vivo T_vivo H_vivo 0 0];
            % Solve ODE
            tend = 1;
            ss = false;
            
            % Update simulation time if steady-state not reached
            while ~ss
                
                tend = tend*10;
                
                [tout,xout] = ode45(dxdt_vivo,[0 tend],x0);
                
                % Check if steady state reached
                dx_fin = xout(end,:)-xout(end-1,:);
                dt_fin = tout(end)-tout(end-1);
                dxdt_fin = abs(dx_fin/dt_fin); % estimate rate of state change at end of simulation
                if sum(dxdt_fin > eps) == 0 % check if any variable is still changing
                    ss = true;
                end
                
            end
            
            % Visualize results
%             figure
%             plot(tout,xout)
%             legend('D','T','H','DT','DH')
            
            % Save results of simulation
            ss_vivo = xout(end,:);
            
            Cave = sum(ss_vivo([1 5])); % in vivo plasma drug concentration
            Cave_free = ss_vivo(1);
            
            %% In vitro simulations (range of drug concentrations)
            
            D_vitro_vec = logspace(0,4,20); % concentration of drug [uM]
            
            DT_noHSA = NaN(size(D_vitro_vec));
            DT_withHSA = NaN(size(D_vitro_vec));
            
            i = 1;

            for D_vitro = D_vitro_vec
                
                %% In vitro -HSA (target and BSA present)
                % state: x = [D T B DT DB]
                dxdt_vitronoHSA = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(4)-kon_B*x(1)*x(3)+koff_B*x(5)
                                          -kon_T*x(1)*x(2)+koff_T*x(4)
                                          -kon_B*x(1)*x(3)+koff_B*x(5)
                                          kon_T*x(1)*x(2)-koff_T*x(4)
                                          kon_B*x(1)*x(3)-koff_B*x(5)];

                % Initial conditions
                x0 = [D_vitro T_vitro B_vitro 0 0];

                % Solve ODE
                tend = 1;
                ss = false;

                % Update simulation time if steady-state not reached
                while ~ss

                    tend = tend*10;

                    [tout,xout] = ode45(dxdt_vitronoHSA,[0 tend],x0);

                    % Check if steady state reached
                    dx_fin = xout(end,:)-xout(end-1,:);
                    dt_fin = tout(end)-tout(end-1);
                    dxdt_fin = abs(dx_fin/dt_fin); % estimate rate of state change at end of simulation
                    if sum(dxdt_fin > eps) == 0 % check if any variable is still changing
                        ss = true;
                    end

                end

                % Visualize results
%                 figure
%                 plot(tout,xout)
%                 legend('D','T','B','DT','DB')
                
                % Save results of simulation
                DT_noHSA(i) = xout(end,4); % steady-state DT value

                %% In vitro +HSA (target, BSA, and HSA present)
                % state: x = [D T B H DT DB DH]
                dxdt_vitrowithHSA = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(5)-kon_B*x(1)*x(3)+koff_B*x(6)-kon_H*x(1)*x(4)+koff_H*x(7)
                                            -kon_T*x(1)*x(2)+koff_T*x(5)
                                            -kon_B*x(1)*x(3)+koff_B*x(6)
                                            -kon_H*x(1)*x(4)+koff_H*x(7)
                                            kon_T*x(1)*x(2)-koff_T*x(5)
                                            kon_B*x(1)*x(3)-koff_B*x(6)
                                            kon_H*x(1)*x(4)-koff_H*x(7)];

                % Initial conditions
                x0 = [D_vitro T_vitro B_vitro H_vitro 0 0 0];

                % Solve ODE
                tend = 1;
                ss = false;

                % Update simulation time if steady-state not reached
                while ~ss

                    tend = tend*10;

                    [tout,xout] = ode45(dxdt_vitrowithHSA,[0 tend],x0);

                    % Check if steady state reached
                    dx_fin = xout(end,:)-xout(end-1,:);
                    dt_fin = tout(end)-tout(end-1);
                    dxdt_fin = abs(dx_fin/dt_fin); % estimate rate of state change at end of simulation
                    if sum(dxdt_fin > eps) == 0 % check if any variable is still changing
                        ss = true;
                    end

                end
                
                % Visualize results
%                 figure
%                 plot(tout,xout)
%                 legend('D','T','B','H','DT','DB','DH')
                
                % Save results of simulation
                DT_withHSA(i) = xout(end,5); % steady-state DT value
                   
                i = i+1;
                disp(i)
                
            end
            
            % Plot serum shift (fraction of target bound to drug)
            figure
            semilogx(D_vitro_vec,[DT_noHSA; DT_withHSA]/T_vitro)
            legend('-HSA','+HSA')
            
            %% Calculate Effective Cave
            
            % Estimate IC50 for dose-response curves with/without HSA
            y_noHSA = DT_noHSA/T_vitro;
            p0 = [0 2 400 1];  %Set initial guess for parameter vector
            p_noHSA = fminsearch(@(p)fun(p,D_vitro_vec',y_noHSA'),p0);  % Call optimizer for fitting
%             Ysim = p_noHSA(4)+(p_noHSA(1)-p_noHSA(4))./(1+(D_vitro_vec/p_noHSA(3)).^p_noHSA(2)); % evaluate final logistic function at measurement points
%             semilogx(D_vitro_vec,y_noHSA,D_vitro_vec,Ysim) %plot measured and fitted values
            IC50_noHSA = p_noHSA(3);
            
            y_withHSA = DT_withHSA/T_vitro;
            p0 = [0 2 400 1];  %Set initial guess for parameter vector
            p_withHSA = fminsearch(@(p)fun(p,D_vitro_vec',y_withHSA'),p0);  % Call optimizer for fitting
%             Ysim = p_noHSA(4)+(p_noHSA(1)-p_noHSA(4))./(1+(D_vitro_vec/p_noHSA(3)).^p_noHSA(2)); % evaluate final logistic function at measurement points
%             semilogx(D_vitro_vec,y_noHSA,D_vitro_vec,Ysim) %plot measured and fitted values
            IC50_withHSA = p_withHSA(3);
            
            shift = IC50_withHSA/IC50_noHSA;
            
            % Calculate effective Cave
            
            Cave_eff = Cave/shift;
            
            %% Calculate in vitro DT at free Cave
            
            % state: x = [D T B DT DB]
            dxdt_vitronoHSA = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(4)-kon_B*x(1)*x(3)+koff_B*x(5)
                                        -kon_T*x(1)*x(2)+koff_T*x(4)
                                        -kon_B*x(1)*x(3)+koff_B*x(5)
                                        kon_T*x(1)*x(2)-koff_T*x(4)
                                        kon_B*x(1)*x(3)-koff_B*x(5)];

            % Initial conditions
            x0 = [Cave_free T_vitro B_vitro 0 0];

            % Solve ODE
            tend = 1;
            ss = false;

            % Update simulation time if steady-state not reached
            while ~ss

                tend = tend*10;

                [tout,xout] = ode45(dxdt_vitronoHSA,[0 tend],x0);

                % Check if steady state reached
                dx_fin = xout(end,:)-xout(end-1,:);
                dt_fin = tout(end)-tout(end-1);
                dxdt_fin = abs(dx_fin/dt_fin); % estimate rate of state change at end of simulation
                if sum(dxdt_fin > eps) == 0 % check if any variable is still changing
                    ss = true;
                end

            end

            % Save results of simulation
            DT_freeCave = xout(end,4); % steady-state DT value at effective Cave

            
            %% Calculate in vitro DT at effective Cave
            
            % state: x = [D T B DT DB]
            dxdt_vitronoHSA = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(4)-kon_B*x(1)*x(3)+koff_B*x(5)
                                        -kon_T*x(1)*x(2)+koff_T*x(4)
                                        -kon_B*x(1)*x(3)+koff_B*x(5)
                                        kon_T*x(1)*x(2)-koff_T*x(4)
                                        kon_B*x(1)*x(3)-koff_B*x(5)];

            % Initial conditions
            x0 = [Cave_eff T_vitro B_vitro 0 0];

            % Solve ODE
            tend = 1;
            ss = false;

            % Update simulation time if steady-state not reached
            while ~ss

                tend = tend*10;

                [tout,xout] = ode45(dxdt_vitronoHSA,[0 tend],x0);

                % Check if steady state reached
                dx_fin = xout(end,:)-xout(end-1,:);
                dt_fin = tout(end)-tout(end-1);
                dxdt_fin = abs(dx_fin/dt_fin); % estimate rate of state change at end of simulation
                if sum(dxdt_fin > eps) == 0 % check if any variable is still changing
                    ss = true;
                end

            end

            % Save results of simulation
            DT_effCave = xout(end,4); % steady-state DT value at effective Cave
            
            %% Compare in vivo and in vitro results
            fracbound_vivo(j,k,m) = ss_vivo(4)/T_vivo; % fraction of target bound in vivo
            fracbound_vitro_freeCave(j,k,m) = DT_freeCave/T_vitro;
            fracbound_vitro_effCave(j,k,m) = DT_effCave/T_vitro; % fraction of target bound in vitro at effective Cave
            
            m = m+1;
        end
        k = k+1;
    end    
    j = j+1;
end

%% Test
% 
% % Simplest model
% % State: [D T DT]
% 
% kon_T = .1;
% koff_T = .001;
% 
% D_vec = logspace(0,4);
% 
% DT_ss = NaN(size(D_vec));
% i = 1;
% 
% for D = D_vec
% 
%     x0 = [D 1001 0];
% 
%     dxdt = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(3)
%                    -kon_T*x(1)*x(2)+koff_T*x(3)
%                    kon_T*x(1)*x(2)-koff_T*x(3)];
% 
%     [tout,xout] = ode45(dxdt,[0 10],x0);
%     
%     DT_ss(i) = xout(end,3);
%     i = i+1;
% 
% end
% 
% semilogx(D_vec,DT_ss)
% % 
% %% Test
% 
% % Model with competition
% % State: [D T B DT DB]
% 
% kon_T = .1;
% koff_T = .001;
% 
% kon_B = .1;
% koff_B = .1;
% 
% D_vec = logspace(0,4);
% 
% DT_ss = NaN(size(D_vec));
% i = 1;
% 
% for D = D_vec
% 
%     x0 = [D T_vitro B_vitro 0 0];
% 
%     dxdt = @(t,x) [-kon_T*x(1)*x(2)+koff_T*x(4) - kon_B*x(1)*x(3)+koff_B*x(5)
%                    -kon_T*x(1)*x(2)+koff_T*x(4)
%                    -kon_B*x(1)*x(3)+koff_B*x(5)
%                    kon_T*x(1)*x(2)-koff_T*x(4)
%                    kon_B*x(1)*x(3)-koff_B*x(5)];
% 
%     [tout,xout] = ode45(dxdt,[0 10],x0);
%     
% %     figure
% %     plot(tout,xout)
% %     legend('D','T','B','DT','DB')
%     
%     DT_ss(i) = xout(end,4);
%     i = i+1;
% 
% end
% 
% semilogx(D_vec,DT_ss)

%%

function obj = fun(p,X,Y)
  Ysim = p(4)+(p(1)-p(4))./(1+(X/p(3)).^p(2));  % evaluate logistic function at measurement points
  obj = (Y-Ysim).'*(Y-Ysim);  % calculate sum of squared differences between measured and fitted values
end





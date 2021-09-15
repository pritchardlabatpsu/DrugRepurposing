
%%% Analytical model of 3-way competitive binding
%%% Scott Leighow - 04/15/21

%% Set up

clear
close all

%% Parameters

% Affinities [uM]
K_T = .01;                     % drug-target binding affinity
K_B_vec = logspace(0,4)*K_T;   % drug-BSA binding affinity
K_H_vec = logspace(0,4)*K_T;   % drug-HSA binding affinity

% Initial species concentrations [uM]
T_tot = .005;  % target concentration
B_tot = 70; %  % BSA concentration
H_tot = 350; % % HSA concentration

%% Loop through drug concentrations
D_vec = logspace(-3,3,250);

IC50_vivo = NaN(length(K_B_vec),length(K_H_vec));
IC50_vitro = IC50_vivo;
IC50_vitroHSA = IC50_vivo;

for i = 1:length(K_B_vec)
    for j = 1:length(K_H_vec)
        
        K_B = K_B_vec(i);
        K_H = K_H_vec(j);
        
        % Identify steady-state DT/T_tot values for range of drug conc
        curve_vivo = SSFinder(K_T,K_B,K_H,T_tot,0,H_tot,D_vec);         % in vivo case - no BSA
        curve_vitro = SSFinder(K_T,K_B,K_H,T_tot,B_tot,0,D_vec);        % in vitro case - no HSA
        curve_vitroHSA = SSFinder(K_T,K_B,K_H,T_tot,B_tot,H_tot,D_vec); % in vitro case with HSA
        
        % Visualize results
%         semilogx(D_vec,curve_vivo)
%         hold on
%         semilogx(D_vec,curve_vitro)
%         semilogx(D_vec,curve_vitroHSA)
%         legend('in vivo','in vitro -HSA','in vitro +HSA')
        
        % Identify IC50s - find value in D_vec that achieves closest to 50%
        % max inhibition
        [~,IC50vivo_idx] = min(abs(curve_vivo-0.5));
        IC50_vivo(i,j) = D_vec(IC50vivo_idx);
        
        [~,IC50vitro_idx] = min(abs(curve_vitro-0.5));
        IC50_vitro(i,j) = D_vec(IC50vitro_idx);
        
        [~,IC50vitroHSA_idx] = min(abs(curve_vitroHSA-0.5));
        IC50_vitroHSA(i,j) = D_vec(IC50vitroHSA_idx);
        
        % Verify that range of D_vec is sufficient to capture IC50
        check = [IC50vivo_idx IC50vitro_idx IC50vitroHSA_idx];
        if sum(check==1)>0 || sum(check==length(D_vec))>0
            disp([i j]);
        end

    end
end

% Export results
writematrix(IC50_vivo,'IC50vivo_041921.csv')
writematrix(IC50_vitro,'IC50vitro_041921.csv')
writematrix(IC50_vitroHSA,'IC50vitroHSA_041921.csv')

writematrix(K_B_vec/K_T,'K_BFoldIncreaseOverK_T.txt')
writematrix(K_H_vec/K_T,'K_HFoldIncreaseOverK_T.txt')

%% Example Curves

% Pick moderate off-target affinities
K_B = 1e3*K_T;
K_H = 1e3*K_T;

curve_vivo = SSFinder(K_T,K_B,K_H,T_tot,0,H_tot,D_vec);         % in vivo case - no BSA
curve_vitro = SSFinder(K_T,K_B,K_H,T_tot,B_tot,0,D_vec);        % in vitro case - no HSA
curve_vitroHSA = SSFinder(K_T,K_B,K_H,T_tot,B_tot,H_tot,D_vec); % in vitro case with HSA

figure
semilogx(D_vec,curve_vivo)
hold on
semilogx(D_vec,curve_vitro)
semilogx(D_vec,curve_vitroHSA)
legend('in vivo','in vitro -HSA','in vitro +HSA')

% Export results
example_curve = [D_vec' curve_vivo' curve_vitro' curve_vitroHSA'];
writematrix(example_curve,'ExampleCurve_041921.csv')

%% Function SSfinder
%  identifies the steady-state fractional occupancy of the drug target for
%  a range of drug concentrations D_vec

function frcbnd = SSFinder(K_T,K_B,K_H,T_tot,B_tot,H_tot,D_vec)

frcbnd = NaN(size(D_vec));
i = 1;

for D_tot = D_vec

    %%% Solution (see notes in OneNote for derivation)

    % Define binomials
    bin_T = [1 K_T];
    bin_B = [1 K_B];
    bin_H = [1 K_H];
    bin_Dtot = [1 -D_tot];

    % Multiply binomials
    poly1 = conv(conv(conv(bin_T,bin_B),bin_H),bin_Dtot);
    poly2 = T_tot*conv(conv(bin_B,bin_H),[1 0]);
    poly3 = B_tot*conv(conv(bin_T,bin_H),[1 0]);
    poly4 = H_tot*conv(conv(bin_T,bin_B),[1 0]);

    % Pad polynomials and add
    poly2_0 = [0 poly2];
    poly3_0 = [0 poly3];
    poly4_0 = [0 poly4];

    poly = poly1+poly2_0+poly3_0+poly4_0;

    % Find roots of polynomial

    rts = roots(poly);

    % Find stable solution

    if sum(rts>0)==1
        D_ss = rts(rts>0);
    else
        D_ss = NaN;
    end

    % Find steady-state values for other species

    DT_ss = D_ss*T_tot/(D_ss+K_T);
    T_ss = T_tot - DT_ss;

    DB_ss = D_ss*B_tot/(D_ss+K_B);
    B_ss = B_tot - DB_ss;

    DH_ss = D_ss*H_tot/(D_ss+K_H);
    H_ss = H_tot - DH_ss;

    %% Analyze solutions

    % Calculate fraction of target bound to drug
    frcbnd(i) = DT_ss/T_tot;
    i = i+1;

end


end
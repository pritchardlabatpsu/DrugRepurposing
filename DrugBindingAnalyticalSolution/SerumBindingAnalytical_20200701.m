clear
close all

%% Parameters

%%% Concentrations %%%

% Drug concentration (Rivera et al ASH 2014 poster)
D_tot = 3.4;      % uM (imatinib)

% Target concentration
BCRABL_actin_ratio = 0.33;  % median value reported in Barnes et al Oncogene 2005
actin_conc = 500;    % [uM] intracellular actin conc (Molecular Cell Biology. 4th edition Section 18.1)
BCRABL_conc = BCRABL_actin_ratio*actin_conc;

blood_vol = 10;     % L (5 L peripheral blood + 5 L marrow)
WBC_count = 98e9;   % [/L] median WBC count at presentation of CML (Qin et al Medicine 2016)
WBC_totcount = WBC_count*blood_vol;   
WBC_vol = 130e-15;  % [L] (from um^3)
WBC_totvol = WBC_totcount*WBC_vol;  % [L]

BCRABL_moles = WBC_totvol*BCRABL_conc;  % micro-moles of BCR-ABL molecules in blood
BCRABL_blood_conc = BCRABL_moles/blood_vol; % uM

T_tot_c = BCRABL_blood_conc;     % uM

% AAG concentration (Gambacorti-Passerini et al CCR 2003)
S_conc = 1;   % g/L
S_MW = 23512; % g/mol
S_tot = S_conc/S_MW*1e6; % uM

% justin estimates
HSA = 340; % [uM]
AAG = 40; % [uM]

%%% Binding affinities %%%

Kd_T_vec = logspace(-5,1,10);  % [uM] ~0.013 uM (Lin et al PNAS 2013)
Kd_S_vec = logspace(-3,2,10);  % [uM] ~5 uM (Gambacorti-Passerini)

%% Find steady-state solutions using analytical solution

DnS_rat = NaN(length(Kd_T_vec),length(Kd_S_vec));
DnS_wc_mat = DnS_rat;
DnS_woc_mat = DnS_rat;

Dss_a = DnS_rat;
Dss_n = DnS_rat;

out = NaN(length(Kd_T_vec)*length(Kd_S_vec),5);
out_row = 1;

for i = 1:length(Kd_T_vec)
    
    Kd_T = Kd_T_vec(i);
    
    for j = 1:length(Kd_S_vec)
        
        Kd_S = Kd_S_vec(j);
        
        %%% With competition %%%
        T_tot = T_tot_c;
        C_1 = 1;
        C_2 = S_tot + T_tot + Kd_T + Kd_S - D_tot;
        C_3 = Kd_T*S_tot + Kd_S*T_tot - D_tot*(Kd_T+Kd_S) + Kd_S*Kd_T;
        C_4 = -Kd_S*Kd_T*D_tot;

        rts_wc = roots([C_1 C_2 C_3 C_4]);
        
        rts_pssbl_wc = rts_wc>0 & imag(rts_wc)==0;
        if sum(rts_pssbl_wc)~=1
            continue
        end
        
        Dss_wc = rts_wc(rts_pssbl_wc);    % steady-state unbound drug
        DS_wc = D_tot - T_tot*Dss_wc/(Dss_wc+Kd_T)-Dss_wc;   % steady-state drug bound to serum
        DnS_wc = (D_tot - DS_wc)/D_tot;     % fraction drug not bound to serum
        
        %%% With competition %%%
        T_tot = 0;
        C_1 = 1;
        C_2 = S_tot + T_tot + Kd_T + Kd_S - D_tot;
        C_3 = Kd_T*S_tot + Kd_S*T_tot - D_tot*(Kd_T+Kd_S) + Kd_S*Kd_T;
        C_4 = -Kd_S*Kd_T*D_tot;

        rts_woc = roots([C_1 C_2 C_3 C_4]);
        
        rts_pssbl_woc = rts_woc>0 & imag(rts_woc)==0;
        if sum(rts_pssbl_woc)~=1
            continue
        end
        
        Dss_woc = rts_woc(rts_pssbl_woc);    % steady-state unbound drug
        DS_woc = D_tot - T_tot*Dss_woc/(Dss_woc+Kd_T)-Dss_woc;   % steady-state drug bound to serum
        DnS_woc = (D_tot - DS_woc)/D_tot;     % fraction drug not bound to serum
        
        Dss_wc_mat(i,j) = Dss_wc/D_tot;
        Dss_woc_mat(i,j) = Dss_woc/D_tot;
        
        out(out_row,:) = [Kd_T Kd_S DnS_wc DnS_woc DnS_wc/DnS_woc];
        out_row = out_row+1;
        
    end
end

%% Save output

csvwrite('DrugNotBoundToSerum_20200701.csv',out)

figure
h1 = heatmap(Dss_wc_mat);
h1.XData = log10(Kd_S_vec);
h1.YData = log10(Kd_T_vec);

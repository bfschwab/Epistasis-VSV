function [defect, defect_struct_cutoff] = VSV_get_defect(types)
% Types should be the number of defects available for change
% As of 4/5/19 (B.Schwab) there are 13 defects
% As of 4/19/19 (B.Schwab) there are 11 or 13 defects. condm and kdp.m may
%   need to be removed
% As of 1/12/23 (B.Schwab) the final paper uses 11 defects

%% Setup factor changes
if types ~= 13 && types ~= 11
    error('n_defects doesnt match current pool')
end

defect = cell(1,types);

% all multiplications after vector are post initial settings modification
eta_n_factor = [0 .2 .4 .6 .8 1];
eta_p_factor = [0 .2 .4 .6 .8 1];
eta_m_factor = [0 .2 .4 .6 .8 1];
eta_g_factor = [0 .2 .4 .6 .8 1];
eta_l_factor = [0 .2 .4 .6 .8 1];

kdp_n_factor = [5 4 3 2 1 0]*480/5 + 20; % B.S 4/19/19 expand range (5x) to see bigger effect
    %
kdp_m_factor = [5 4 3 2 1 0]*9/5 + 1;
kd_m_factor = [5 4 3 2 1 0]*9/5 + 1;
kd_nc_factor = [5 4 3 2 1 0]*96.5/5 + 3.5; % B.S 4/19/19 (x10) expand range to see bigger effect

% ke_pol_factor = [0 .2 .4 .6 .8 1]; restrict range, ke_pol taking to much
% fitness
ke_pol_factor = [0 1 2 3 4 5]*.7/5 + .3;
sprom_factor = [0 .2 .4 .6 .8 1]*.85;
condm_factor = [5 4 3 2 1 0]*9/5 + 1; % B.S 4/19/19 avoid w_poor > w_poor_wild
k_int_factor = [0 .2 .4 .6 .8 1];

%% Fill defect with mutation information

% Eta
defect{1} = struct();
defect{1}.par = ["eta", "n"];
defect{1}.factor = eta_n_factor;
defect{1}.light = linspace(defect{1}.factor(4),defect{1}.factor(6),6);

defect{2} = struct();
defect{2}.par = ["eta", "p"];
defect{2}.factor = eta_p_factor;
defect{2}.light = linspace(defect{2}.factor(4),defect{2}.factor(6),6);

defect{3} = struct();
defect{3}.par = ["eta", "m"];
defect{3}.factor = eta_m_factor;
defect{3}.light = linspace(defect{3}.factor(4),defect{3}.factor(6),6);

defect{4} = struct();
defect{4}.par = ["eta", "g"];
defect{4}.factor = eta_g_factor;
defect{4}.light = linspace(defect{4}.factor(4),defect{4}.factor(6),6);

defect{5} = struct();
defect{5}.par = ["eta", "l"];
defect{5}.factor = eta_l_factor;
defect{5}.light = linspace(defect{5}.factor(4),defect{5}.factor(6),6);


% Degredation constants

defect{6} = struct();
defect{6}.par = ["kdp", "n"];
defect{6}.factor = kdp_n_factor;
defect{6}.light = linspace(defect{6}.factor(4),defect{6}.factor(6),6);

defect{7} = struct();
defect{7}.par = ["kdp", "m"];
defect{7}.factor = kdp_m_factor;
defect{7}.light = linspace(defect{7}.factor(4),defect{7}.factor(6),6);

defect{8} = struct();
defect{8}.par = 'kd_m';
defect{8}.factor = kd_m_factor;
defect{8}.light = linspace(defect{8}.factor(4),defect{8}.factor(6),6);

defect{9} = struct();
defect{9}.par = 'kd_nc';
defect{9}.factor = kd_nc_factor;
defect{9}.light = linspace(defect{9}.factor(4),defect{9}.factor(6),6);

% Other

defect{10} = struct();
defect{10}.par = 'ke_pol';
defect{10}.factor = ke_pol_factor;
defect{10}.light = linspace(defect{10}.factor(4),defect{10}.factor(6),6);

defect{11} = struct();
defect{11}.par = 'sprom';
defect{11}.factor = sprom_factor;
defect{11}.light = linspace(defect{11}.factor(4),defect{11}.factor(6),6);

defect{12} = struct();
defect{12}.par = 'condm';
defect{12}.factor = condm_factor;
defect{12}.light = linspace(defect{12}.factor(4),defect{12}.factor(6),6);

defect{13} = struct();
defect{13}.par = 'k_int';
defect{13}.factor = k_int_factor;
defect{13}.light = linspace(defect{13}.factor(4),defect{13}.factor(6),6);


if types == 11

    defect(12) = []; %remove condm as a changeable defect
    defect(7) = []; %remove kdp.m as a changeable defect
    defect_struct_cutoff = 6; %indicies 6 and below are structure variables

elseif types == 13

    defect_struct_cutoff = 7; %indicies 7 and below are structure variables
end

end

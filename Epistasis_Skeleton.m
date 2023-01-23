function Epistasis_Skeleton(n_reps, n_defects, severity, ...
    cluster_num, process_num, light)
% n_reps is the number of repititions at each number of defects tested
% n_defects is the max number of defects that will be tested
%   n_defects should be in [1,defect_types] (see VSV_get_defect)
% tlenght is the length of each simulation.
%   As of 4/5/19 (B.Schwab) tlength = 25*3600; %25 hr
% dt is the time step for RK4 in the VSV toolbox
%   As of 4/5/19 (B.Schwab) dt = 15; %sec
% severity is the severity of defects desired
%   severity should be in [1,5], (see VSV_get_defect for intervals)
% defect_types is the number of defects in the current pool
% process_num and cluster_num were used in the distributed computing to
% uniquely seed the random number generator and save the results,
% any integers can be used for these values
% light is a boolean where false uses the original five parameter defect
% ranges, and true splits the least severe two ranges into five new ones
% to produce less deleterious mutations 

if ischar(n_reps)
    n_reps = str2double(n_reps);
end

if ischar(n_defects)
    n_defects = str2double(n_defects);
end

if ischar(severity)
    severity = str2double(severity);
end

if ischar(cluster_num)
    cluster_num = str2double(cluster_num);
end

if ischar(process_num)
    process_num = str2double(process_num);
end

%% setup pars with default parameters
format long g
defect_types = 11;
tlength = 25*3600; %25 hr
dt = 15; %seconds per step

defpars = Default_VSV_toolbox_RK4();
[defect, defect_struct_cutoff] = VSV_get_defect(defect_types);

seed1 = str2double(strcat(num2str(cluster_num),num2str(process_num)));

%% Not plotting from this toolbox
NO_PLOTS = 1;



%% Set severity of defects

% severity should be in [1, 2, 3, 4, 5]
% corresponding to [(.8,1) (.6,.8) (.4,.6) (.2,.4) (0,.2)]
inter1 = 6 - severity; % min effect (5,6), max effect (1,2)
inter2 = 7 - severity;

% setup random numbers (avoid same strings of numbers between trials)
index_store = cell(n_defects,n_reps);
fac_store = cell(n_defects,n_reps);

rng(seed1,'twister')
for q = 1:n_defects
    for r = 1:n_reps
        index_store{q,r} = randperm(defect_types,q);
    end
end

rng('shuffle','twister')
for q = 1:n_defects
    for r = 1:n_reps
        fac_store{q,r} = rand(1,q);
    end
end

% initialize output
store = cell(1,n_defects);


%% Run

for j = 1:n_defects

    store{j} = struct();

    w_rich = NaN(n_reps,1);
    w_poor = NaN(n_reps,1);
    w_info = strings(n_reps,j);


    for i = 1:n_reps

        input = defpars;

%         index = randperm(defect_types,j);
%         fac = rand(1,j);
        index = index_store{j,i};
        fac = fac_store{j,i};

        scale = NaN(1,j);

        for m = 1:j

            if ~light
                scale(m) = fac(m)*(defect{index(m)}.factor(inter2)-defect{index(m)}.factor(inter1)) + ...
                    defect{index(m)}.factor(inter1);
            elseif light
                scale(m) = fac(m)*(defect{index(m)}.light(inter2)-defect{index(m)}.light(inter1)) + ...
                    defect{index(m)}.light(inter1);
            else
                error('Input light needs to be a logical')
            end


            if index(m) <= defect_struct_cutoff
                new_par = input.(defect{index(m)}.par(1)).(defect{index(m)}.par(2))*scale(m);
                input.(defect{index(m)}.par(1)).(defect{index(m)}.par(2)) = new_par;
                w_info(i,m) = join([defect{index(m)}.par " = " num2str(new_par) ...
                    "| scale" " = " num2str(scale(m))],"");
            else
                new_par = input.(defect{index(m)}.par)*scale(m);
                input.(defect{index(m)}.par) = new_par;
                w_info(i,m) = join([defect{index(m)}.par " = " num2str(new_par) ...
                    "| scale" " = " num2str(scale(m))],"");
            end
        end

        output = VSV_toolbox_RK4(input, tlength, dt, NO_PLOTS);


        w_poor(i) = max(output.progen);
        w_rich(i) = max((output.progen).^(1./output.tt));
    end

    store{j}.poor = w_poor;
    store{j}.rich = w_rich;
    store{j}.info = w_info;


%% Status update
disp(j)

disp(store{j}.poor)
disp(store{j}.rich)
disp(store{j}.info)

end

if ~light
    save(strcat('epistasis_output_',num2str(cluster_num),'_',num2str(process_num),...
        '_severity',num2str(severity),'.mat'), 'store');
elseif light
    save(strcat('epistasis_output_',num2str(cluster_num),'_',num2str(process_num),...
    '_light',num2str(severity),'.mat'), 'store');
end

end

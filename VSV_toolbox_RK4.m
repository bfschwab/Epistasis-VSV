function output = VSV_toolbox_RK4(pars, tlength, dt, NO_PLOTS)



    %use logic to pass parameters to ff functions to check and allocate

    %% Constants
    vc = 1.7*1e-12;     %L
    navo = 6.023*1e23;
    
    q = tlength/dt; %number of steps

    nrib = 5*1e6/navo/vc;
    bound_inc = 10^-7;

    n_n = 1258;
    n_p = 466;
    n_m = 1826;
    n_g = 1205;
    n_l = 50;

    half_host = 5379.3;
    kd_h = 1/half_host*log(2);

    te = .5;
    dec_index = 1;

    %% Defaults, Changeable Parameters
    isub = {'n', 'p', 'm', 'g', 'l'};

    % kdp input or defaults
    kdpisub = {'n', 'p', 'm', 'g', 'l1', 'l2'};
    kdp_emb = struct(); %kdp_emb is an embedded structure with default values
    kdp_emb.n = 3.45*1e-5;   % s-1
    kdp_emb.p = 1.4*1e-6;   % s-1
    kdp_emb.m = 1.5*1e-4;   % s-1
    kdp_emb.g = 5.7*1e-5;   % s-1  
    kdp_emb.l1 = 1.2*1e-5;   % s-1
    kdp_emb.l2 = 4.3*1e-5;   % s-1

    if isfield(pars, 'kdp')
        for i = 1:length(kdpisub)       
            if isfield(pars.kdp, kdpisub{i})
                kdp.(kdpisub{i}) = pars.kdp.(kdpisub{i});
            else
                kdp.(kdpisub{i}) = kdp_emb.(kdpisub{i});
            end
        end
    else 
        kdp = kdp_emb;
    end


    % eta input or defaults
    eta_emb = struct();
    eta_emb.n = 1.0;
    eta_emb.p = 0.75;
    eta_emb.m = 0.75;
    eta_emb.g = 0.75;
    eta_emb.l = 0.05;
    
    if isfield(pars, 'eta')
        for i = 1:length(isub)       
            if isfield(pars.eta, isub{i})
                eta.(isub{i}) = pars.eta.(isub{i});
            else
                eta.(isub{i}) = eta_emb.(isub{i});
            end
        end
    else 
        eta = eta_emb;
    end
    

    % lm inputs or defaults
    lmisub = {'le', 'tr', 'tr', 't'};

    lm_emb = struct();
    lm_emb.le = 47;          %nt, length of leader sequence
    lm_emb.lec = lm_emb.le;       %nt, length of complementary leader region of VSV genome unused !!!!

    lm_emb.tr = 59;          %nt, length of trailer, this partitioning will be used to create delays in the transcription and replication of the genome
    lm_emb.trc = lm_emb.tr;       %nt, length of anti-genomic promoter
    lm_emb.t = 11150;        %nt, sum of all parts. we are missing 11 nt from the total, who knows why?

    lm_emb.n = 1326;        %nt, length of N mRNA
    lm_emb.p = 814;         %nt, length of P mRNA
    lm_emb.m = 831;         %nt, length of M mRNA
    lm_emb.g = 1665;        %nt, length of G mRNA
    lm_emb.l = 6373;        %nt, length of L mRNA
     
    if isfield(pars, 'lm')
        for i = 1:length(isub)       
            if isfield(pars.lm, isub{i})
                lm.(isub{i}) = pars.lm.(isub{i});
            else
                lm.(isub{i}) = lm_emb.(isub{i});
            end
        end
        
        for i = 1:length(lmisub)       
            if isfield(pars.lm, lmisub{i})
                lm.(lmisub{i}) = pars.lm.(lmisub{i});
            else
                lm.(lmisub{i}) = lm_emb.(lmisub{i});
            end
        end
    else 
        lm = lm_emb;
    end
    

    % lp input or defaults
    lp_emb = struct();
    lp_emb.n = 422;  %aa, length of protein N
    lp_emb.p = 265;  %aa, length of protein P
    lp_emb.m = 229;  %aa, length of protein M
    lp_emb.g = 511;  %aa, length of protein G
    lp_emb.l = 2109; %aa, length of protein L
    
    if isfield(pars, 'lp')
        for i = 1:length(isub)       
            if isfield(pars.lp, isub{i})
                lp.(isub{i}) = pars.lp.(isub{i});
            else
                lp.(isub{i}) = lp_emb.(isub{i});
            end
        end
    else 
        lp = lp_emb;
    end

    % nsec input or defaults
    nsecisub = {'le', 'tr', 't'};
    
    nsec_emb = struct();
    nsec_emb.le = 1; %count, number of segments in gene leader
    nsec_emb.tr = 1; %count, number of segments in gene trailer
    nsec_emb.t = 40; %count, number of segments in genome

    nsec_emb.n = 5;  %count, number of segments in gene N
    nsec_emb.p = 3;  %count, number of segments in gene P
    nsec_emb.m = 3;  %count, number of segments in gene M
    nsec_emb.g = 6;  %count, number of segments in gene G
    nsec_emb.l = 23; %count, number of segments in gene L

    if isfield(pars, 'nsec')
        for i = 1:length(isub)       
            if isfield(pars.nsec, isub{i})
                nsec.(isub{i}) = pars.nsec.(isub{i});
            else
                nsec.(isub{i}) = nsec_emb.(isub{i});
            end
        end
        
        for i = 1:length(nsecisub)       
            if isfield(pars.nsec, nsecisub{i})
                nsec.(nsecisub{i}) = pars.nsec.(nsecisub{i});
            else
                nsec.(nsecisub{i}) = nsec_emb.(nsecisub{i});
            end
        end
    else 
        nsec = nsec_emb;
    end

    % other variables not stored in structures

    if isfield(pars, 'kd_m')
        kd_m = pars.kd_m;
    else
        kd_m = 1.9*1e-4;    % s-1
    end

    if isfield(pars, 'kd_nc')
        kd_nc = pars.kd_nc;
    else
        kd_nc = 1.9*1e-5;    % s-1
    end

    if isfield(pars, 'kd_nt')
        kd_nt = pars.kd_nt;
    else
        kd_nt = 1.9*1e-4;  % s-1    
    end

    if isfield(pars, 'ke_pol')
        ke_pol = pars.ke_pol;
    else
        ke_pol = 3.7; % nt/s
    end

    if isfield(pars, 'ke_rib')
        ke_rib = pars.ke_rib;
    else
        ke_rib = 6.0*1.0;       % aa/s
    end

    if isfield(pars, 'p_lratio')
        p_lratio = pars.p_lratio;
    else
        p_lratio = 3.6;    % mol ratio
    end

    if isfield(pars, 'ini_a')
        ini_a = pars.ini_a;
    else
        ini_a = 0.99;     % initial affinity of host mrna to ribosome  
    end

    if isfield(pars, 'sprom')
        sprom = pars.sprom;
    else
        sprom = 5.43;
    end

    if isfield(pars, 'spol')
        spol = pars.spol;
    else
        spol = 170;
    end    

    if isfield(pars, 'srib')
        srib = pars.srib;
    else
        srib = 238.5; 
    end    

    if isfield(pars, 'scond')
        scond = pars.scond;
    else
        scond = 2.77*1e-4*0.2*1.125;   
    end    

    if isfield(pars, 'condm')
        condm = pars.condm;
    else
        condm = 0.1;
    end    

    if isfield(pars, 'k_int')
        k_int = pars.k_int;
    else
        k_int = 0.7194/3600;
    end    


    %%

    % 
    % kd_m = 1.9*1e-4;    % s-1
    % kd_nc = 1.9*1e-5;   % s-1   
    % kd_nt = 1.9*1e-4;  % s-1    
    % ke_pol = 3.7;       % nt/s
    % ke_rib = 6.0*1.0;       % aa/s   
    % p_lratio = 3.6;    % mol ratio
    % ini_a = 0.99;     % initial affinity of host mrna to ribosome  
    % dec_cri = ini_a*0.999999;      
    % dec_index = 1;
    % sprom = 5.43;
    % spol = 170;
    % srib = 238.5; 
    % scond = 2.77*1e-4*0.2*1.125;     
    % condm = 0.1;
    % k_int = 0.7194;
    % k_int = k_int/3600; 
    % binding_t = 1; %hr  
    % moi = 3;
    % bound_inc = 10^-7; 


    %% Initial Condition

    pn = NaN(1,q);
    pp = NaN(1,q);
    pm = NaN(1,q);
    pg = NaN(1,q);
    pl = NaN(1,q);
    pn(1) = 0;    
    pp(1) = n_p/vc/navo*bound_inc;    
    pm(1) = n_m/vc/navo*bound_inc;    
    pg(1) = n_g/vc/navo*bound_inc;    
    pl(1) = n_l/vc/navo*bound_inc;  

    n_nc = NaN(1,q);
    p_nc = NaN(1,q);
    n_nc(1) = 1/vc/navo*bound_inc;   
    p_nc(1) = 0;

    n_cap = NaN(1,q);
    cap = NaN(1,q);
    n_cap(1) = 1;
    cap(1) = 1 - n_cap(1);

    mrna_n = NaN(1,q);
    mrna_p = NaN(1,q);
    mrna_m = NaN(1,q);
    mrna_g = NaN(1,q);
    mrna_l = NaN(1,q);
    nt_new = NaN(1,q);
    mrna_n(1) = 0;
    mrna_p(1) = 0;
    mrna_m(1) = 0;
    mrna_g(1) = 0;
    mrna_l(1) = 0;
    nt_new(1) = 0;

    rib = NaN(1,q);
    rib_n = NaN(1,q);
    rib_p = NaN(1,q);
    rib_m = NaN(1,q);
    rib_g = NaN(1,q);
    rib_l = NaN(1,q);
    rib(1) = 0;
    rib_n(1) = 0;
    rib_p(1) = 0;
    rib_m(1) = 0;
    rib_g(1) = 0;
    rib_l(1) = 0;

    nt = NaN(1,q);
    proge = NaN(1,q);
    progen = NaN(1,q);
    n_nc_m = NaN(1,q);
    rel = NaN(1,q);
    relm = NaN(1,q);
    n_nc_pol = NaN(1,q);
    d_pol_le = NaN(1,q);
    d_pol_term = NaN(1,q);
    d_pol_trc = NaN(1,q);
    pol_r_t = NaN(1,q);

    nt(1) = 0;
    proge(1) = 0;
    progen(1) = 0;
    n_nc_m(1) = 0;
    rel(1) = 0;
    relm(1) = 0;
    d_pol_trc(1) = 0;

    pol_n = NaN(nsec.n,q);
    pol_p = NaN(nsec.p,q);
    pol_m = NaN(nsec.m,q);
    pol_g = NaN(nsec.g,q);
    pol_l = NaN(nsec.l,q);
    pol_n(:,1) = zeros(nsec.n,1);
    pol_p(:,1) = zeros(nsec.p,1);
    pol_m(:,1) = zeros(nsec.m,1);
    pol_g(:,1) = zeros(nsec.g,1);
    pol_l(:,1) = zeros(nsec.l,1);

    poln_n = NaN(nsec.n,q);
    poln_p = NaN(nsec.p,q);
    poln_m = NaN(nsec.m,q);
    poln_g = NaN(nsec.g,q);
    poln_l = NaN(nsec.l,q);
    poln_n(:,1) = zeros(nsec.n,1);
    poln_p(:,1) = zeros(nsec.p,1);
    poln_m(:,1) = zeros(nsec.m,1);
    poln_g(:,1) = zeros(nsec.g,1);
    poln_l(:,1) = zeros(nsec.l,1);

    pol_tn = NaN(1,q);
    pol_tp = NaN(1,q);
    pol_tm = NaN(1,q);
    pol_tg = NaN(1,q);
    pol_tl = NaN(1,q);
    pol_tn(1) = 0;
    pol_tp(1) = 0;
    pol_tm(1) = 0;
    pol_tg(1) = 0;
    pol_tl(1) = 0;

    pol_le = NaN(1,q);
    pol_trc = NaN(1,q);
    polp_lec = NaN(1,q);
    pol_t = NaN(1,q);
    pol_f = NaN(1,q);
    pol_le(1) = 0;
    pol_trc(1) = 0;
    polp_lec(1) = 0;
    pol_t(1) = 0;
    pol_f(1) = 0;

    poln_tr = NaN(1,q);
    polp = NaN(nsec.t,q);
    poln_tr(1) = 0;
    polp(:,1) = zeros(nsec.t,1);

    pol_le(1) = 0/vc/navo;
    pol_trc(1) = 0/vc/navo;
    polp_lec(1) = 0/vc/navo;
    poln_tr(1) = 0/vc/navo;
    pol_t(1) = 0/vc/navo;

    rib_host = NaN(1,q);
    f_dec = NaN(1,q);
    v_ex = NaN(1,q);
    rib_host(1) = ini_a;
    v_ex(1) = .3;  

    %%

    min_pm = pm(1);

    if pp(1)/pl(1) > p_lratio
        pol_f(1) = 1.0;
    else
        pol_f(1) = pp(1)/pl(1)/p_lratio;
    end

    n_nc_pol(1) = (n_nc(1)-pol_le(1)*spol/lm.le)/(scond*(1.-condm)*pm(1)/(pl(1)*pol_f(1)...
       -pol_t(1))+1.);

    if (n_nc_pol(1)*lm.le/spol+(p_nc(1)*lm.tr/spol*sprom-pol_trc(1))) < (pl(1)*pol_f(1)-pol_t(1))
       d_pol_term(1) = n_nc_pol(1)*lm.le/spol+(p_nc(1)*lm.tr/spol*sprom-pol_trc(1));
    else
       d_pol_term(1) = (pl(1)*pol_f(1)-pol_t(1));
    end

    d_pol_le(1) = n_nc_pol(1)*lm.le/spol/(n_nc_pol(1)*lm.le/spol+(p_nc(1)*lm.tr/spol*sprom-pol_trc(1)))...
       *d_pol_term(1);
    pol_le(1) = pol_le(1)+d_pol_le(1);
    pol_trc(1) = pol_trc(1)+(d_pol_term(1)-d_pol_le(1));
    pol_r_t(1) = pol_trc(1)+pol_le(1);

    for i = 1:q
       if i*dt/3600 > 1    
            k_int = 0;      
       end 

       rib_host(1) = ini_a;
       mx = (ini_a*100/1.88*1e-5)^(-1./0.9582);

       if (pm(i)*navo*vc < (v_ex(1)-v_ex(i))*n_m) && (i > 1) && pm(i) < pm(i-1)   
          min_pm = pm(i);
          rib_host(i) = ini_a;
       end

       if  (i > 1) && (pm(i) > pm(i-1))
          if (pm(i)-min_pm)*(1-condm)*navo*vc < mx
             rib_host(i) = ini_a;
          else
             rib_host(i) = 1.88*1e5*((pm(i)-min_pm)*(1-condm)*navo*vc)^-0.9582/100;
          end
       end

       if i > 1 && (pm(i-1) > pm(i))
          rib_host(i) = rib_host(i-1);
       end

       f_dec(i) = exp(-1*kd_h*(i-dec_index)/q*tlength );


       if pp(i) > pl(i)*10    % For considering the stability of L depending on the level of pp
          kdp.l = kdp.l1;
       else
          kdp.l = kdp.l2;
       end

       k1_ent = dt*f_ent_emb(k_int,v_ex(i));   
       k2_ent = dt*f_ent_emb(k_int,v_ex(i)+k1_ent/2);   
       k3_ent = dt*f_ent_emb(k_int,v_ex(i)+k2_ent/2);   
       k4_ent = dt*f_ent_emb(k_int,v_ex(i)+k2_ent/2);   
       v_ex(i+1) = v_ex(i) + 1/6*(k1_ent + 2*k2_ent + 2*k3_ent + k4_ent); 

       if v_ex(i+1) < 0   
           v_ex(i+1) = 0;
       end    

       bound_inc = -1*1/6*(k1_ent + 2*k2_ent + 2*k3_ent + k4_ent)/vc/navo;  

       k11 = dt*f1_emb(ke_pol,eta.n,n_cap(i),pol_le(i),lm.le,nsec.n,lm.n,pol_n(1,i));
       k21 = dt*f1_emb(ke_pol,eta.n,n_cap(i),pol_le(i),lm.le,nsec.n,lm.n,pol_n(1,i)+k11/2);
       k31 = dt*f1_emb(ke_pol,eta.n,n_cap(i),pol_le(i),lm.le,nsec.n,lm.n,pol_n(1,i)+k21/2);
       k41 = dt*f1_emb(ke_pol,eta.n,n_cap(i),pol_le(i),lm.le,nsec.n,lm.n,pol_n(1,i)+k31);
       pol_n(1,i+1) = pol_n(1,i) +1/6*(k11 + 2*k21 + 2*k31 + k41);
       k12(1) = k11;
       k22(1) = k21;
       k32(1) = k31;
       k42(1) = k41;


       for j = 2:nsec.n

          k12(j) = dt*f2_emb(ke_pol,nsec.n,lm.n,pol_n(j-1,i),pol_n(j,i));
          k22(j) = dt*f2_emb(ke_pol,nsec.n,lm.n,pol_n(j-1,i)+k12(j-1)/2,pol_n(j,i)+k12(j)/2);
          k32(j) = dt*f2_emb(ke_pol,nsec.n,lm.n,pol_n(j-1,i)+k22(j-1)/2,pol_n(j,i)+k22(j)/2);
          k42(j) = dt*f2_emb(ke_pol,nsec.n,lm.n,pol_n(j-1,i)+k32(j-1),pol_n(j,i)+k32(j));
          pol_n(j,i+1) = pol_n(j,i) +1/6*(k12(j) + 2*k22(j) + 2*k32(j) + k42(j));     

       end

       k13 = dt*f3_emb(ke_pol,eta.p,nsec.n,lm.n,pol_n(nsec.n,i),nsec.p,lm.p,pol_p(1,i));
       k23 = dt*f3_emb(ke_pol,eta.p,nsec.n,lm.n,pol_n(nsec.n,i)+k12(nsec.n)/2,...
             nsec.p,lm.p,pol_p(1,i)+k13/2);
       k33 = dt*f3_emb(ke_pol,eta.p,nsec.n,lm.n,pol_n(nsec.n,i)+k22(nsec.n)/2,...
             nsec.p,lm.p,pol_p(1,i)+k23/2);
       k43 = dt*f3_emb(ke_pol,eta.p,nsec.n,lm.n,pol_n(nsec.n,i)+k32(nsec.n),...
             nsec.p,lm.p,pol_p(1,i)+k33);
       pol_p(1,i+1) = pol_p(1,i) +1/6*(k13 + 2*k23 + 2*k33 + k43);
       k12(1) = k13;
       k22(1) = k23;
       k32(1) = k33;
       k42(1) = k43;


       for j = 2:nsec.p
            k12(j) = dt*f2_emb(ke_pol,nsec.p,lm.p,pol_p(j-1,i),pol_p(j,i));
            k22(j) = dt*f2_emb(ke_pol,nsec.p,lm.p,pol_p(j-1,i)+k12(j-1)/2,pol_p(j,i)+k12(j)/2);
            k32(j) = dt*f2_emb(ke_pol,nsec.p,lm.p,pol_p(j-1,i)+k22(j-1)/2,pol_p(j,i)+k22(j)/2);
            k42(j) = dt*f2_emb(ke_pol,nsec.p,lm.p,pol_p(j-1,i)+k32(j-1),pol_p(j,i)+k32(j));
            pol_p(j,i+1) = pol_p(j,i) +1/6*(k12(j) + 2*k22(j) + 2*k32(j) + k42(j));  
       end

       k13 = dt*f3_emb(ke_pol,eta.m,nsec.p,lm.p,pol_p(nsec.p,i),nsec.m,lm.m,pol_m(1,i));
       k23 = dt*f3_emb(ke_pol,eta.m,nsec.p,lm.p,pol_p(nsec.p,i)+k12(nsec.p)/2,...
             nsec.m,lm.m,pol_m(1,i)+k13/2);
       k33 = dt*f3_emb(ke_pol,eta.m,nsec.p,lm.p,pol_p(nsec.p,i)+k22(nsec.p)/2,...
             nsec.m,lm.m,pol_m(1,i)+k23/2);
       k43 = dt*f3_emb(ke_pol,eta.m,nsec.p,lm.p,pol_p(nsec.p,i)+k32(nsec.p),...
             nsec.m,lm.m,pol_m(1,i)+k33);
       pol_m(1,i+1) = pol_m(1,i) +1/6*(k13 + 2*k23 + 2*k33 + k43);
       k12(1) = k13;
       k22(1) = k23;
       k32(1) = k33;
       k42(1) = k43;


       for j = 2:nsec.m
        k12(j) = dt*f2_emb(ke_pol,nsec.m,lm.m,pol_m(j-1,i),pol_m(j,i));
        k22(j) = dt*f2_emb(ke_pol,nsec.m,lm.m,pol_m(j-1,i)+k12(j-1)/2,pol_m(j,i)+k12(j)/2);
        k32(j) = dt*f2_emb(ke_pol,nsec.m,lm.m,pol_m(j-1,i)+k22(j-1)/2,pol_m(j,i)+k22(j)/2);
        k42(j) = dt*f2_emb(ke_pol,nsec.m,lm.m,pol_m(j-1,i)+k32(j-1),pol_m(j,i)+k32(j));
        pol_m(j,i+1) = pol_m(j,i) +1/6*(k12(j) + 2*k22(j) + 2*k32(j) + k42(j));  
       end

       k13 = dt*f3_emb(ke_pol,eta.g,nsec.m,lm.m,pol_m(nsec.m,i),nsec.g,lm.g,pol_g(1,i));
       k23 = dt*f3_emb(ke_pol,eta.g,nsec.m,lm.m,pol_m(nsec.m,i)+k12(nsec.m)/2,...
             nsec.g,lm.g,pol_g(1,i)+k13/2);
       k33 = dt*f3_emb(ke_pol,eta.g,nsec.m,lm.m,pol_m(nsec.m,i)+k22(nsec.m)/2,...
             nsec.g,lm.g,pol_g(1,i)+k23/2);
       k43 = dt*f3_emb(ke_pol,eta.g,nsec.m,lm.m,pol_m(nsec.m,i)+k32(nsec.m),...
             nsec.g,lm.g,pol_g(1,i)+k33);
       pol_g(1,i+1) = pol_g(1,i) +1/6*(k13 + 2*k23 + 2*k33 + k43);
       k12(1) = k13;
       k22(1) = k23;
       k32(1) = k33;
       k42(1) = k43;


       for j = 2:nsec.g
        k12(j) = dt*f2_emb(ke_pol,nsec.g,lm.g,pol_g(j-1,i),pol_g(j,i));
        k22(j) = dt*f2_emb(ke_pol,nsec.g,lm.g,pol_g(j-1,i)+k12(j-1)/2,pol_g(j,i)+k12(j)/2);
        k32(j) = dt*f2_emb(ke_pol,nsec.g,lm.g,pol_g(j-1,i)+k22(j-1)/2,pol_g(j,i)+k22(j)/2);
        k42(j) = dt*f2_emb(ke_pol,nsec.g,lm.g,pol_g(j-1,i)+k32(j-1),pol_g(j,i)+k32(j));
        pol_g(j,i+1) = pol_g(j,i) +1/6*(k12(j) + 2*k22(j) + 2*k32(j) + k42(j)); 
       end

       k13 = dt*f3_emb(ke_pol,eta.l,nsec.g,lm.g,pol_g(nsec.g,i),nsec.l,lm.l,pol_l(1,i));
       k23 = dt*f3_emb(ke_pol,eta.l,nsec.g,lm.g,pol_g(nsec.g,i)+k12(nsec.g)/2,...
             nsec.l,lm.l,pol_l(1,i)+k13/2);
       k33 = dt*f3_emb(ke_pol,eta.l,nsec.g,lm.g,pol_g(nsec.g,i)+k22(nsec.g)/2,...
             nsec.l,lm.l,pol_l(1,i)+k23/2);
       k43 = dt*f3_emb(ke_pol,eta.l,nsec.g,lm.g,pol_g(nsec.g,i)+k32(nsec.g),...
             nsec.l,lm.l,pol_l(1,i)+k33);
       pol_l(1,i+1) = pol_l(1,i) +1/6*(k13 + 2*k23 + 2*k33 + k43);
       k12(1) = k13;
       k22(1) = k23;
       k32(1) = k33;
       k42(1) = k43;


       for j = 2:nsec.l
        k12(j) = dt*f2_emb(ke_pol,nsec.l,lm.l,pol_l(j-1,i),pol_l(j,i));
        k22(j) = dt*f2_emb(ke_pol,nsec.l,lm.l,pol_l(j-1,i)+k12(j-1)/2,pol_l(j,i)+k12(j)/2);
        k32(j) = dt*f2_emb(ke_pol,nsec.l,lm.l,pol_l(j-1,i)+k22(j-1)/2,pol_l(j,i)+k22(j)/2);
        k42(j) = dt*f2_emb(ke_pol,nsec.l,lm.l,pol_l(j-1,i)+k32(j-1),pol_l(j,i)+k32(j));
        pol_l(j,i+1) = pol_l(j,i) +1/6*(k12(j) + 2*k22(j) + 2*k32(j) + k42(j)); 
       end

       k14 = dt*f4_emb(ke_pol,cap(i),pol_le(i),lm.le,nsec.n,lm.n,poln_n(1,i));
       k24 = dt*f4_emb(ke_pol,cap(i),pol_le(i),lm.le,nsec.n,lm.n,poln_n(1,i)+k14/2);
       k34 = dt*f4_emb(ke_pol,cap(i),pol_le(i),lm.le,nsec.n,lm.n,poln_n(1,i)+k24/2);
       k44 = dt*f4_emb(ke_pol,cap(i),pol_le(i),lm.le,nsec.n,lm.n,poln_n(1,i)+k34);
       poln_n(1,i+1) = poln_n(1,i) +1/6*(k14 + 2*k24 + 2*k34 + k44);
       k15(1) = k14;
       k25(1) = k24;
       k35(1) = k34;
       k45(1) = k44;


       for j = 2:nsec.n
        k15(j) = dt*f5_emb(ke_pol,nsec.n,lm.n,poln_n(j-1,i),poln_n(j,i));
        k25(j) = dt*f5_emb(ke_pol,nsec.n,lm.n,poln_n(j-1,i)+k15(j-1)/2,poln_n(j,i)+k15(j)/2);
        k35(j) = dt*f5_emb(ke_pol,nsec.n,lm.n,poln_n(j-1,i)+k25(j-1)/2,poln_n(j,i)+k25(j)/2);
        k45(j) = dt*f5_emb(ke_pol,nsec.n,lm.n,poln_n(j-1,i)+k35(j-1),poln_n(j,i)+k35(j));
        poln_n(j,i+1) = poln_n(j,i) +1/6*(k15(j) + 2*k25(j) + 2*k35(j) + k45(j)); 
       end

       k16 = dt*f6_emb(ke_pol,cap(i),nsec.n,lm.n,poln_n(nsec.n,i),nsec.p,lm.p,poln_p(1,i));
       k26 = dt*f6_emb(ke_pol,cap(i),nsec.n,lm.n,poln_n(nsec.n,i)+k15(nsec.n)/2,...
             nsec.p,lm.p,poln_p(1,i)+k16/2);
       k36 = dt*f6_emb(ke_pol,cap(i),nsec.n,lm.n,poln_n(nsec.n,i)+k25(nsec.n)/2,...
             nsec.p,lm.p,poln_p(1,i)+k26/2);
       k46 = dt*f6_emb(ke_pol,cap(i),nsec.n,lm.n,poln_n(nsec.n,i)+k35(nsec.n),...
             nsec.p,lm.p,poln_p(1,i)+k36);
       poln_p(1,i+1) = poln_p(1,i) +1/6*(k16 + 2*k26 + 2*k36 + k46);
       k15(1) = k16;
       k25(1) = k26;
       k35(1) = k36;
       k45(1) = k46;


       for j = 2:nsec.p
        k15(j) = dt*f5_emb(ke_pol,nsec.p,lm.p,poln_p(j-1,i),poln_p(j,i));
        k25(j) = dt*f5_emb(ke_pol,nsec.p,lm.p,poln_p(j-1,i)+k15(j-1)/2,poln_p(j,i)+k15(j)/2);
        k35(j) = dt*f5_emb(ke_pol,nsec.p,lm.p,poln_p(j-1,i)+k25(j-1)/2,poln_p(j,i)+k25(j)/2);
        k45(j) = dt*f5_emb(ke_pol,nsec.p,lm.p,poln_p(j-1,i)+k35(j-1),poln_p(j,i)+k35(j));
        poln_p(j,i+1) = poln_p(j,i) +1/6*(k15(j) + 2*k25(j) + 2*k35(j) + k45(j)); 
       end

       k16 = dt*f6_emb(ke_pol,cap(i),nsec.p,lm.p,poln_p(nsec.p,i),nsec.m,lm.m,poln_m(1,i));
       k26 = dt*f6_emb(ke_pol,cap(i),nsec.p,lm.p,poln_p(nsec.p,i)+k15(nsec.p)/2,...
             nsec.m,lm.m,poln_m(1,i)+k16/2);
       k36 = dt*f6_emb(ke_pol,cap(i),nsec.p,lm.p,poln_p(nsec.p,i)+k25(nsec.p)/2,...
             nsec.m,lm.m,poln_m(1,i)+k26/2);
       k46 = dt*f6_emb(ke_pol,cap(i),nsec.p,lm.p,poln_p(nsec.p,i)+k35(nsec.p),...
             nsec.m,lm.m,poln_m(1,i)+k36);
       poln_m(1,i+1) = poln_m(1,i) +1/6*(k16 + 2*k26 + 2*k36 + k46);
       k15(1) = k16;
       k25(1) = k26;
       k35(1) = k36;
       k45(1) = k46;

       for j = 2:nsec.m
       k15(j) = dt*f5_emb(ke_pol,nsec.m,lm.m,poln_m(j-1,i),poln_m(j,i));
       k25(j) = dt*f5_emb(ke_pol,nsec.m,lm.m,poln_m(j-1,i)+k15(j-1)/2,poln_m(j,i)+k15(j)/2);
       k35(j) = dt*f5_emb(ke_pol,nsec.m,lm.m,poln_m(j-1,i)+k25(j-1)/2,poln_m(j,i)+k25(j)/2);
       k45(j) = dt*f5_emb(ke_pol,nsec.m,lm.m,poln_m(j-1,i)+k35(j-1),poln_m(j,i)+k35(j));
       poln_m(j,i+1) = poln_m(j,i) +1/6*(k15(j) + 2*k25(j) + 2*k35(j) + k45(j)); 
       end

       k16 = dt*f6_emb(ke_pol,cap(i),nsec.m,lm.m,poln_m(nsec.m,i),nsec.g,lm.g,poln_g(1,i));
       k26 = dt*f6_emb(ke_pol,cap(i),nsec.m,lm.m,poln_m(nsec.m,i)+k15(nsec.m)/2,...
             nsec.g,lm.g,poln_g(1,i)+k16/2);
       k36 = dt*f6_emb(ke_pol,cap(i),nsec.m,lm.m,poln_m(nsec.m,i)+k25(nsec.m)/2,...
             nsec.g,lm.g,poln_g(1,i)+k26/2);
       k46 = dt*f6_emb(ke_pol,cap(i),nsec.m,lm.m,poln_m(nsec.m,i)+k35(nsec.m),...
             nsec.g,lm.g,poln_g(1,i)+k36);
       poln_g(1,i+1) = poln_g(1,i) +1/6*(k16 + 2*k26 + 2*k36 + k46);
       k15(1) = k16;
       k25(1) = k26;
       k35(1) = k36;
       k45(1) = k46;

       for j = 2:nsec.g
           k15(j) = dt*f5_emb(ke_pol,nsec.g,lm.g,poln_g(j-1,i),poln_g(j,i));
           k25(j) = dt*f5_emb(ke_pol,nsec.g,lm.g,poln_g(j-1,i)+k15(j-1)/2,poln_g(j,i)+k15(j)/2);
           k35(j) = dt*f5_emb(ke_pol,nsec.g,lm.g,poln_g(j-1,i)+k25(j-1)/2,poln_g(j,i)+k25(j)/2);
           k45(j) = dt*f5_emb(ke_pol,nsec.g,lm.g,poln_g(j-1,i)+k35(j-1),poln_g(j,i)+k35(j));
           poln_g(j,i+1) = poln_g(j,i) +1/6*(k15(j) + 2*k25(j) + 2*k35(j) + k45(j)); 
       end

       k16 = dt*f6_emb(ke_pol,cap(i),nsec.g,lm.g,poln_g(nsec.g,i),nsec.l,lm.l,poln_l(1,i));
       k26 = dt*f6_emb(ke_pol,cap(i),nsec.g,lm.g,poln_g(nsec.g,i)+k15(nsec.g)/2,...
             nsec.l,lm.l,poln_l(1,i)+k16/2);
       k36 = dt*f6_emb(ke_pol,cap(i),nsec.g,lm.g,poln_g(nsec.g,i)+k25(nsec.g)/2,...
             nsec.l,lm.l,poln_l(1,i)+k26/2);
       k46 = dt*f6_emb(ke_pol,cap(i),nsec.g,lm.g,poln_g(nsec.g,i)+k35(nsec.g),...
             nsec.l,lm.l,poln_l(1,i)+k36);
       poln_l(1,i+1) = poln_l(1,i) +1/6*(k16 + 2*k26 + 2*k36 + k46);
       k15(1) = k16;
       k25(1) = k26;
       k35(1) = k36;
       k45(1) = k46;

       for j = 2:nsec.l
           k15(j) = dt*f5_emb(ke_pol,nsec.l,lm.l,poln_l(j-1,i),poln_l(j,i));
           k25(j) = dt*f5_emb(ke_pol,nsec.l,lm.l,poln_l(j-1,i)+k15(j-1)/2,poln_l(j,i)+k15(j)/2);
           k35(j) = dt*f5_emb(ke_pol,nsec.l,lm.l,poln_l(j-1,i)+k25(j-1)/2,poln_l(j,i)+k25(j)/2);
           k45(j) = dt*f5_emb(ke_pol,nsec.l,lm.l,poln_l(j-1,i)+k35(j-1),poln_l(j,i)+k35(j));
           poln_l(j,i+1) = poln_l(j,i) +1/6*(k15(j) + 2*k25(j) + 2*k35(j) + k45(j)); 
       end

       k17 = dt*f7_emb(ke_pol,nsec.l,lm.l,poln_l(nsec.l,i),poln_tr(i),lm.tr);
       k27 = dt*f7_emb(ke_pol,nsec.l,lm.l,poln_l(nsec.l,i)+k15(nsec.l)/2,...
             poln_tr(i)+k17/2,lm.tr);
       k37 = dt*f7_emb(ke_pol,nsec.l,lm.l,poln_l(nsec.l,i)+k25(nsec.l)/2,...
              poln_tr(i)+k27/2,lm.tr);
       k47 = dt*f7_emb(ke_pol,nsec.l,lm.l,poln_l(nsec.l,i)+k35(nsec.l),...
             poln_tr(i)+k37,lm.tr);
       poln_tr(i+1) = poln_tr(i) +1/6*(k17 + 2*k27 + 2*k37 + k47);

        k18 = dt*f8_emb(ke_pol,cap(i),pol_trc(i),lm.tr,nsec.t,lm.t,lm.le,polp(1,i));
        k28 = dt*f8_emb(ke_pol,cap(i),pol_trc(i),lm.tr,nsec.t,lm.t,lm.le,polp(1,i)+k18/2);
        k38 = dt*f8_emb(ke_pol,cap(i),pol_trc(i),lm.tr,nsec.t,lm.t,lm.le,polp(1,i)+k28/2);
        k48 = dt*f8_emb(ke_pol,cap(i),pol_trc(i),lm.tr,nsec.t,lm.t,lm.le,polp(1,i)+k38);
        polp(1,i+1) = polp(1,i) +1/6*(k18 + 2*k28 + 2*k38 + k48);
        k19(1) = k18;
        k29(1) = k28;
        k39(1) = k38;
        k49(1) = k48;

    for j = 2:nsec.t
      k19(j) = dt*f9_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(j-1,i),polp(j,i));
      k29(j) = dt*f9_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(j-1,i)+k19(j-1)/2,polp(j,i)+k19(j)/2);
      k39(j) = dt*f9_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(j-1,i)+k29(j-1)/2,polp(j,i)+k29(j)/2);
      k49(j) = dt*f9_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(j-1,i)+k39(j-1),polp(j,i)+k39(j));
      polp(j,i+1) = polp(j,i) +1/6*(k19(j) + 2*k29(j) + 2*k39(j) + k49(j));
    end

     k19_2 = dt*f9_2_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(nsec.t,i),polp_lec(i));
       k29_2 = dt*f9_2_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(nsec.t,i)+k19(nsec.t)/2,polp_lec(i)+k19_2/2);
       k39_2 = dt*f9_2_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(nsec.t,i)+k29(nsec.t)/2,polp_lec(i)+k29_2/2);
       k49_2 = dt*f9_2_emb(ke_pol,nsec.t,lm.t,lm.tr,lm.le,polp(nsec.t,i)+k39(nsec.t),polp_lec(i)+k39_2);
       polp_lec(i+1) = polp_lec(i) +1/6*(k19_2 + 2*k29_2 + 2*k39_2 + k49_2);


      k110 = dt*f10_emb(ke_pol,cap(i),poln_tr(i),lm.tr,kd_nc,p_nc(i));
      k210 = dt*f10_emb(ke_pol,cap(i),poln_tr(i)+k17/2,lm.tr,kd_nc,p_nc(i)+k110/2);
      k310 = dt*f10_emb(ke_pol,cap(i),poln_tr(i)+k27/2,lm.tr,kd_nc,p_nc(i)+k210/2);
      k410 = dt*f10_emb(ke_pol,cap(i),poln_tr(i)+k37,lm.tr,kd_nc,p_nc(i)+k310);
      p_nc(i+1) = p_nc(i) +1/6*(k110 + 2*k210 + 2*k310 + k410);

    k111 = dt*f11_emb(ke_pol,cap(i),polp_lec(i),lm.le,kd_nc,n_nc(i));
    k211 = dt*f11_emb(ke_pol,cap(i),polp_lec(i)+k19_2/2,lm.le,kd_nc,n_nc(i)+k111/2);
    k311 = dt*f11_emb(ke_pol,cap(i),polp_lec(i)+k29_2/2,lm.le,kd_nc,n_nc(i)+k211/2);
    k411 = dt*f11_emb(ke_pol,cap(i),polp_lec(i)+k39_2,lm.le,kd_nc,n_nc(i)+k311);
    n_nc(i+1) = n_nc(i) +1/6*(k111 + 2*k211 + 2*k311 + k411);
    n_nc(i+1) = n_nc(i+1) + bound_inc;  



       pol_tn(i+1) = 0/vc/navo;
       pol_tp(i+1) = 0/vc/navo;
       pol_tm(i+1) = 0/vc/navo;
       pol_tg(i+1) = 0/vc/navo;
       pol_tl(i+1) = 0/vc/navo;

       for j = 1:nsec.n
        pol_tn(i+1) = pol_tn(i+1)+pol_n(j,i+1) ;
       end
       for j = 1:nsec.p
        pol_tp(i+1) = pol_tp(i+1)+pol_p(j,i+1) ;
       end
       for j = 1:nsec.m
        pol_tm(i+1) = pol_tm(i+1)+pol_m(j,i+1) ;
       end
       for j = 1:nsec.g
        pol_tg(i+1) = pol_tg(i+1)+pol_g(j,i+1) ;
       end
       for j = 1:nsec.l
        pol_tl(i+1) = pol_tl(i+1)+pol_l(j,i+1) ;
       end

    k112 = dt*f12_emb(ke_pol,lm.n,pol_tn(i),kd_m,mrna_n(i));
    k212 = dt*f12_emb(ke_pol,lm.n,pol_tn(i),kd_m,mrna_n(i)+k112/2);
    k312 = dt*f12_emb(ke_pol,lm.n,pol_tn(i),kd_m,mrna_n(i)+k212/2);
    k412 = dt*f12_emb(ke_pol,lm.n,pol_tn(i),kd_m,mrna_n(i)+k312);
    mrna_n(i+1) = mrna_n(i) +1/6*(k112 + 2*k212 + 2*k312 + k412);

    k112 = dt*f12_emb(ke_pol,lm.p,pol_tp(i),kd_m,mrna_p(i));
    k212 = dt*f12_emb(ke_pol,lm.p,pol_tp(i),kd_m,mrna_p(i)+k112/2);
    k312 = dt*f12_emb(ke_pol,lm.p,pol_tp(i),kd_m,mrna_p(i)+k212/2);
    k412 = dt*f12_emb(ke_pol,lm.p,pol_tp(i),kd_m,mrna_p(i)+k312);
    mrna_p(i+1) = mrna_p(i) +1/6*(k112 + 2*k212 + 2*k312 + k412);

    k112 = dt*f12_emb(ke_pol,lm.m,pol_tm(i),kd_m,mrna_m(i));
    k212 = dt*f12_emb(ke_pol,lm.m,pol_tm(i),kd_m,mrna_m(i)+k112/2);
    k312 = dt*f12_emb(ke_pol,lm.m,pol_tm(i),kd_m,mrna_m(i)+k212/2);
    k412 = dt*f12_emb(ke_pol,lm.m,pol_tm(i),kd_m,mrna_m(i)+k312);
    mrna_m(i+1) = mrna_m(i) +1/6*(k112 + 2*k212 + 2*k312 + k412);

    k112 = dt*f12_emb(ke_pol,lm.g,pol_tg(i),kd_m,mrna_g(i));
    k212 = dt*f12_emb(ke_pol,lm.g,pol_tg(i),kd_m,mrna_g(i)+k112/2);
    k312 = dt*f12_emb(ke_pol,lm.g,pol_tg(i),kd_m,mrna_g(i)+k212/2);
    k412 = dt*f12_emb(ke_pol,lm.g,pol_tg(i),kd_m,mrna_g(i)+k312);
    mrna_g(i+1) = mrna_g(i) +1/6*(k112 + 2*k212 + 2*k312 + k412);

    k112 = dt*f12_emb(ke_pol,lm.l,pol_tl(i),kd_m,mrna_l(i));
    k212 = dt*f12_emb(ke_pol,lm.l,pol_tl(i),kd_m,mrna_l(i)+k112/2);
    k312 = dt*f12_emb(ke_pol,lm.l,pol_tl(i),kd_m,mrna_l(i)+k212/2);
    k412 = dt*f12_emb(ke_pol,lm.l,pol_tl(i),kd_m,mrna_l(i)+k312);
    mrna_l(i+1) = mrna_l(i) +1/6*(k112 + 2*k212 + 2*k312 + k412);


       dete = (lm.n*mrna_n(i+1)+lm.p*mrna_p(i+1)+lm.m*mrna_m(i+1)+lm.g*mrna_g(i+1)...
          +lm.l*mrna_l(i+1))/srib;
       if nrib*(1-rib_host(i))*f_dec(i) < dete
          rib(i+1) = nrib*(1-rib_host(i))*f_dec(i);
       else
          rib(i+1) = dete;
       end
       if dete == 0 
          rib_n(i+1) = 0/vc/navo;
          rib_p(i+1) = 0/vc/navo;
          rib_m(i+1) = 0/vc/navo;
          rib_g(i+1) = 0/vc/navo;
          rib_l(i+1) = 0/vc/navo;
        else   
          rib_n(i+1) = lm.n*mrna_n(i+1)/(dete*srib)*rib(i+1); 
          rib_p(i+1) = lm.p*mrna_p(i+1)/(dete*srib)*rib(i+1); 
          rib_m(i+1) = lm.m*mrna_m(i+1)/(dete*srib)*rib(i+1); 
          rib_g(i+1) = lm.g*mrna_g(i+1)/(dete*srib)*rib(i+1); 
          rib_l(i+1) = lm.l*mrna_l(i+1)/(dete*srib)*rib(i+1); 
        end

      k113 = dt*f13_emb(ke_rib,lp.p,rib_p(i),kdp.p,pp(i));
      k213 = dt*f13_emb(ke_rib,lp.p,rib_p(i),kdp.p,pp(i)+k113/2);
      k313 = dt*f13_emb(ke_rib,lp.p,rib_p(i),kdp.p,pp(i)+k213/2);
      k413 = dt*f13_emb(ke_rib,lp.p,rib_p(i),kdp.p,pp(i)+k313);
      pp(i+1) = pp(i) +1/6*(k113 + 2*k213 + 2*k313 + k413);
      pp(i+1) = pp(i+1) + n_p*bound_inc; 

      k113 = dt*f13_emb(ke_rib,lp.m,rib_m(i),kdp.m,pm(i));
      k213 = dt*f13_emb(ke_rib,lp.m,rib_m(i),kdp.m,pm(i)+k113/2);
      k313 = dt*f13_emb(ke_rib,lp.m,rib_m(i),kdp.m,pm(i)+k213/2);
      k413 = dt*f13_emb(ke_rib,lp.m,rib_m(i),kdp.m,pm(i)+k313);
      pm(i+1) = pm(i) +1/6*(k113 + 2*k213 + 2*k313 + k413);
      pm(i+1) = pm(i+1) + n_m*bound_inc; 

      k113 = dt*f13_emb(ke_rib,lp.g,rib_g(i),kdp.g,pg(i));
      k213 = dt*f13_emb(ke_rib,lp.g,rib_g(i),kdp.g,pg(i)+k113/2);
      k313 = dt*f13_emb(ke_rib,lp.g,rib_g(i),kdp.g,pg(i)+k213/2);
      k413 = dt*f13_emb(ke_rib,lp.g,rib_g(i),kdp.g,pg(i)+k313);
      pg(i+1) = pg(i) +1/6*(k113 + 2*k213 + 2*k313 + k413);
      pg(i+1) = pg(i+1) + n_g*bound_inc; 

      k113 = dt*f13_emb(ke_rib,lp.l,rib_l(i),kdp.l,pl(i));
      k213 = dt*f13_emb(ke_rib,lp.l,rib_l(i),kdp.l,pl(i)+k113/2);
      k313 = dt*f13_emb(ke_rib,lp.l,rib_l(i),kdp.l,pl(i)+k213/2);
      k413 = dt*f13_emb(ke_rib,lp.l,rib_l(i),kdp.l,pl(i)+k313);
      pl(i+1) = pl(i) +1/6*(k113 + 2*k213 + 2*k313 + k413);  
      pl(i+1) = pl(i+1) + n_l*bound_inc; 

      k116 = dt*f16_emb(kd_nt,nt(i));
      k216 = dt*f16_emb(kd_nt,nt(i)+k116/2);
      k316 = dt*f16_emb(kd_nt,nt(i)+k216/2);
      k416 = dt*f16_emb(kd_nt,nt(i)+k316);
      nt(i+1) = nt(i) +1/6*(k116 + 2*k216 + 2*k316 + k416);

      if nt(i+1) < 0
         nt(i+1) = 0.;
      end

      k114 = dt*f14_emb(ke_pol,lm.t,pol_r_t(i));
      k214 = dt*f14_emb(ke_pol,lm.t,pol_r_t(i));
      k314 = dt*f14_emb(ke_pol,lm.t,pol_r_t(i));
      k414 = dt*f14_emb(ke_pol,lm.t,pol_r_t(i));
      nt_new(i+1) = nt_new(i) +1/6*(k114 + 2*k214 + 2*k314 + k414);

      nt_new(i+1) = nt_new(i+1)+nt(i+1);

      k115 = dt*f15_emb(ke_rib,lp.n,rib_n(i),kdp.n,pn(i),kd_nc,n_n,p_nc(i),...
         n_nc(i),cap(i),ke_pol,lm.tr,poln_tr(i),lm.le,polp_lec(i));
      k215 = dt*f15_emb(ke_rib,lp.n,rib_n(i),kdp.n,pn(i)+k115/2,kd_nc,n_n,p_nc(i)+k110/2,...
         n_nc(i)+k111/2,cap(i),ke_pol,lm.tr,poln_tr(i)+k17/2,lm.le,polp_lec(i)+k19_2/2);
      k315 = dt*f15_emb(ke_rib,lp.n,rib_n(i),kdp.n,pn(i)+k215/2,kd_nc,n_n,p_nc(i)+k210/2,...
         n_nc(i)+k211/2,cap(i),ke_pol,lm.tr,poln_tr(i)+k27/2,lm.le,polp_lec(i)+k29_2/2);
      k415 = dt*f15_emb(ke_rib,lp.n,rib_n(i),kdp.n,pn(i)+k315,kd_nc,n_n,p_nc(i)+k310,...
         n_nc(i)+k311,cap(i),ke_pol,lm.tr,poln_tr(i)+k37,lm.le,polp_lec(i)+k39_2);
      pn(i+1) = pn(i) +1/6*(k115 + 2*k215 + 2*k315 + k415);
      pn(i+1) = pn(i+1) + n_n*bound_inc*0; 

      if pn(i+1) < 0
          pn(i+1) = 0;
      end
      if pp(i+1) < 0
          pp(i+1) = 0;
      end
      if pm(i+1) < 0
          pm(i+1) = 0;
      end
      if pg(i+1) < 0
          pg(i+1) = 0;
      end
      if pl(i+1) < 0
          pl(i+1) = 0;
      end

      if (1-pn(i+1)/n_n/nt_new(i+1)) > 0
         n_cap(i+1) = (1-pn(i+1)/n_n/nt_new(i+1));
      else
         n_cap(i+1) = 0.;
      end

      cap(i+1) = 1 - n_cap(i+1);
      nt(i+1) = nt_new(i+1)*n_cap(i+1);
      nt_new(i+1) = 0.0;

      if pp(i+1)/pl(i+1) > p_lratio
        pol_f(i+1) = 1.0;
      else
        pol_f(i+1) = pp(i+1)/pl(i+1)/p_lratio;
      end

      k117 = dt*f17_emb(ke_pol,pol_le(i),lm.le);
      k217 = dt*f17_emb(ke_pol,pol_le(i)+k117/2,lm.le);
      k317 = dt*f17_emb(ke_pol,pol_le(i)+k217/2,lm.le);
      k417 = dt*f17_emb(ke_pol,pol_le(i)+k317,lm.le);
      pol_le(i+1) = pol_le(i) +1/6*(k117 + 2*k217+ 2*k317 + k417);

      k118 = dt*f18_emb(ke_pol,pol_trc(i),lm.tr);
      k218 = dt*f18_emb(ke_pol,pol_trc(i)+k118/2,lm.tr);
      k318 = dt*f18_emb(ke_pol,pol_trc(i)+k218/2,lm.tr);
      k418 = dt*f18_emb(ke_pol,pol_trc(i)+k318,lm.tr);
      pol_trc(i+1) = pol_trc(i) +1/6*(k118 + 2*k218+ 2*k318 + k418);

      a11 = n_nc(i+1);
      a12 = spol/lm.le;
      a13 = scond*(1.-condm)*pm(i+1);
      a14 = pl(i+1)*pol_f(i+1);
      a15 = pol_le(i+1);
      a21 = lm.le/spol;
      a22 = p_nc(i+1)*lm.tr/spol*sprom;
      a23 = pol_trc(i+1);
      a24 = sprom/sprom;
      a25 = pl(i+1)*pol_f(i+1);
      a31 = lm.le/spol;
      a32 = lm.le/spol;
      a33 = p_nc(i+1)*lm.tr/spol*sprom;
      a34 = pol_trc(i+1);
      a35 = sprom/sprom;
      a53 = pol_trc(i+1);
      a54 = pol_le(i+1);
      a55 = poln_tr(i+1);
      a56 = polp_lec(i+1);

      a51 = 0;
      for j = 1:nsec.t
         a51 = a51+polp(j,i+1);
      end
      a51 = a51+polp_lec(i+1);

      a52 = 0;
      for j = 1:nsec.n
        a52 = a52 + poln_n(j,i+1) ;
      end
      for j = 1:nsec.p
        a52 = a52 + poln_p(j,i+1) ;
      end
      for j = 1:nsec.m
        a52 = a52 + poln_m(j,i+1) ;
      end
      for j = 1:nsec.g
        a52 = a52 + poln_g(j,i+1) ;
      end
      for j = 1:nsec.l
        a52 = a52 + poln_l(j,i+1) ;
      end
      a52 = a52 + poln_tr(i+1);

      a61 = pol_tn(i+1)+pol_tp(i+1)+pol_tm(i+1)+pol_tg(i+1)+pol_tl(i+1);

      x(1,1) = n_nc_pol(i);
      x(2,1) = d_pol_le(i);
      x(3,1) = pol_t(i);
      x(4,1) = d_pol_term(i);
      x(5,1) = d_pol_term(i)-d_pol_le(i);
      x(6,1) = pol_r_t(i);
      xx = x;
      crite = 100;  
      criteo = 1;
      k = 1;
      while (crite > 0.001) && (k < 100)

      k = k+1;
      A = zeros(6);
      A(1,1) = 1.;
      A(1,2) = a12/(a13/(a14-xx(3))+1);
      A(1,3) = (a11-a12*(a15+xx(2)))*a13/(a14-xx(3)+a13)^2.;
      if (xx(1)*a21+(a22-(a23+xx(5)))*a24) < (a25-xx(3))
       A(2,1) = -1*a21;
       A(2,4) = 1;
       A(2,5) = a24;
      else
       A(2,3) = 1;
       A(2,4) = 1;
      end
      A(3,1) = -1*a31*a35*xx(4)*(a33-(a34+xx(5)))/(a32*xx(1)+a35*(a33-(a34+xx(5))))^2.;
      A(3,2) = 1;
      A(3,4) = -1*a31*xx(1)/(a32*xx(1)+a35*(a33-(a34+xx(5))));
      A(3,5) = -1*a31*xx(1)*xx(4)*a35/(a32*xx(1)+a35*(a33-(a34+xx(5))))^2.;
      A(4,2) = 1;
      A(4,4) = -1;
      A(4,5) = 1;
      A(5,2) = -1;
      A(5,5) = -1;
      A(5,6) = 1;
      A(6,3) =1;
      A(6,6) = -1;

      f(1) = xx(1)-(a11-a12*(a15+xx(2)))/(a13/(a14-xx(3))+1);
      if (xx(1)*a21+a24*(a22-(a23+xx(5)))) < (a25-xx(3))
       f(2) = xx(4)-(xx(1)*a21+a24*(a22-(a23+xx(5))));
      else
       f(2) = xx(4)-a25+xx(3);
      end
      f(3) = xx(2) - a31*xx(1)/(a32*xx(1)+a35*(a33-(a34+xx(5))))*xx(4);
      f(4) = xx(5)-xx(4)+xx(2);
      f(5) = xx(6)-(a53+xx(5))-a51-(a54+xx(2))-a52-a55-a56;
      f(6) = xx(3)-xx(6)-a61;
      ff = f';
      delt = -A\ff;

      criten = sum(delt);
      crite = abs(criteo-criten)/abs(criteo);
      criteo = criten;
      xx = xx + te*delt;
      end

      n_nc_pol(i+1) = xx(1);
      d_pol_le(i+1) = xx(2);
      pol_t(i+1) = xx(3);
      d_pol_term(i+1) = xx(4);
      d_pol_trc(i+1) = xx(5);
      pol_r_t(i+1) = xx(6);

      pol_le(i+1) = pol_le(i+1)+d_pol_le(i+1);
      pol_trc(i+1) = pol_trc(i+1)+d_pol_trc(i+1);

      
%       if (n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1))/n_nc(i+1) < -0.001
%           test = abs((n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1))/n_nc(i+1))
%       end

        if (n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1)) < 0
            delta_cor = (n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1));
            n_nc_pol(i+1) = n_nc_pol(i+1)+1.0001*delta_cor*(n_nc_pol(i+1)/(n_nc_pol(i+1)+pol_le(i+1)*spol/lm.le));
            pol_le(i+1) = pol_le(i+1)+1.0001*delta_cor*(pol_le(i+1)*spol/lm.le/(n_nc_pol(i+1)+pol_le(i+1)*spol/lm.le))...
              *lm.le/spol;
            relm(i+1) = 0;
        else


            if (pm(i+1)*(1-condm)/n_m/(n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1))) < 1
            relm(i+1) = (pm(i+1)*(1-condm)/n_m/(n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1)));
            else
            relm(i+1) = 1;
            end 
        end

      n_nc_m(i+1) = relm(i+1)*(n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1))+n_nc_m(i);
      n_nc(i+1) = n_nc(i+1)-relm(i+1)*(n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1));


    if (pg(i+1)/n_g/n_nc_m(i+1)) < 1
        rel(i+1) = (pg(i+1)/n_g/n_nc_m(i+1));
    else
        rel(i+1) = 1;
    end

        proge(i+1) = proge(i)+rel(i+1)*n_nc_m(i+1);
        progen(i+1) = (proge(i)+rel(i+1)*n_nc_m(i+1))*navo*vc;

        n_nc_m(i+1) = n_nc_m(i+1) - rel(i+1)*n_nc_m(i+1);


        pp(i+1) = pp(i+1) - n_p*rel(i+1)*n_nc_m(i+1);
        pm(i+1) = pm(i+1) - n_m*relm(i+1)*(n_nc(i+1)-pol_le(i+1)*spol/lm.le-n_nc_pol(i+1));
        pg(i+1) = pg(i+1) - n_g*rel(i+1)*n_nc_m(i+1);
        pl(i+1) = pl(i+1) - n_l*rel(i+1)*n_nc_m(i+1);



        % Analysis routine %
    %     n_p_ratio(i+1) = n_nc(i+1)/p_nc(i+1);
    %     r_virus(i) = (proge(i+1)-proge(i))*navo*vc/dt;
    %     r_n_nc(i) = (n_nc(i+1)-n_nc(i))*navo*vc/dt;
    %     r_p_nc(i) = (p_nc(i+1)-p_nc(i))*navo*vc/dt;
    %     r_pn(i) = (pn(i+1)-pn(i))*navo*vc/dt;
    %     r_pp(i) = (pp(i+1)-pp(i))*navo*vc/dt;
    %     r_pm(i) = (pm(i+1)-pm(i))*navo*vc/dt;
    %     r_pg(i) = (pg(i+1)-pg(i))*navo*vc/dt;
    %     r_pl(i) = (pl(i+1)-pl(i))*navo*vc/dt;
    %     rr_pn(i) = (pn(i+1)-pn(i)+n_n*(n_nc(i+1)+p_nc(i+1)+n_nc_m(i+1)+proge(i+1)-...
    %      n_nc(i)-p_nc(i)-n_nc_m(i)-proge(i)))*navo*vc/dt;
    %     rr_pp(i) = (pp(i+1)-pp(i)+n_p*(proge(i+1)-proge(i)))*navo*vc/dt;
    %     rr_pm(i) = (pm(i+1)-pm(i)+n_m*(n_nc_m(i+1)+proge(i+1)-...
    %      n_nc_m(i)-proge(i)))*navo*vc/dt;
    %     rr_pg(i) = (pg(i+1)-pg(i)+n_g*(proge(i+1)-proge(i)))*navo*vc/dt;
    %     rr_pl(i) = (pl(i+1)-pl(i)+n_l*(proge(i+1)-proge(i)))*navo*vc/dt;
    %     r_mrna_n(i) = (mrna_n(i+1)-mrna_n(i))*navo*vc/dt;
    %     r_mrna_p(i) = (mrna_p(i+1)-mrna_p(i))*navo*vc/dt;
    %     r_mrna_m(i) = (mrna_m(i+1)-mrna_m(i))*navo*vc/dt;
    %     r_mrna_g(i) = (mrna_g(i+1)-mrna_g(i))*navo*vc/dt;
    %     r_mrna_l(i) = (mrna_l(i+1)-mrna_l(i))*navo*vc/dt;
    %     ratioo(i+1) = (proge(i+1)+n_nc(i+1))/p_nc(i+1);

    end


    tt = (1:q+1)*dt/3600;


    if ~NO_PLOTS
        figure
        plot(tt,progen)
        xlabel('Time(hr)')
        ylabel('# of Progeny/Cell')

        uo = progen(q+1)

        figure
        plot(tt,mrna_n,'-k')
        hold on
        plot(tt,mrna_p,'-b')
        plot(tt,mrna_m,'-g')
        plot(tt,mrna_g,'-r')
        plot(tt,mrna_l,'--k')
        %axis([0 25 0 2.5E5])
        xlabel('Time, hr')
        ylabel('mRNA production (#/cell)')

        pro_index = find(progen >= 1,1);
        figure
        semilogy(tt(pro_index:end),progen(pro_index:end))
        axis([0 25 1E-3 1E4])
        xlabel('Time, hr')
        ylabel('Virion production (#/cell)')

    end


    %% Package them up
    output = struct();

    varbls = who;
    %only output state variables i.e. that have a recorded value for every time
    %point

    for i = 1:length(varbls)
        if length(eval(varbls{i})) > q
            output.(varbls{i}) = eval(varbls{i});
        end
    end

    output.uo = progen(q+1);



end


%% %%%%%%%%%%%%%%%%%%%%%%%

function w1 = f1_emb(ke_pol, eta_n, n_cap, pol_le, l_le, nsec_n, lm_n, pol_n)

    w1 = ke_pol*(eta_n*n_cap*pol_le/l_le - nsec_n/lm_n*pol_n);

    %Conc of polymerase at gene N, segment 1
    %ke_pol = elongation rate of polymerase (nt/s)
end

function w2 = f2_emb(ke_pol, nsec_n, lm_n, pol_n1, pol_n2)
    
    w2 = ke_pol*nsec_n/lm_n*(pol_n1 - pol_n2);
    
end

function w3 = f3_emb(ke_pol, eta_p, nsec_n, lm_n, pol_n, nsec_p, lm_p, pol_p)

    w3 = ke_pol*(eta_p*nsec_n/lm_n*pol_n - nsec_p/lm_p*pol_p);

end

function w4 = f4_emb(ke_pol, cap, pol_le, l_le, nsec_n, lm_n, poln_n)

    w4 = ke_pol*(cap*pol_le/l_le - nsec_n/lm_n*poln_n); 
    
end

function w5 = f5_emb(ke_pol, nsec_n, lm_n, poln_n1, poln_n2)

    w5 = ke_pol*nsec_n/lm_n*(poln_n1 - poln_n2);

end

function w6 = f6_emb(ke_pol, cap, nsec_p, lm_p, poln_p, nsec_m, lm_m, poln_m)

    w6 = ke_pol*(cap*nsec_p/lm_p*poln_p - nsec_m/lm_m*poln_m);    
    
end

function w7 = f7_emb(ke_pol, nsec_l, lm_l, poln_l, poln_tr, l_tr)

    w7 = ke_pol*(nsec_l/lm_l*poln_l - poln_tr/l_tr);    

end

function w8 = f8_emb(ke_pol, cap, pol_trc, l_tr, nsec_t, l_t, l_le, polp)

    w8 = ke_pol*(cap*pol_trc/l_tr - nsec_t/(l_t - l_tr - l_le)*polp);

end

function w9_2 = f9_2_emb(ke_pol, nsec_t, l_t, l_tr, l_le, polp, pol_lec)

    w9_2 = ke_pol*(nsec_t/(l_t - l_tr - l_le)*polp - pol_lec/l_le);

end

function w9 = f9_emb(ke_pol, nsec_t, l_t, l_tr, l_le, polp1, polp2)

    w9 = ke_pol*nsec_t/(l_t - l_tr - l_le)*(polp1 - polp2);

end

function w10 = f10_emb(ke_pol, cap, poln_tr, l_tr, kd_nc, p_nc)

    w10 = (ke_pol*cap*poln_tr/l_tr - kd_nc*p_nc);

end

function w11 = f11_emb(ke_pol, cap, polp, l_le, kd_nc, n_nc)

    w11 = (ke_pol*cap*polp/l_le - kd_nc*n_nc);     
    
end

function w12 = f12_emb(ke_pol, lm_n, pol_tn, kd_m, mrna_n)

    w12 = (ke_pol/lm_n*pol_tn - kd_m*mrna_n);  

end

function w13 = f13_emb(ke_rib, lp_p, rib_p, kdp_p, pp)

    w13 = (ke_rib/lp_p*rib_p - kdp_p*pp);
 
end

function w14 = f14_emb(ke_pol, l_t, pol_r_t)

    w14 = (ke_pol/l_t*pol_r_t);

end

function w15 = f15_emb(ke_rib, lp_n, rib_n, kdp_n, pn, kd_nc, n_n, p_nc, ...
    n_nc, cap, ke_pol, l_tr, poln_tr, l_le, polp_lec)

    w15 = ke_rib/lp_n*rib_n - kdp_n*pn + kd_nc*n_n*(p_nc + n_nc)...
        - cap*ke_pol*n_n*(1/l_tr*poln_tr + polp_lec/l_le);
   
end

function w16 = f16_emb(kd_nt, nt)

    w16 = -1*kd_nt*nt;
 
end

function w17 = f17_emb(ke_pol, pol_le, l_le)

    w17 = -ke_pol*pol_le/l_le;
     
end

function w18 = f18_emb(ke_pol, pol_trc, l_tr)

    w18 = -ke_pol*pol_trc/l_tr;
  
end

function w_ent = f_ent_emb(k_int, v_ex)

    w_ent = -k_int*v_ex;

end


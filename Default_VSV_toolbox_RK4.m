function pars = Default_VSV_toolbox_RK4()


    %% setup pars with default parameters
    % ie this would yield same result on toolbox if pars = ()
    pars = struct();

    pars.kdp = struct();
    pars.kdp.n = 3.45*1e-5;   
    pars.kdp.p = 1.4*1e-6;   
    pars.kdp.m = 1.5*1e-4;   
    pars.kdp.g = 5.7*1e-5;    
    pars.kdp.l1 = 1.2*1e-5;  
    pars.kdp.l2 = 4.3*1e-5; 

    pars.eta = struct();
    pars.eta.n = 1;
    pars.eta.p = 0.75;
    pars.eta.m = 0.75;
    pars.eta.g = 0.75;
    pars.eta.l = 0.05;

    pars.lm = struct();
    pars.lm.le = 47;          %nt, length of leader sequence
    pars.lm.lec = pars.lm.le;       %nt, length of complementary leader region of VSV genome unused !!!!

    pars.lm.tr = 59;          %nt, length of trailer, this partitioning will be used to create delays in the transcription and replication of the genome
    pars.lm.trc = pars.lm.tr;       %nt, length of anti-genomic promoter
    pars.lm.t = 11150;        %nt, sum of all parts. we are missing 11 nt from the total, who knows why?

    pars.lm.n = 1326;        %nt, length of N mRNA
    pars.lm.p = 814;         %nt, length of P mRNA
    pars.lm.m = 831;         %nt, length of M mRNA
    pars.lm.g = 1665;        %nt, length of G mRNA
    pars.lm.l = 6373;        %nt, length of L mRNA

    pars.lp = struct();
    pars.lp.n = 422;  %aa, length of protein N
    pars.lp.p = 265;  %aa, length of protein P
    pars.lp.m = 229;  %aa, length of protein M
    pars.lp.g = 511;  %aa, length of protein G
    pars.lp.l = 2109; %aa, length of protein L    

    pars.nsec = struct();
    pars.nsec.le = 1; %count, number of segments in gene leader
    pars.nsec.tr = 1; %count, number of segments in gene trailer
    pars.nsec.t = 40; %count, number of segments in genome

    pars.nsec.n = 5;  %count, number of segments in gene N
    pars.nsec.p = 3;  %count, number of segments in gene P
    pars.nsec.m = 3;  %count, number of segments in gene M
    pars.nsec.g = 6;  %count, number of segments in gene G
    pars.nsec.l = 23; %count, number of segments in gene L

    pars.kd_m = 1.9*1e-4; 
    pars.kd_nc = 1.9*1e-5; 
    pars.kd_nt = 1.9*1e-4;   
    pars.ke_pol = 3.7;      
    pars.ke_rib = 6.0*1.0;  
    pars.p_lratio = 3.6; 
    pars.ini_a = 0.99; 
    pars.sprom = 5.43;
    pars.spol = 170;
    pars.srib = 238.5; 
    pars.scond = 2.77*1e-4*0.2*1.125;     
    pars.condm = 0.1;
    pars.k_int = 0.7194/3600;

end





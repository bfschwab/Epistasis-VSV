t = 0:.1:25;
maxx = 2500;
W_poor = zeros(1,length(t));
W_poor(61:201) = maxx/14*(t(61:201)-6);
W_poor(202:end) = maxx;

W_rich = W_poor.^(1./t);

xsize = 30;
ysize = 25;
width = 3;
ticksize = 18;

%%
%% setup pars with default parameters
% ie this would yield same result on toolbox if pars = ()

defpars = Default_VSV_toolbox_RK4();

NO_PLOTS = 1;
tlength = 25*3600; %25 hr
dt = 15; %seconds per step

% default results
output = VSV_toolbox_RK4(defpars, tlength, dt, NO_PLOTS);
t_wt = 0:(25/6000):25;
W_poor_wt = output.progen;
W_rich_wt = output.progen.^(1./t_wt);

%%
figure
subplot(2,1,1)
plot(t, W_poor,'--k','LineWidth',width)
hold on
plot(t_wt, W_poor_wt,'k','LineWidth',width)
axis([0 25 0 3000])
ax = gca; % current axes
ax.FontSize = ticksize;
ax.TickDir = 'both';
xlabel('$t$','Interpreter','latex','FontSize',xsize)
ylabel('$N(t)$','Interpreter','latex','FontSize',ysize)


subplot(2,1,2)
plot(t, W_rich,'--k','LineWidth',width)
hold on
plot(t_wt, W_rich_wt,'k','LineWidth',width)
axis([0 25 0 3])
ax = gca; % current axes
ax.FontSize = ticksize;
ax.TickDir = 'both';
xlabel('$t$','Interpreter','latex','FontSize',xsize)
ylabel('$[N(t)]^{1/t}$','Interpreter','latex','FontSize',ysize)


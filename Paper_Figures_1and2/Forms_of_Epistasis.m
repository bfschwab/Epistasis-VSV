n = linspace(0,11,111);

alpha1 = .00418;
beta1 = 1;

alpha2 = .0001;
beta2 = 2.75;

alpha3 = .015;
beta3 = .4;

y1 = -alpha1*n.^beta1;
y2 = -alpha2*n.^beta2;
y3 = -alpha3*n.^beta3;

xsize = 30;
ysize = 30;
lsize = 15;
width = 3;
ticksize = 18;

figure
plot(n,y1,'-k','LineWidth',width)
hold on
plot(n,y2,':k','LineWidth',width)
plot(n,y3,'--k','LineWidth',width)
axis([0 11 -.1 0])
ax = gca; % current axes
ax.FontSize = ticksize;
ax.TickDir = 'both';
xlabel('$n$','Interpreter','latex','FontSize',xsize)
ylabel('log($w$)','Interpreter','latex','FontSize',ysize)
legend('multiplicative','synergistic','antagonistic',...
    'FontSize',lsize,'Location','southwest')
yticks(-.1:.02:0)
xticks(0:2:10)
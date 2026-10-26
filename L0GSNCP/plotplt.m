%% One figure for each data set: objective values and certificate residual.
function plt0=plotplt(loss,timerun,R,name)
plt0=figure('Color','w','Name',name,'NumberTitle','off','Position',[100,100,1500,650]);
subplot(1,2,1);
semilogy(timerun,loss(:,1),'-','LineWidth',2.6);
hold on;
semilogy(timerun,loss(:,2),'--','LineWidth',2.6);
prettyAxes().gbase2();
set(gca,'FontSize',25);
xlabel('Time (seconds)','FontSize',25);
ylabel('Objective function value','FontSize',25);
h=legend({'J','J+C'},'Location','best');
set(h,'FontSize',25);
title(name,'Interpreter','none');

subplot(1,2,2);
semilogy(timerun,R,'-','LineWidth',2.6);
prettyAxes().gbase2();
set(gca,'FontSize',25);
xlabel('Time (seconds)','FontSize',25);
ylabel('Certificate residual','FontSize',25);
title(name,'Interpreter','none');
end

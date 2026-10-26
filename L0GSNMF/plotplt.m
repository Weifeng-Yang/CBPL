function plt=plotplt(lossdata,trdata,residual,index)
names={'lp_ship12l','BASEHOCK'};
plt=figure('Name',names{index},'NumberTitle','off');
subplot(1,2,1);
semilogy(trdata,lossdata(:,1),'-','LineWidth',2);
hold on;
semilogy(trdata,lossdata(:,2),'--','LineWidth',2);
prettyAxes().gbase2();
legend({'$J$','$J+C$'},'Interpreter','latex');
xlabel('Time (seconds)');
ylabel('Objective function value');
title(names{index},'Interpreter','none');
subplot(1,2,2);
semilogy(trdata,residual,'-','LineWidth',2);
prettyAxes().gbase2();
xlabel('Time (seconds)');
ylabel('Certificate residual');
title(names{index},'Interpreter','none');
end

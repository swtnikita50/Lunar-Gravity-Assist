% Created on: 14/09/22 18:40
% Function: Plot Orbit from keplerian orbital elements with initial and
% final true anomaly

function plotOrbit(kepOE, TAinit, TAfinal, mu_E, nNodes, dispName,clr)
TA = TAinit:(TAfinal-TAinit)/(nNodes-1):TAfinal;
currOE(1:5) = kepOE(1:5);

for i = 1:nNodes
    currOE(6) = TA(i);
    [rr(i,:),~] = po2pv(currOE,mu_E);
end

plot3(rr(:,1),rr(:,2),rr(:,3),'DisplayName',dispName,'Color',clr,'LineWidth',2);hold on;grid on;
% scatter3(rr(1,1),rr(1,2),rr(1,3),clr);
% scatter3(rr(end,1),rr(end,2),rr(end,3),clr);



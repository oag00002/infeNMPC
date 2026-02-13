close all;
clear all;
term = importdata('term_check_nom.tab') ; 
simlen=98;
figure 
buffer=15;
plot(0:simlen,term(2:simlen+2,2),0:simlen,term(1:simlen+1,3),'--',0:simlen,term(1:simlen+1,3)+buffer)
axis([0 100 0 40])
legend({'NMPC $S_T$','asNMPC $S_T$','Recommended N'},'Location','Northwest','interpreter','latex')
xlabel('k')
ylabel('i')

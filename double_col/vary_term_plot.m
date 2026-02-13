close all;
clear all;

figure
r2n10 = importdata('r2n10.tab');
r100n10 = importdata('r100n10.tab');
r1000n10 = importdata('r1000n10.tab');
adaptiven10 = importdata('adaptiven10.tab');
plot(adaptiven10(:,1),adaptiven10(:,2),r2n10(:,1),r2n10(:,2),'--',r100n10(:,1),r100n10(:,2),r1000n10(:,1),r1000n10(:,2))
legend('adaptive','\rho_L=\rho_\Psi=2','100','1000')
xlabel('k')
ylabel('$|x_k|$')

r1p1 = importdata('rho1p1.tab');
r5 = importdata('rho5.tab');
r10 = importdata('rho10.tab');
r50 = importdata('rho50.tab');
r200 = importdata('rho200.tab');
r1000 = importdata('rho1000.tab');
r10000 = importdata('rho10000.tab') ; 
r100000 = importdata('rho100000.tab') ; 

figure
plot(r1p1(:,1),r1p1(:,2),r5(:,1),r5(:,2),r10(:,1),r10(:,2),r50(:,1),r50(:,2),r200(:,1),r200(:,2),r1000(:,1),r1000(:,2),r10000(:,1),r10000(:,2),r100000(:,1),r100000(:,2))
legend('$\rho=1.1$','$\rho=5$','$\rho=10$','$\rho=50$','$\rho=200$','1000','10000','100000')
xlabel('k')
ylabel('$||x(k)||$')

n10rho10=importdata('n10rho10.tab');
n10rho100=importdata('n10rho100.tab');
n10rho1000=importdata('n10rho1000.tab');
n10rho10000=importdata('n10rho10000.tab');
n10rho100000=importdata('n10rho100000.tab');

figure
plot(n10rho10(:,1),n10rho10(:,2),n10rho100(:,1),n10rho100(:,2),'--',n10rho1000(:,1),n10rho1000(:,2),'-.',n10rho10000(:,1),n10rho10000(:,2),':',n10rho100000(:,1),n10rho100000(:,2))
legend('\rho_\psi=10','\rho_\psi=100','\rho_\psi=1000','\rho_\psi=10000','\rho_\psi=100000')
xlabel('k')
ylabel('$||x(k)||$')
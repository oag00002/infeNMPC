clear all;
close all;
K=149;

qiahp1=importdata('qiahnmpc_p1_nom_50.tab');
qiahp2=importdata('qiahnmpc_p2_nom_50.tab');
qiahp3=importdata('qiahnmpc_p3_nom_50.tab');
nmpc=importdata('standard_nom_50.tab');

ussp1=[.5621708 ,.8478292,4.296642,4.858812 ,.2839724,.5638568,2.410230,2.694203];
ussp2=[ 0.553094, 0.856906, 5.29703, 5.85012, 0.580305, 0.276602, 2.19216, 2.77246];
ussp3=[ 0.267924, 1.14208, 3.99199, 4.25992, 0.582946, 0.55913, 3.07372,3.65667] ;

xasetp1=.9839 ;
xasetp2= .9806;
xasetp3= .9749;
xbsetp1= .95;
xbsetp2= .95; 
xbsetp3 = .95;
xcsetp1= .9943;
xcsetp2= .9929;
xcsetp3= .9938;

qiah_states=[qiahp1(:,2);qiahp2(:,2);qiahp3(:,2)];
qiah_N=[qiahp1(:,5);qiahp2(:,5);qiahp3(:,5)];
qiah_time=[qiahp1(:,6);qiahp2(:,6);qiahp3(:,6)];
term_step=[qiahp1(:,3);qiahp2(:,3);qiahp3(:,3)];
sens_term_step=[qiahp1(:,4);qiahp2(:,4);qiahp3(:,4)];
nmpc_xa = nmpc(:,5);
nmpc_xb = nmpc(:,6) ; 
nmpc_xc = nmpc(:,7) ;
nmpc_u=nmpc(:,254:261);

for i = 1:50
nmpc_u_norm(i) = norm(nmpc_u(i,:) - ussp1) ;  
end
for i = 51:100
nmpc_u_norm(i) = norm(nmpc_u(i,:) - ussp2) ;  
end
for i = 101:150
nmpc_u_norm(i) = norm(nmpc_u(i,:) - ussp3) ;  
end

qiah_xa = [qiahp1(:,7);qiahp2(:,7);qiahp3(:,7)] ;
qiah_xb = [qiahp1(:,8);qiahp2(:,8);qiahp3(:,8)] ;
qiah_xc = [qiahp1(:,9);qiahp2(:,9);qiahp3(:,9)] ;
qiah_u = [qiahp1(:,256:263);qiahp2(:,256:263);qiahp3(:,256:263)];

for i = 1:50
qiah_u_norm(i) = norm(qiah_u(i,:) - ussp1) ;  
end
for i = 51:100
qiah_u_norm(i) = norm(qiah_u(i,:) - ussp2) ;  
end
for i = 101:150
qiah_u_norm(i) = norm(qiah_u(i,:) - ussp3) ;  
end

states=nmpc(:,2) ;
time=nmpc(:,4) ; 

time(1)=104.684 ; 
qiah_time(1) = time(1);

figure
hold on 
plot(0:K,states,0:K,qiah_states,'--')
axis([0 150 0 3]) 
stairs([0,50,100,150],[.41,.36,.38,.38],'Linewidth',1.75)
ylabel('$||x_k-x_{ss}||$')
xlabel('$k$')
legend('NMPC','AH-NMPC','Terminal Region')
box on
print('nom_states','-dpdf')

figure
hold on 
plot(0:K,nmpc_u_norm,0:K,qiah_u_norm,'--')
axis([0 150 0 4]) 
ylabel('$||u_k-u_{ss}||$')
xlabel('$k$')
legend('NMPC','AH-NMPC')
box on
print('nom_cont','-dpdf')


figure
plot(0:K-1,25*ones(1,K),0:K,qiah_N)
axis([0 150 0 35]) 
ylabel('$N_k$')
xlabel('$k$')
legend('NMPC','AH-NMPC')

figure
plot(0:K,time,0:K,qiah_time)
axis([0 150 0 200]) 
ylabel('CPU time (s)')
xlabel('$k$')
legend('NMPC','AH-NMPC')

figure
hold on
plot(0:K,nmpc_xa,0:K,nmpc_xb,0:K,nmpc_xc)
stairs([0,50,100,150],[xasetp1,xasetp2,xasetp3,xasetp3],'Linewidth',1.75,'Color',[ 0    0.4470    0.7410],'LineStyle','--')
stairs([0,50,100,150],[xbsetp1,xbsetp2,xbsetp3,xbsetp3],'Linewidth',1.75,'Color',[0.8500    0.3250    0.0980],'LineStyle','--')
stairs([0,50,100,150],[xcsetp1,xcsetp2,xcsetp3,xcsetp3],'Linewidth',1.75,'Color',[    0.9290    0.6940    0.1250],'LineStyle','--')
legend('xa','xb','xc','Location','southeast')
xlabel('k')
ylabel('$x_{out}$')
%get(gca,'colororder')
axis([0 150 .95 1])

figure
hold on
plot(0:K,nmpc_u(:,1),0:K,nmpc_u(:,5),0:K,nmpc_u(:,6))
stairs([0,50,100,150],[ussp1(1),ussp2(1),ussp3(1),ussp3(1)],'Linewidth',1.75,'Color',[ 0    0.4470    0.7410],'LineStyle','--')
stairs([0,50,100,150],[ussp1(5),ussp2(5),ussp3(5),ussp3(5)],'Linewidth',1.75,'Color',[0.8500    0.3250    0.0980],'LineStyle','--')
stairs([0,50,100,150],[ussp1(6),ussp2(6),ussp3(6),ussp3(6)],'Linewidth',1.75,'Color',[    0.9290    0.6940    0.1250],'LineStyle','--')
legend('D1','D2','B2','Location','northeast')
xlabel('k')
ylabel('$u_{out}$')
%get(gca,'colororder')
axis([0 150 0 2])
%breakyaxis([2 5])

figure
hold on
plot(0:K,qiah_xa,0:K,qiah_xb,0:K,qiah_xc)
stairs([0,50,100,150],[xasetp1,xasetp2,xasetp3,xasetp3],'Linewidth',1.75,'Color',[ 0    0.4470    0.7410],'LineStyle','--')
stairs([0,50,100,150],[xbsetp1,xbsetp2,xbsetp3,xbsetp3],'Linewidth',1.75,'Color',[0.8500    0.3250    0.0980],'LineStyle','--')
stairs([0,50,100,150],[xcsetp1,xcsetp2,xcsetp3,xcsetp3],'Linewidth',1.75,'Color',[    0.9290    0.6940    0.1250],'LineStyle','--')
legend('xa','xb','xc','Location','southeast')
xlabel('k')
ylabel('$x_{out}$')
axis([0 150 .95 1])
xlabel('sens')

figure
hold on
plot(0:K,qiah_u(:,1),0:K,qiah_u(:,5),0:K,qiah_u(:,6))
stairs([0,50,100,150],[ussp1(1),ussp2(1),ussp3(1),ussp3(1)],'Linewidth',1.75,'Color',[ 0    0.4470    0.7410],'LineStyle','--')
stairs([0,50,100,150],[ussp1(5),ussp2(5),ussp3(5),ussp3(5)],'Linewidth',1.75,'Color',[0.8500    0.3250    0.0980],'LineStyle','--')
stairs([0,50,100,150],[ussp1(6),ussp2(6),ussp3(6),ussp3(6)],'Linewidth',1.75,'Color',[    0.9290    0.6940    0.1250],'LineStyle','--')
legend('D1','D2','B2','Location','northeast')
xlabel('k')
ylabel('$u_{out}$')
%get(gca,'colororder')
axis([0 150 0 2])
%breakyaxis([2 5])


close all;
ind=2
time = 1 : 81; 

economic=importdata('economic.tab') ;
economic = economic(time,ind) ; 

tracking=importdata('tracking.tab') ;
tracking = tracking(time,ind) ; 

gersh100 = importdata('gersh_100.tab') ;
gersh100 = gersh100(time,ind);

sigmap01 = importdata('desc_con_p01.tab') ;
sigmap01 = sigmap01(time,ind);

gersh100mis = importdata('gersh_100_mismatch.tab') ;
gersh100mis = gersh100mis(time,ind);

sigmap01mis = importdata('desc_con_p01_mismatch.tab') ;
sigmap01mis = sigmap01mis(time,ind);



timeplot = 0:80

figure
plot(timeplot,gersh100,time,sigmap01)
legend('100 \% Regularization','Stab con, 1 norm, $\delta$=.01','location','southeast')
ylabel('$l^{ec}(x_k,u_k)$')
xlabel('k')
%axis([0 80 -10 2])
set(findall(gcf,'type','text'),'fontSize',14)





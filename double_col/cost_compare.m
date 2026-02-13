close all;
horizon1 = 1;
horizon2=10;
horizon_unc=10;
lss=-0.223548 ; 

nominal

economic=importdata('economic.tab') ;
economic = economic(horizon1:horizon2,3) ; 
economic = sum(economic+0.223548) ;
economic

tracking=importdata('tracking.tab') ;
tracking = tracking(horizon1:horizon2,3) ; 
tracking = sum(tracking+0.223548) ;
tracking

sigmap99 = importdata('desc_con_p99.tab') ;
sigmap99 = sigmap99(horizon1:horizon2,3);
sigmap99 = sum(sigmap99+0.223548);
sigmap99 

sigmap9 = importdata('desc_con_p9.tab') ;
sigmap9 = sigmap9(horizon1:horizon2,3);
sigmap9 = sum(sigmap9+0.223548);
sigmap9 

sigmap7 = importdata('desc_con_p7.tab') ;
sigmap7 = sigmap7(horizon1:horizon2,3);
sigmap7 = sum(sigmap7+0.223548);
sigmap7 

sigmap5 = importdata('desc_con_p5.tab') ;
sigmap5 = sigmap5(horizon1:horizon2,3);
sigmap5 = sum(sigmap5+0.223548);
sigmap5 

sigmap3 = importdata('desc_con_p3.tab') ;
sigmap3 = sigmap3(horizon1:horizon2,3);
sigmap3 = sum(sigmap3+0.223548);
sigmap3 

sigmap1 = importdata('desc_con_p1.tab') ;
sigmap1 = sigmap1(horizon1:horizon2,3);
sigmap1 = sum(sigmap1+0.223548);
sigmap1 

sigmap01 = importdata('desc_con_p01.tab') ;
sigmap01 = sigmap01(horizon1:horizon2,3);
sigmap01 = sum(sigmap01+0.223548);
sigmap01 



%{
sigmap9s = importdata('desc_con_p9_seq.tab') ;
sigmap9s = sigmap9s(horizon1:horizon2,3);
sigmap9s = sum(sigmap9s+0.223548);
sigmap9s 

sigmap7s = importdata('desc_con_p7_seq.tab') ;
sigmap7s = sigmap7s(horizon1:horizon2,3);
sigmap7s = sum(sigmap7s+0.223548);
sigmap7s 

sigmap5s = importdata('desc_con_p5_seq.tab') ;
sigmap5s = sigmap5s(horizon1:horizon2,3);
sigmap5s = sum(sigmap5s+0.223548);
sigmap5s 

sigmap3s = importdata('desc_con_p3_seq.tab') ;
sigmap3s = sigmap3s(horizon1:horizon2,3);
sigmap3s = sum(sigmap3s+0.223548);
sigmap3s 

sigmap1s = importdata('desc_con_p1_seq.tab') ;
sigmap1s = sigmap1s(horizon1:horizon2,3);
sigmap1s = sum(sigmap1s+0.223548);
sigmap1s 

sigmap01s = importdata('desc_con_p01_seq.tab') ;
sigmap01s = sigmap01s(horizon1:horizon2,3);
sigmap01s = sum(sigmap01s+0.223548);
sigmap01s 
%}

nominal

sigmap99n1 = importdata('desc_con_1norm_p99.tab') ;
sigmap99n1 = sigmap99n1(horizon1:horizon2,3);
sigmap99n1 = sum(sigmap99n1+0.223548);
sigmap99n1 

sigmap9n1 = importdata('desc_con_1norm_p9.tab') ;
sigmap9n1 = sigmap9n1(horizon1:horizon2,3);
sigmap9n1 = sum(sigmap9n1+0.223548);
sigmap9n1 

sigmap7n1 = importdata('desc_con_1norm_p7.tab') ;
sigmap7n1 = sigmap7n1(horizon1:horizon2,3);
sigmap7n1 = sum(sigmap7n1+0.223548);
sigmap7n1 

sigmap5n1 = importdata('desc_con_1norm_p5.tab') ;
sigmap5n1 = sigmap5n1(horizon1:horizon2,3);
sigmap5n1 = sum(sigmap5n1+0.223548);
sigmap5n1 

sigmap3n1 = importdata('desc_con_1norm_p3.tab') ;
sigmap3n1 = sigmap3n1(horizon1:horizon2,3);
sigmap3n1 = sum(sigmap3n1+0.223548);
sigmap3n1 

sigmap1n1 = importdata('desc_con_1norm_p1.tab') ;
sigmap1n1 = sigmap1n1(horizon1:horizon2,3);
sigmap1n1 = sum(sigmap1n1+0.223548);
sigmap1n1 

sigmap01n1 = importdata('desc_con_1norm_p01.tab') ;
sigmap01n1 = sigmap01n1(horizon1:horizon2,3);
sigmap01n1 = sum(sigmap01n1+0.223548);
sigmap01n1 

sigmap99n50 = importdata('desc_con_n50_p99.tab') ;
sigmap99n50 = sigmap99n50(horizon1:horizon2,3);
sigmap99n50 = sum(sigmap99n50+0.223548);
sigmap99n50 

sigmap9n50 = importdata('desc_con_n50_p9.tab') ;
sigmap9n50 = sigmap9n50(horizon1:horizon2,3);
sigmap9n50 = sum(sigmap9n50+0.223548);
sigmap9n50 

sigmap7n50 = importdata('desc_con_n50_p7.tab') ;
sigmap7n50 = sigmap7n50(horizon1:horizon2,3);
sigmap7n50 = sum(sigmap7n50+0.223548);
sigmap7n50 

sigmap5n50 = importdata('desc_con_n50_p5.tab') ;
sigmap5n50 = sigmap5n50(horizon1:horizon2,3);
sigmap5n50 = sum(sigmap5n50+0.223548);
sigmap5n50 

sigmap3n50= importdata('desc_con_n50_p3.tab') ;
sigmap3n50= sigmap3n50(horizon1:horizon2,3);
sigmap3n50= sum(sigmap3n50+0.223548);
sigmap3n50

sigmap1n50= importdata('desc_con_n50_p1.tab') ;
sigmap1n50 = sigmap1n50(horizon1:horizon2,3);
sigmap1n50 = sum(sigmap1n50+0.223548);
sigmap1n50 

sigmap01n50 = importdata('desc_con_n50_p01.tab') ;
sigmap01n50 = sigmap01n50(horizon1:horizon2,3);
sigmap01n50 = sum(sigmap01n50+0.223548);
sigmap01n50 

gersh100 = importdata('gersh_100.tab') ;
gersh100 = gersh100(horizon1:horizon2,3) ;
gersh100 =sum(gersh100+0.223548);
gersh100 

gersh75 = importdata('gersh_75.tab') ;
gersh75 = gersh75(horizon1:horizon2,3) ;
gersh75=sum(gersh75+0.223548);
gersh75  

gersh50 = importdata('gersh_50.tab');
gersh50 = gersh50(horizon1:horizon2,3) ;
gersh50=sum(gersh50+0.223548);
gersh50  

gersh25 = importdata('gersh_25.tab') ;
gersh25 = gersh25(horizon1:horizon2,3) ;
gersh25=sum(gersh25+0.223548);
gersh25  
%}

%uncertian rusn
economic_mismatch=importdata('economic_mismatch.tab') ;
economic_mismatch = economic_mismatch(horizon1:horizon_unc,3) ; 
economic_mismatch = sum(economic_mismatch+0.223548) ;
economic_mismatch

tracking_mismatch=importdata('tracking_mismatch.tab') ;
tracking_mismatch = tracking_mismatch(horizon1:horizon_unc,3) ; 
tracking_mismatch = sum(tracking_mismatch+0.223548) ;
tracking_mismatch

sigmap99_mismatch = importdata('desc_con_p99_mismatch.tab') ;
sigmap99_mismatch = sigmap99_mismatch(horizon1:horizon_unc,3);
sigmap99_mismatch = sum(sigmap99_mismatch+0.223548);
sigmap99_mismatch 

sigmap9_mismatch = importdata('desc_con_p9_mismatch.tab') ;
sigmap9_mismatch = sigmap9_mismatch(horizon1:horizon_unc,3);
sigmap9_mismatch = sum(sigmap9_mismatch+0.223548);
sigmap9_mismatch 

sigmap7_mismatch = importdata('desc_con_p7_mismatch.tab') ;
sigmap7_mismatch = sigmap7_mismatch(horizon1:horizon_unc,3);
sigmap7_mismatch = sum(sigmap7_mismatch+0.223548);
sigmap7_mismatch 

sigmap5_mismatch = importdata('desc_con_p5_mismatch.tab') ;
sigmap5_mismatch = sigmap5_mismatch(horizon1:horizon_unc,3);
sigmap5_mismatch = sum(sigmap5_mismatch+0.223548);
sigmap5_mismatch 

sigmap3_mismatch = importdata('desc_con_p3_mismatch.tab') ;
sigmap3_mismatch = sigmap3_mismatch(horizon1:horizon_unc,3);
sigmap3_mismatch = sum(sigmap3_mismatch+0.223548);
sigmap3_mismatch 

sigmap1_mismatch = importdata('desc_con_p1_mismatch.tab') ;
sigmap1_mismatch = sigmap1_mismatch(horizon1:horizon_unc,3);
sigmap1_mismatch = sum(sigmap1_mismatch+0.223548);
sigmap1_mismatch 

sigmap01_mismatch = importdata('desc_con_p01_mismatch.tab') ;
sigmap01_mismatch = sigmap01_mismatch(horizon1:horizon_unc,3);
sigmap01_mismatch = sum(sigmap01_mismatch+0.223548);
sigmap01_mismatch 

gersh100_mismatch = importdata('gersh_100_mismatch.tab') ;
gersh100_mismatch = gersh100_mismatch(horizon1:horizon_unc,3) ;
gersh100_mismatch =sum(gersh100_mismatch+0.223548);
gersh100_mismatch 

gersh75_mismatch = importdata('gersh_75_mismatch.tab') ;
gersh75_mismatch = gersh75_mismatch(horizon1:horizon_unc,3) ;
gersh75_mismatch=sum(gersh75_mismatch+0.223548);
gersh75_mismatch  

gersh50_mismatch = importdata('gersh_50_mismatch.tab');
gersh50_mismatch = gersh50_mismatch(horizon1:horizon_unc,3) ;
gersh50_mismatch=sum(gersh50_mismatch+0.223548);
gersh50_mismatch  

gersh25_mismatch = importdata('gersh_25_mismatch.tab') ;
gersh25_mismatch = gersh25_mismatch(horizon1:horizon_unc,3) ;
gersh25_mismatch=sum(gersh25_mismatch+0.223548);
gersh25_mismatch  

%end uncertain


costs = [tracking,gersh100,gersh75,gersh50,gersh25,sigmap99,sigmap9,sigmap7,sigmap5,sigmap3,sigmap1,sigmap01,sigmap99n50,sigmap9n50,sigmap7n50,sigmap5n50,sigmap3n50,sigmap1n50,sigmap01n50,sigmap99n1,sigmap9n1,sigmap7n1,sigmap5n1,sigmap3n1,sigmap1n1,sigmap01n1,economic]
gershcost = [economic,gersh25,gersh50,gersh75,gersh100];
gershper=[0,25,50,75,100] ;
n2costs=[sigmap01,sigmap1,sigmap3,sigmap5,sigmap7,sigmap9,sigmap99];
sigma1n=[.01,.1,.3,.5,.7,.9,.99 ]; 
n1costs=[sigmap01n1,sigmap1n1,sigmap3n1,sigmap5n1,sigmap7n1,sigmap9n1,sigmap99n1];
sigma2n=[.01,.1,.3,.5,.7,.9,.99] ; 

mismatch_costs=[tracking_mismatch,sigmap99_mismatch,sigmap9_mismatch,sigmap7_mismatch,sigmap5_mismatch,sigmap3_mismatch,sigmap1_mismatch,sigmap01_mismatch,gersh100_mismatch,gersh75_mismatch,gersh50_mismatch,gersh25_mismatch,economic_mismatch] ;


figure
scatter(1:length(costs),costs,'.')
xlabel('Case') ;
ylabel('$\sum l(x_k,u_k) - l_{ss}$');
title('Cost comparison')
%axis([0 29 -3.5 0]) ;

figure
scatter(1:length(mismatch_costs),mismatch_costs,'.')
xlabel('Case') ;
ylabel('$\sum l(x_k,u_k) - l_{ss}$');
title('Cost comparison')
%axis([0 29 -3.5 0]) ;
%{
figure
seqcosts = [sigmap99n50(1),sigmap9s(1),sigmap7s(1),sigmap5s(1),sigmap3s(1),sigmap1s(1),sigmap01s(1)];
scatter(1:length(seqcosts),seqcosts,'.')
xlabel('Case') ;
ylabel('cost');
title('sig compare')
%}
%{
figure
scatter(gershper,gershcost,'.')
xlabel('$\%$ of Gershgorin weight') ;
ylabel('ave cost');
title('Regularized Objective')
axis([0 100 -3.5 -2]) ;

figure
scatter(sigma1n,n1costs,'.')
xlabel('$\sigma$') ;
ylabel('ave cost');
title('Stabilizing constraint, 1 norm')
axis([0 1 -3.5 -2]) ;

figure
scatter(sigma2n,n2costs,'.')
xlabel('$\sigma$') ;
ylabel('ave cost');
title('Stabilizing constraint, 2 norm')
axis([0 1 -3.5 -2]) ;
%}

%{
economic=importdata('economic_sin.tab') ;
economic = economic(horizon1:horizon2,3) ; 
economic = sum(economic) ;
economic

tracking=importdata('tracking_sin.tab') ;
tracking = tracking(horizon1:horizon2,3) ; 
tracking = sum(tracking) ;
tracking

sigmap9 = importdata('desc_con_p9_sin.tab') ;
sigmap9 = sigmap9(horizon1:horizon2,3);
sigmap9 = sum(sigmap9);
sigmap9 

sigmap7 = importdata('desc_con_p7_sin.tab') ;
sigmap7 = sigmap7(horizon1:horizon2,3);
sigmap7 = sum(sigmap7);
sigmap7 

sigmap5 = importdata('desc_con_p5_sin.tab') ;
sigmap5 = sigmap5(horizon1:horizon2,3);
sigmap5 = sum(sigmap5);
sigmap5 

sigmap3 = importdata('desc_con_p3_sin.tab') ;
sigmap3 = sigmap3(horizon1:horizon2,3);
sigmap3 = sum(sigmap3);
sigmap3 

sigmap1 = importdata('desc_con_p1_sin.tab') ;
sigmap1 = sigmap1(horizon1:horizon2,3);
sigmap1 = sum(sigmap1);
sigmap1 

sigmap01 = importdata('desc_con_p01_sin.tab') ;
sigmap01 = sigmap01(horizon1:horizon2,3);
sigmap01 = sum(sigmap01);
sigmap01 

gersh100 = importdata('gersh_100_sin.tab') ;
gersh100 = gersh100(horizon1:horizon2,3) ;
gersh100=sum(gersh100);
gersh100 

gersh75 = importdata('gersh_75_sin.tab') ;
gersh75 = gersh75(horizon1:horizon2,3) ;
gersh75=sum(gersh75);
gersh75  

gersh50 = importdata('gersh_50_sin.tab');
gersh50 = gersh50(horizon1:horizon2,3) ;
gersh50=sum(gersh50);
gersh50  

gersh25 = importdata('gersh_25_sin.tab') ;
gersh25 = gersh25(horizon1:horizon2,3) ;
gersh25=sum(gersh25);
gersh25  

costs = [lss,tracking,gersh100,gersh75,gersh50,gersh25,sigmap9,sigmap7,sigmap5,sigmap3,sigmap1,sigmap01,economic]
figure
scatter(1:length(costs),costs,'.')
xlabel('case') ;
ylabel('ave cost');
title('sin feed');
%}
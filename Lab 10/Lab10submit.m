%% Digital Signal Processing |[Lab-10]|                                       
% * Authors: _Kshitij Srivastava(1510110200)_ and _Nilambar Saha(partner)(1510110246)_
% * Lab Instructor: _Dr. Ravi Kant Saini_
% * Batch: ECE32, Friday 1 to 3 PM
% * Email: ks435@snu.edu.in
%% Objective: 
% Design of IIR filter  (In this lab, we found the 
% freqency response of various kinds of IIR filters such as LFP and HPF o 
% and Bandpass filter and Bandstp filters)
%% Program: 

% * |*Matlab Commands for First-order Low pass IIR digital filter *|

wc=pi/4;
alpha=(1-sin(wc) )/(cos(wc));

diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_lp(i)=((1-alpha)/2)* ( 1+exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_lp_abs=abs(H_lp);
H_lp_phase=phase(H_lp);
H_lp_max=max(H_lp_abs);
H_lp_3db=H_lp_max/sqrt(2);

%%
% * |*Matlab Commands for First-order Low pass filter for alpha=0.5 *|

alpha=0.5;

diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_lp_point5(i)=((1-alpha)/2)* ( 1+exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_lp_abs_point5=abs(H_lp_point5);
H_lp_phase_point5=phase(H_lp_point5);
H_lp__point5_max=max(H_lp_abs_point5);
level=H_lp__point5_max/sqrt(2);

count=1;
level1=level+0.05;
level2=level-0.05;
for i=1:61
    if level1>H_lp_abs_point5(i) && level2<H_lp_abs_point5(i)
        H_lp_point5_3db(count)=H_lp_abs_point5(i);
        H_lp_point5_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for First-order Low pass filter for alpha=0.7 *|

alpha=0.7;
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_lp_point7(i)=((1-alpha)/2)* ( 1+exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_lp_abs_point7=abs(H_lp_point7);
H_lp_phase_point7=phase(H_lp_point7);
H_lp_point7_max=max(H_lp_abs_point7);
level=H_lp_point7_max/sqrt(2);

count=1;
level1=level+0.05;
level2=level-0.05;

for i=1:61
    if level1>H_lp_abs_point7(i) && level2<H_lp_abs_point7(i)
        H_lp_point7_3db(count)=H_lp_abs_point7(i);
        H_lp_point7_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for First-order Low pass filter for alpha=0.8 *|
alpha=0.8;
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_lp_point8(i)=((1-alpha)/2)* ( 1+exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_lp_abs_point8=abs(H_lp_point8);
H_lp_phase_point8=phase(H_lp_point8);
H_lp_point8_max=max(H_lp_abs_point8);
level=H_lp_point8_max/sqrt(2);

level1=level+0.05;
level2=level-0.05;

for i=1:61
    if level1>H_lp_abs_point8(i) && level2<H_lp_abs_point8(i)
        H_lp_point8_3db(count)=H_lp_abs_point8(i);
        H_lp_point8_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for First-order High pass IIR digital filter *|
wc=pi/4;
alpha=(1-sin(wc) )/(cos(wc));

diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_hp(i)=((1+alpha)/2)* ( 1- exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_hp_abs=abs(H_hp);
H_hp_phase=phase(H_hp);
H_hp_max=max(H_hp_abs);
H_hp_3db=H_hp_max/sqrt(2);

%%
% * |*Matlab Commands for First-order High pass filter for alpha=0.5 *|

alpha=0.5;

diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_hp_point5(i)=((1+alpha)/2)* ( 1- exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_hp_abs_point5=abs(H_hp_point5);
H_hp_phase_point5=phase(H_hp_point5);
H_hp__point5_max=max(H_hp_abs_point5);
level=H_hp__point5_max/sqrt(2);

count=1;
level1=level+0.05;
level2=level-0.05;
for i=1:61
    if level1>H_hp_abs_point5(i) && level2<H_hp_abs_point5(i)
        H_hp_point5_3db(count)=H_hp_abs_point5(i);
        H_hp_point5_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for First-order High pass filter for alpha=0.7 *|

alpha=0.7;
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_hp_point7(i)=((1+alpha)/2)* ( 1- exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_hp_abs_point7=abs(H_hp_point7);
H_hp_phase_point7=phase(H_hp_point7);
H_hp_point7_max=max(H_hp_abs_point7);
level=H_hp_point7_max/sqrt(2);

count=1;
level1=level+0.05;
level2=level-0.05;

for i=1:61
    if level1>H_hp_abs_point7(i) && level2<H_hp_abs_point7(i)
        H_hp_point7_3db(count)=H_hp_abs_point7(i);
        H_hp_point7_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for First-order High pass filter for alpha=0.8 *|

alpha=0.8;
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    H_hp_point8(i)=((1+alpha)/2)* ( 1- exp(-1*1j*w) )/(1-(alpha*exp(-1*1j*w)) );
end
H_hp_abs_point8=abs(H_hp_point8);
H_hp_phase_point8=phase(H_hp_point8);
H_hp_point8_max=max(H_hp_abs_point8);
level=H_hp_point8_max/sqrt(2);

level1=level+0.05;
level2=level-0.05;

for i=1:61
    if level1>H_hp_abs_point8(i) && level2<H_hp_abs_point8(i)
        H_hp_point8_3db(count)=H_hp_abs_point8(i);
        H_hp_point8_index(count)=i;
        count=count+1;
    end
end

%%
% * |*Matlab Commands for Second order Band Pass filter *|

wo=0.4*pi;
bandw_3db=0.1*pi;
beta=cos(wo);
% cos(bandw_3db)=0.9511=(2*alpha)/(1+(alpha^2))
% alpha=1.37 or alpha=0.72
alpha=0.72;

h_bandpass=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass(i)=h_bandpass(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs=abs(h_bandpass);
h_bandpass_phase=phase(h_bandpass);
level=max(h_bandpass_abs)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;
for i=1:61
    if level1>h_bandpass_abs(i) && level2<h_bandpass_abs(i)
        h_bandpass_3db(count)=h_bandpass_abs(i);
        H_bandpass_index(count)=i;
        count=count+1;
    end
end
h_bandpass_3db_bandwidth=41-20;
Q=wo/(h_bandpass_3db_bandwidth*diff);

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.5 and alpha=0.2 *|

beta=0.5;
alpha=0.2;
h_bandpass_a2b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass_a2b5(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) ));
    h_bandpass_a2b5(i)=h_bandpass_a2b5(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs_a2b5=abs(h_bandpass_a2b5);
h_bandpass_phase_a2b5=phase(h_bandpass_a2b5);
level=max(h_bandpass_abs_a2b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;
for i=1:61
    if level1>h_bandpass_abs_a2b5(i) && level2<h_bandpass_abs_a2b5(i)
        h_bandpass_a2b5_3db(count)=h_bandpass_abs_a2b5(i);
        H_bandpass_a2b5_index(count)=i;
        count=count+1;
    end
end
h_bandpass_a2b5_3db_band=(42-19)*diff;
Q_bandpass_a2b5=acos(beta)/h_bandpass_a2b5_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.5 and alpha=0.5 *|

beta=0.5;
alpha=0.5;

h_bandpass_a5b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass_a5b5(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass_a5b5(i)=h_bandpass_a5b5(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs_a5b5=abs(h_bandpass_a5b5);
h_bandpass_phase_a5b5=phase(h_bandpass_a5b5);
level=max(h_bandpass_abs_a5b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandpass_abs_a5b5(i) && level2<h_bandpass_abs_a5b5(i)
        h_bandpass_a5b5_3db(count)=h_bandpass_abs_a5b5(i);
        H_bandpass_a5b5_index(count)=i;
        count=count+1;
    end
end
h_bandpass_a5b5_3db_band=(42-19)*diff;
Q_bandpass_a5b5=acos(beta)/h_bandpass_a5b5_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.5 and alpha=0.8 *|

beta=0.5;
alpha=0.8;

h_bandpass_a8b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass_a8b5(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass_a8b5(i)=h_bandpass_a8b5(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs_a8b5=abs(h_bandpass_a8b5);
h_bandpass_phase_a8b5=phase(h_bandpass_a8b5);
level=max(h_bandpass_abs_a8b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandpass_abs_a8b5(i) && level2<h_bandpass_abs_a8b5(i)
        h_bandpass_a8b5_3db(count)=h_bandpass_abs_a8b5(i);
        h_bandpass_a8b5_index(count)=i;
        count=count+1;
    end
end
h_bandpass_a8b5_3db_band=(42-19)*diff;
Q_bandpass_a8b5=acos(beta)/h_bandpass_a8b5_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.1 and alpha=0.6 *|

beta=0.1;
alpha=0.6;

h_bandpass_b1a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass_b1a6(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass_b1a6(i)=h_bandpass_b1a6(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2); 
end

h_bandpass_abs_b1a6=abs(h_bandpass_b1a6);
h_bandpass_phase_b1a6=phase(h_bandpass_b1a6);
level=max(h_bandpass_abs_b1a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandpass_abs_b1a6(i) && level2<h_bandpass_abs_b1a6(i)
        h_bandpass_b1a6_3db(count)=h_bandpass_abs_b1a6(i);
        h_bandpass_b1a6_index(count)=i;
        count=count+1;
    end
end
h_bandpass_b1a6_3db_band=(47-14)*diff;
Q_bandpass_b1a6=acos(beta)/h_bandpass_b1a6_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.5 and alpha=0.6 *|

beta=0.5;
alpha=0.6;

h_bandpass_b5a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff); 
    h_bandpass_b5a6(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass_b5a6(i)=h_bandpass_b5a6(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs_b5a6=abs(h_bandpass_b5a6);
h_bandpass_phase_b5a6=phase(h_bandpass_b5a6);
level=max(h_bandpass_abs_b5a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandpass_abs_b5a6(i) && level2<h_bandpass_abs_b5a6(i)
        h_bandpass_b5a6_3db(count)=h_bandpass_abs_b5a6(i);
        h_bandpass_b5a6_index(count)=i;
        count=count+1;
    end
end
h_bandpass_b5a6_3db_band=(43-18)*diff;
Q_bandpass_b5a6=acos(beta)/h_bandpass_b5a6_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Pass filter for beta=0.8 and alpha=0.6 *|

beta=0.8;
alpha=0.5;

h_bandpass_b8a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandpass_b8a6(i)=(1/(1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w))));
    h_bandpass_b8a6(i)=h_bandpass_b8a6(i)*( 1-exp(-2*1j*w) )*((1-alpha)/2);
end

h_bandpass_abs_b8a6=abs(h_bandpass_b8a6);
h_bandpass_phase_b8a6=phase(h_bandpass_b8a6);
level=max(h_bandpass_abs_b8a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandpass_abs_b8a6(i) && level2<h_bandpass_abs_b8a6(i)
        h_bandpass_b8a6_3db(count)=h_bandpass_abs_b8a6(i);
        h_bandpass_b8a6_index(count)=i;
        count=count+1;
    end
end
h_bandpass_b8a6_3db_band=(40-20)*diff;
Q_bandpass_b8a6=acos(beta)/h_bandpass_b8a6_3db_band;

%%
% * |*Matlab Commands for Second-order Band stop IIR digital filter *|

wo=0.4*pi;
threedB_bandw=pi*0.1;

beta=acos(wo);
% cos(threedB_bandw)=0.9511=2*alpha/(1+alpha^2)
% alpha=0.72,1.37 so alpha=0.72
alpha=0.72;
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop(i)=h_bandstop(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop(i)= h_bandstop(i)*((1+alpha)/2);
end

h_bandstop_abs=abs(h_bandstop);
h_bandstop_phase=phase(h_bandstop);
level=max(h_bandstop_abs)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs(i) && level2<h_bandstop_abs(i)
        h_bandstop_3db(count)=h_bandstop_abs(i);
        h_bandstop_index(count)=i;
        count=count+1;
    end
end

h_bandstop_3db_band=(23-8)*diff;
Q_bandstop=acos(beta)/h_bandstop_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Stop filter for beta=0.5 and alpha=0.2 *|

beta=0.5;
alpha=0.2;
h_bandstop_a2b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop_a2b5(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_a2b5(i)=h_bandstop_a2b5(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_a2b5(i)= h_bandstop_a2b5(i)*((1+alpha)/2);
end
h_bandstop_abs_a2b5=abs(h_bandstop_a2b5);
h_bandstop_phase_a2b5=phase(h_bandstop_a2b5);
level=max(h_bandstop_abs_a2b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;
for i=1:61
    if level1>h_bandstop_abs_a2b5(i) && level2<h_bandstop_abs_a2b5(i)
        h_bandstop_a2b5_3db(count)=h_bandstop_abs_a2b5(i);
        H_bandstop_a2b5_index(count)=i;
        count=count+1;
    end
end
h_bandstop_a2b5_3db_band=(36-25)*diff;
Q_bandstop_a2b5=acos(beta)/h_bandstop_a2b5_3db_band;
%%

% * |*Matlab Commands for 2nd order Band Stop filter for beta=0.5 and alpha=0.5 *|

beta=0.5;
alpha=0.5;

h_bandstop_a5b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop_a5b5(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_a5b5(i)=h_bandstop_a5b5(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_a5b5(i)= h_bandstop_a5b5(i)*((1+alpha)/2);
end

h_bandstop_abs_a5b5=abs(h_bandstop_a5b5);
h_bandstop_phase_a5b5=phase(h_bandstop_a5b5);
level=max(h_bandstop_abs_a5b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs_a5b5(i) && level2<h_bandstop_abs_a5b5(i)
        h_bandstop_a5b5_3db(count)=h_bandstop_abs_a5b5(i);
        H_bandstop_a5b5_index(count)=i;
        count=count+1;
    end
end

h_bandstop_a5b5_3db_band=(44-17)*diff;
Q_bandstop_a5b5=acos(beta)/h_bandstop_a5b5_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Stop filter for beta=0.5 and alpha=0.8 *|

beta=0.5;
alpha=0.8;

h_bandstop_a8b5=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop_a8b5(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_a8b5(i)=h_bandstop_a8b5(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_a8b5(i)= h_bandstop_a8b5(i)*((1+alpha)/2);
end

h_bandstop_abs_a8b5=abs(h_bandstop_a8b5);
h_bandstop_phase_a8b5=phase(h_bandstop_a8b5);
level=max(h_bandstop_abs_a8b5)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs_a8b5(i) && level2<h_bandstop_abs_a8b5(i)
        h_bandstop_a8b5_3db(count)=h_bandstop_abs_a8b5(i);
        h_bandstop_a8b5_index(count)=i;
        count=count+1;
    end
end
h_bandstop_a8b5_3db_band=(42-19)*diff;
Q_bandstop_a8b5=acos(beta)/h_bandstop_a8b5_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Stop filter for beta=0.1 and alpha=0.6 *|

beta=0.1;
alpha=0.6;

h_bandstop_b1a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop_b1a6(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_b1a6(i)=h_bandstop_b1a6(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_b1a6(i)= h_bandstop_b1a6(i)*((1+alpha)/2); 
end

h_bandstop_abs_b1a6=abs(h_bandstop_b1a6);
h_bandstop_phase_b1a6=phase(h_bandstop_b1a6);
level=max(h_bandstop_abs_b1a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs_b1a6(i) && level2<h_bandstop_abs_b1a6(i)
        h_bandstop_b1a6_3db(count)=h_bandstop_abs_b1a6(i);
        h_bandstop_b1a6_index(count)=i;
        count=count+1;
    end
end
h_bandstop_b1a6_3db_band=(47-14)*diff;
Q_bandstop_b1a6=acos(beta)/h_bandstop_b1a6_3db_band;

%%
% * |*Matlab Commands for 2nd order Band stop filter for beta=0.5 and alpha=0.6 *|

beta=0.5;
alpha=0.6;

h_bandstop_b5a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff); 
    h_bandstop_b5a6(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_b5a6(i)=h_bandstop_b5a6(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_b5a6(i)=h_bandstop_b5a6(i)*((1+alpha)/2);
end

h_bandstop_abs_b5a6=abs(h_bandstop_b5a6);
h_bandstop_phase_b5a6=phase(h_bandstop_b5a6);
level=max(h_bandstop_abs_b5a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs_b5a6(i) && level2<h_bandstop_abs_b5a6(i)
    h_bandstop_b5a6(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_b5a6(i)=h_bandstop_b5a6(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_b5a6(i)=h_bandstop_b5a6(i)*((1+alpha)/2);
    count=count+1;
    end
end
h_bandstop_b5a6_3db_band=(43-18)*diff;
Q_bandstop_b5a6=acos(beta)/h_bandstop_b5a6_3db_band;

%%
% * |*Matlab Commands for 2nd order Band Stop filter for beta=0.8 and alpha=0.6 *|

beta=0.8;
alpha=0.5;

h_bandstop_b8a6=zeros(1,61);
diff=(2*pi)/61;
for i=1:61
    w=(-1*pi)+ (i*diff);
    h_bandstop_b8a6(i)=1/( 1-(beta*(alpha+1)*exp(-1*1j*w))+(alpha*exp(-2*1j*w)) );
    h_bandstop_b8a6(i)=h_bandstop_b8a6(i)*(1-(2*beta*exp(-1*1j*w))+exp(-2*1j*w));
    h_bandstop_b8a6(i)=h_bandstop_b8a6(i)*((1+alpha)/2);
end

h_bandstop_abs_b8a6=abs(h_bandstop_b8a6);
h_bandstop_phase_b8a6=phase(h_bandstop_b8a6);
level=max(h_bandstop_abs_b8a6)/sqrt(2);

level1=level+0.05;
level2=level-0.05;
count=1;

for i=1:61
    if level1>h_bandstop_abs_b8a6(i) && level2<h_bandstop_abs_b8a6(i)
        h_bandstop_b8a6_3db(count)=h_bandstop_abs_b8a6(i);
        h_bandstop_b8a6_index(count)=i;
        count=count+1;
    end
end
h_bandstop_b8a6_3db_band=(40-20)*diff;
Q_bandstop_b8a6=acos(beta)/h_bandstop_b8a6_3db_band;

%% Result
% * |*Plot for the Question No 1(a)*|
figure;plot(H_lp_abs);
title('First order LPF Magnitude Response');
xlabel('Frequencies');ylabel('Amplitude');

figure;plot(H_lp_phase);
title('First order LPF Phase Response');
xlabel('Frequencies');ylabel('Angle');

%%
% * |*Plot for the Question No 1(b)*|
figure;plot(H_lp_abs_point5);
title('First order LPF Magnitude Response plot','r');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(H_lp_abs_point7,'g');
hold on
plot(H_lp_abs_point8,'b');
legend('\color{red}alpha=0.7 ','\color{green}alpha=0.5','\color{blue}alpha=0.8');

%%
% * |*Plot for the Question No 1(b)*|
figure;plot(H_lp_phase_point5,'r');
title('First order LPF Phase Response plot');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(H_lp_phase_point7,'g');
hold on
plot(H_lp_phase_point8,'b');
legend('\color{red}alpha=0.5 ','\color{green}alpha=0.7','\color{blue}alpha=0.8');

%%
% * |*Plot for the Question No 2(a)*|
figure;plot(H_hp_abs);
title('First order HPF Magnitude Response');
xlabel('Frequencies');ylabel('Amplitude');

figure;plot(H_hp_phase);
title('First order HPF Phase Response');
xlabel('Frequencies');ylabel('Angle');

%%
% * |*Plot for the Question No 2(b)*|

figure;plot(H_hp_abs_point5);
title('First order LPF Magnitude Response plot','r');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(H_hp_abs_point7,'g');
hold on
plot(H_hp_abs_point8,'b');
legend('\color{red}alpha=0.5 ','\color{green}alpha=0.7','\color{blue}alpha=0.8');


%%
% * |*Plot for the Question No 2(b)*|
figure;plot(H_hp_phase_point5);
title('First order LPF Phase Response plot','r');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(H_hp_phase_point7,'g');
hold on
plot(H_hp_phase_point8,'b');
legend('\color{red}alpha=0.7 ','\color{green}alpha=0.7','\color{blue}alpha=0.8');

%%
% * |*Plot for the Question No 3(a)*|
figure;plot(h_bandpass_abs);
title('Second-order Band pass filter Magnitude Response');
xlabel('Frequencies');ylabel('Amplitude');

figure;plot(h_bandpass_phase);
title('Second-order Band pass filter Phase Response');
xlabel('Frequencies');ylabel('Angle');

%%
% * |*Plot for the Question No 3(b)*|
figure;plot(h_bandpass_abs_a2b5,'r');
title(' Second-order Band pass filter magnitude Response plot');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(h_bandpass_abs_a5b5,'g');
hold on
plot(h_bandpass_abs_a8b5,'b');
legend('\color{red}alpha=0.2','\color{green}alpha=0.5','\color{blue}alph=0.8');

%%
% * |*Plot for the Question No 3(b)*|
figure;plot(h_bandpass_phase_a2b5,'r');
title(' Second-order Band pass filter Phase Response plot');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(h_bandpass_phase_a5b5,'g');
hold on
plot(h_bandpass_phase_a8b5,'b');
legend('\color{red}alpha=0.2','\color{green}alpha=0.5','\color{blue}alpha=0.8');

%%
% * |*Result for the Question No 3(b)*|
% 3db Bandwidth for a=0.2 and beta=0.5 
h_bandpass_a2b5_3db_band
% Quality factor for a=0.2 and beta=0.5 
Q_bandpass_a2b5
% 3db Bandwidth for a=0.5 and beta=0.5 
h_bandpass_a5b5_3db_band
% Quality factor for a=0.5 and beta=0.5 
Q_bandpass_a5b5
% 3db Bandwidth for a=0.8 and beta=0.5 
h_bandpass_a8b5_3db_band
% Quality factor for a=0.8 and beta=0.5 
Q_bandpass_a8b5

%%
% * |*Plot for the Question No 3(d)*|
figure;plot(h_bandpass_abs_b1a6,'r');
title(' Second-order Band pass filter Amplitude Response plot');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(h_bandpass_abs_b5a6,'g');
hold on
plot(h_bandpass_abs_b8a6,'b');
legend('\color{red}beta=0.1','\color{green}beta=0.5','\color{blue}beta=0.8');

%%
% * |*Plot for the Question No 3(d)*|
figure;plot(h_bandpass_phase_b1a6,'r');
title(' Second-order Band pass filter Phase Response plot');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(h_bandpass_phase_b5a6,'g');
hold on
plot(h_bandpass_phase_b8a6,'b');
legend('\color{red}beta=0.1','\color{green}beta=0.5','\color{blue}beta=0.8');

%%
% * |*Result for the Question No 3(d)*|
% 3db Bandwidth for a=0.6 and beta=0.1 
h_bandpass_b1a6_3db_band
% Quality factor for a=0.6 and beta=0.1 
Q_bandpass_b1a6
% 3db Bandwidth for a=0.6 and beta=0.5
h_bandpass_b5a6_3db_band
% Quality factor for a=0.6 and beta=0.5 
Q_bandpass_b5a6
% 3db Bandwidth for a=0.6 and beta=0.8
h_bandpass_b8a6_3db_band
% Quality factor for a=0.6 and beta=0.8 
Q_bandpass_b8a6

%%
% * |*Plot for the Question No 4(a)*|
figure;plot(h_bandstop_abs);
title('Second-order Band stop filter Magnitude Response');
xlabel('Frequencies');ylabel('Amplitude');

figure;plot(h_bandstop_phase);
title('Second-order Band stop filter Phase Response');
xlabel('Frequencies');ylabel('Angle');

%%
% * |*Plot for the Question No 4(b)*|
figure;plot(h_bandstop_abs_a2b5,'r');
title(' Second-order Band stop filter magnitude Response plot');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(h_bandstop_abs_a5b5,'g');
hold on
plot(h_bandstop_abs_a8b5,'b');
legend('\color{red}alpha=0.2','\color{green}alpha=0.5','\color{blue}alpha=0.8');

%%
% * |*Plot for the Question No 4(b)*|
figure;plot(h_bandstop_phase_a2b5,'r');
title(' Second-order Band stop filter Phase Response plot');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(h_bandstop_phase_a5b5,'g');
hold on
plot(h_bandstop_phase_a8b5,'b');
legend('\color{red}alpha=0.2','\color{green}alpha=0.5','\color{blue}alpha=0.8');

%%
% * |*Result for the Question No 4(b)*|
% 3db Bandwidth for a=0.2 and beta=0.5 
h_bandstop_a2b5_3db_band
% Quality factor for a=0.2 and beta=0.5 
Q_bandstop_a2b5
% 3db Bandwidth for a=0.5 and beta=0.5 
h_bandstop_a5b5_3db_band
% Quality factor for a=0.5 and beta=0.5 
Q_bandstop_a5b5
% 3db Bandwidth for a=0.8 and beta=0.5 
h_bandstop_a8b5_3db_band
% Quality factor for a=0.8 and beta=0.5 
Q_bandstop_a8b5

%%
% * |*Plot for the Question No 3(d)*|
figure;plot(h_bandstop_abs_b1a6,'r');
title(' Second-order Band stop filter Amplitude Response plot');
xlabel('Frequencies');ylabel('Amplitude');
hold on
plot(h_bandstop_abs_b5a6,'g');
hold on
plot(h_bandstop_abs_b8a6,'y');
legend('\color{red}beta=0.1','\color{green}beta=0.5','\color{blue}beta=0.8');

%%
% * |*Plot for the Question No 3(d)*|
figure;plot(h_bandstop_phase_b1a6,'r');
title(' Second-order Band stop filter Phase Response plot');
xlabel('Frequencies');ylabel('Phase');
hold on
plot(h_bandstop_phase_b5a6,'g');
hold on
plot(h_bandstop_phase_b8a6,'b');
legend('\color{red}beta=0.1','\color{green}beta=0.5','\color{blue}beta=0.8');

%%
% * |*Result for the Question No 3(d)*|
% 3db Bandwidth for a=0.6 and beta=0.1 
h_bandstop_b1a6_3db_band
% Quality factor for a=0.6 and beta=0.1 
Q_bandstop_b1a6
% 3db Bandwidth for a=0.6 and beta=0.5
h_bandstop_b5a6_3db_band
% Quality factor for a=0.6 and beta=0.5 
Q_bandstop_b5a6
% 3db Bandwidth for a=0.6 and beta=0.8
h_bandstop_b8a6_3db_band
% Quality factor for a=0.6 and beta=0.8 
Q_bandstop_b8a6

%%
%initial
clc;clear;
close all;
warning off;
data=xlsread('ATR of skin (1).xlsx');

figure;
plot(data(:,1),data(:,2))
title('ATR-IR spectrum');
xlabel('wavenumber(cm-1)');
ylabel('absorbance(OD)');
grid;
hale=xlsread('Hale.xlsx');
%parameters
absorbance=data(:,2)./10;%2.5 bounces
R_TE=(10.^(-absorbance));
lamdba=1e7./data(:,1); %nm
figure;
plot(lamdba,R_TE)
title('Square modulus of the complex amplitude reflection coefficient')
xlabel('wavelength(nm)');
ylabel('R_T_E');
grid;
figure;
plot(lamdba,data(:,2))
title('ATR-IR spectrum');
xlabel('wavenlength(nm)');
ylabel('absorbance(OD)');
grid;
theta=pi/4;   %angle of incidence
n_1=2.38;     %refraction index for prism
n_2=1.32;     %refraction index for sample %1.3 when in
%reflection coefficient
n_r=n_2/n_1;
step=1
%newton

t0=tic;
for i=1:length(R_TE)
%reflection coefficient
syms n_i;
func=@(n_i) R_TE(i)-abs((cos(theta)-sqrt((n_r^2-n_i^2-sin(theta)^2)+2*1i*n_r*n_i))/(cos(theta)+sqrt((n_r^2-n_i^2-sin(theta)^2)+2*1i*n_r*n_i)))^2;
%f=feval(fun,0);
dfun=diff(func(n_i));
%f22=double(subs(dfun,n_i,0));
it_max=150;
ep=1e-20;
x=0;k=0;
while k<it_max
x1=x;f1=feval(func,x);
f2=real(double(subs(dfun,n_i,x)));
x=x-f1/f2;
%误差区间
err = abs(x-x1);
%相对误差区间，可以避免发散
relerr = 2*err/(abs(x1)+ep);
y=feval(func,x);
epsilon=1e-10;
    if (err<ep)||(relerr<ep)||(abs(y)<epsilon)
        index=1;break;
    end
    k=k+1;
end
nn(i)=double(x*n_1);
i
end
toc(t0);
%{
%math
step=2
syms n_img
re_eq=2*(n_img^4 + (22873*n_img^2)/14161 + 29691601/802135684)^(1/2) - 1;
img_eq=(n_img^2 + (n_img^4 + (22873*n_img^2)/14161 + 29691601/802135684)^(1/2) + 5449/28322)^(1/2)*2;
mu_eq=(2*((n_img^4 + (22873*n_img^2)/14161 + 29691601/802135684)^(1/2) - n_img^2 - 5449/28322)^(1/2) + 2*(n_img^4 + (22873*n_img^2)/14161 + 29691601/802135684)^(1/2) + 1);
equ_2=(re_eq/mu_eq)^2+(img_eq/mu_eq)^2;
step=3
t1=tic;
for j=1:length(R_TE)
sola=solve(equ_2==R_TE(j),n_img);  %待求解的变量是n_img
n_i(j)=double(abs(sola(1))*n_1);
run=j
end
toc(t1);
figure;
plot(lamdba,n_i);hold on
plot(lamdba,n_i,'--');
plot(hale(1:124,1)*1000,hale(1:124,2));hold off
legend('Netwen-Raphson method with 10 bounchs','Symbolic maths method with 10 bounchs','Hale & Querry data');
title('n_i');
xlabel('wavelength(nm)');
ylabel('n_i');
grid;
%}
step=4

figure;
plot(lamdba,nn);
title('n_i');
xlabel('wavelength(nm)');
ylabel('n_i');
grid;
%scattering coefficient
lamdba_1=970; %nm
mu_s_1=7.1;%mm-1
lamdba_2=1600;%nm
mu_s_2=6.3;%mm-1
syms  aa nnn
eq1=aa/(lamdba_1^nnn)==mu_s_1;
eq2=aa/(lamdba_2^nnn)==mu_s_2;
A=solve(eq1,eq2,aa,nnn);

a=double(A.aa);
n=double(A.nnn);
lamdba=1e7./data(:,1); %nm
mu_s=a./lamdba.^n;%mm-1

%Beer-lambert law
%absorption coefficient 
for j=1:127
mu_a_H(j)=(hale(j,2)/n_1)*4*pi/(hale(j,1)*1e3*1e-6);
end
for j=1:length(R_TE)
mu_a(j)=((nn(j)/n_1)*4*pi/(lamdba(j)*1e-6));
end
mu_a=mu_a';

figure; 
%subplot 221;
plot(lamdba,mu_a);
xlabel('wavelength(nm)');
ylabel('mu_a(mm-1)');
title('absorption coefficient');grid;
%subplot 222;
figure;
plot(lamdba,mu_s)
xlabel('wavelength(nm)');
ylabel('mu_s(mm-1)');
title('scattering coefficient');grid;
%total attenuation
mu_t=mu_a+mu_s;
%subplot 223;
figure;
plot(lamdba,mu_t)
xlabel('wavelength(nm)');
ylabel('mu_t(mm-1)');
title('total attenuation coefficient');grid;
%backscatter coefficient
mu_b=0.05*mu_s;
%subplot 224;
figure
plot(lamdba,mu_b)
title('backscatter coefficient');grid;
xlabel('wavelength(nm)');
ylabel('mu_b(mm-1)');


step=5
syms z

for j=1:length(R_TE)
Iz=exp(-mu_a(j)*z);
I(j)=double(int(Iz,z,0,inf));
Io(j)=(mu_b(j)*I(j));
OD_ocs(j)=double(-log10(Io(j)));
run1=j
end
figure;
plot(data(:,1),OD_ocs)
title('OCS spectrum');
xlabel('wavenumber(cm-1)');
ylabel('absorbance(OD)');
grid;
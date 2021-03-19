# Fuel-cell

clear all 

Tfc=1073.15;%K 

f=5.5; 

uf=0.44; 

 

  

%statheres 

Pfc=1;%bar 

R=8.314e3;%kJ/kmol K 

iL=5;%A/cm2 

io=1;%A/cm2 

Rw=0.1;%? cm2 

T=298.15 

F=96485.3329e3;%J/kmol V 

 

%revma 3 

T3=575+273.15; 

N3=4.82840102273176;%kmol/s 

y3=[0.043658292668462, 0.302243554851983, 3.88854556685154e-002, 0.490648498787850, 9.34980225073141e-002,	0, 3.10661755158753e-002]; 

n3=N3*y3; 

 

% 1.CH4 2.H2O 3.CO 4.H2 5.CO2 6.O2 7.N2 

r=uf*n3(4);%kmol/s 

 

%revma 4 

T4=T3; 

N4=(n3(1)+n3(3)+n3(4))*f;%kmol/s 

y4=[0,0, 0, 0, 0, 0.21, 0.79]; 

n4=N4*y4; 

% 1.CH4 2.H2O 3.CO 4.H2 5.CO2 6.O2 7.N2 

 

%revma 5 

T5=Tfc; 

n5=[n3(1),n3(2)+r,n3(3),n3(4)-r,n3(5),0,n3(7)]; 

N5=sum(n5); 

y5=n5./N5; 

% 1.CH4 2.H2O 3.CO 4.H2 5.CO2 6.O2 7.N2 

 

%revma 6 

T6=Tfc; 

n6=[0,0, 0, 0, 0 ,n4(6)-0.5*r, n4(7)]; 

N6=sum(n6); 

y6=n6./N6; 

% 1.CH4 2.H2O 3.CO 4.H2 5.CO2 6.O2 7.N2 

 

Nout=N5+N6; 

nout=n5+n6; 

yout=nout/Nout; 

 

y=[y3; y4; yout]; 

Ti=[T3, T3, T5]; 

 

for j=1:3 

    

h(j)=enthalpy(y(j,:),Ti(j));%kJ/kmol 

 

end 

 

Wel=(N3*h(1)+N4*h(2)-(N5+N6)*h(3));%kJ/s 

 

 

I=2*F*r;%W V 

 

V=(Wel*1000/I);%V 

 

syms i 

g=1.2528-(2.307e-4)*Tfc+((R*Tfc)/(2*F))*log((y5(5)*sqrt(y6(7)))/y5(4))-(Rw*i+((2*R*Tfc)/(F))*asinh(i/(io*2))+(R*Tfc/(2*F))*log(iL/(iL-i)))-V; 

ii=solve(g,i); 

 

A=I/ii; 

      Aa=double(A) 

 

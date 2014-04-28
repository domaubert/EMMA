

nH=0.2;
//x=1e-6;
temp=spanl(1e1,1e9,2048);
CLUMPF=1.;
z=20.;
nh2=nH*1e-6*(1+z)^3;// ! m-3 ==> cm-3

//case A recombination rate
lambda=2e0*157807e0/temp;
alpha_a=1.269e-13*pow(lambda,1.503)/pow(1e0+pow(lambda/0.522,0.470),1.923); //cm3/s


// collisional ionisation rate
T5=temp/1e5;
beta=5.85e-11*sqrt(temp)/(1+sqrt(T5))*exp(-(157809e0/temp));

//equilibrium ionisation ratio

x=1.-1./(1.+beta/alpha_a);


// Collisional Ionization Cooling

c1=exp(-157809.1e0/temp)*1.27e-21*sqrt(temp)/(1.+sqrt(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;


// Case A Recombination Cooling

c2=1.778e-29*temp*pow(2e0*157807e0/temp,1.965e0)/pow(1.+pow(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2*CLUMPF;
  

// Collisional excitation cooling

c3=exp(-118348e0/temp)*7.5e-19/(1+sqrt(temp/1e5))*x*(1.-x)*nh2*nh2*CLUMPF;


// Bremmsstrahlung

c4=1.42e-27*1.5e0*sqrt(temp)*x*x*nh2*nh2*CLUMPF;

// Compton Cooling

c5=1.017e-37*pow(2.727*(1+z),4)*(temp-2.727*(1+z))*nh2*x;
//c5=(c5>0.?c5:0); // to be checked

// Overall Cooling

lambda=c1+c2+c3+c4+c5;// ! erg*cm-3*s-1
  

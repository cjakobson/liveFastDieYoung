%predict population fraction of BIG+/naive given growth rate and death rate
%differences 

function [t,BIGfrac]=oscillate(params)



x0=zeros(params.N,1);
x1=zeros(params.N,1);
t=zeros(params.N,1);

%set initial population ratio
x0(1)=params.x00;
x1(1)=params.x10;
t(1)=0;

for i=2:params.N %number of oscillations
    
    %first do nutrient phase
    temp1=x0(i-1)*exp(params.u0*params.t1);
    temp2=x1(i-1)*exp(params.u1*params.t1);
    
    %then death phase
    x0(i)=temp1*exp(-params.d0*params.t2);
    x1(i)=temp2*exp(-params.d1*params.t2);
    
    t(i)=t(i-1)+params.t1+params.t2;
    
end

BIGfrac=x1./x0;


function x0=getInitialGuess(n,Wi,Aphi,RHSphi,parameters)

Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions
nadd=0; % number of additional variables
if(parameters.bcType==3)
    nadd=3;
end
x0=zeros(2*n+nadd,1); % holder for initial guess


%initial guess
x0(Nw,1)=Wi;

x0(Nphi,1)= Aphi\RHSphi(x0(Nw,1),0.*x0(Nw,1)); %use 1st step from picard iteration to get initial guess for phi





end
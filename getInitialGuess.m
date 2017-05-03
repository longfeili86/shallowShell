function x0=getInitialGuess(n,W0,PHI0,Aphi,RHSphi,parameters)

Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions
nadd=0; % number of additional variables
if(parameters.bcType==3)
    nadd=3;
end
x0=zeros(2*n+nadd,1); % holder for initial guess



%initial guess is (PHI0,W0)

PHI0= Aphi\RHSphi(W0,PHI0); %use 1st step from picard iteration as PHI0 for fsolve


x0(Nphi,1)=PHI0;
x0(Nw,1)=W0;


end
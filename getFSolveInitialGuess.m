function x0=getFSolveInitialGuess(W0,PHI0,bcType)

x0=[PHI0;W0];
if(bcType==3)
    x0=[x0;0;0;0]; % add 3 additional initial guess for lambdas
end

end
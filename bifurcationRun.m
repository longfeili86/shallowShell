% bifurcation run

clear
runShell -noplot -case=coupledSystem -solver=newton -xi=-1090 -funcDefFile=bifurcationFuncDef1 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults -saveICFile=bifurIC1 -nonlinear -tol=1e-6 -maxIter=200
runShell -noplot -case=coupledSystem -solver=newton -xi=-1090 -funcDefFile=bifurcationFuncDef2 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults -saveICFile=bifurIC2 -nonlinear -tol=1e-6 -maxIter=200
runShell -noplot -case=coupledSystem -solver=newton -xi=-1090 -funcDefFile=bifurcationFuncDef3 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults -saveICFile=bifurIC3 -nonlinear -tol=1e-6 -maxIter=200

Xi=-1090:1:-1080;

for i=1:length(Xi)
    
xi=Xi(i);  


cmd=sprintf('runShell -noplot -case=coupledSystem -solver=imPicard -funcDefFile=bifurcationFuncDef1 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults_Xi%iimPicard -xi=%f -saveICFile=bifurICimPicard -readICFile=bifurIC1 -nonlinear -tol=1e-6 -maxIter=200',xi,xi);
fprintf('%s\n',cmd);
eval(cmd);


cmd=sprintf('runShell -noplot -case=coupledSystem -solver=newton -funcDefFile=bifurcationFuncDef1 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults_Xi%ia -xi=%f -saveICFile=bifurIC1 -readICFile=bifurICimPicard -nonlinear -tol=1e-6 -maxIter=200',xi,xi);
fprintf('%s\n',cmd);
eval(cmd);


cmd=sprintf('runShell -noplot -case=coupledSystem -solver=newton -funcDefFile=bifurcationFuncDef2 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults_Xi%ib -xi=%f -saveICFile=bifurIC2 -readICFile=bifurIC2 -nonlinear -tol=1e-6 -maxIter=200',xi,xi);
fprintf('%s\n',cmd);
eval(cmd);

cmd=sprintf('runShell -noplot -case=coupledSystem -solver=newton -funcDefFile=bifurcationFuncDef3 -nx=40 -ny=40 -bcType=1 -f=bifurcationTestResults_Xi%ic -xi=%f -saveICFile=bifurIC3 -readICFile=bifurIC3 -nonlinear -tol=1e-6 -maxIter=200',xi,xi);
fprintf('%s\n',cmd);
eval(cmd);


end



for i=1:length(Xi)
    xi=Xi(i); 
    results=sprintf('bifurcationTestResults_Xi%ia/results.mat',xi);
    load(results);
    wa(i)=Wplot(20,20);
  
    results=sprintf('bifurcationTestResults_Xi%ib/results.mat',xi);
    load(results);
    wb(i)=Wplot(20,20);
    
    results=sprintf('bifurcationTestResults_Xi%ic/results.mat',xi);
    load(results);
    wc(i)=Wplot(20,20);    
    
    results=sprintf('bifurcationTestResults_Xi%iimPicard/results.mat',xi);
    load(results);
    wp(i)=Wplot(20,20);  
end

figure
hold on
plot(Xi,wa,'r>')
plot(Xi,wb,'g*')
plot(Xi,wc,'bo')
plot(Xi,wp,'k+')
hold off
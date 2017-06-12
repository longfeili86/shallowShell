% this script contains the commands for convRate with Free BC. 
% since Free requres more resolution to converge, we do it
% separately here rather in the convRate.m


% exPicard
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=exPicardNonlinearcoupledSystemTestFreeG8 -nonlinear;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=exPicardNonlinearcoupledSystemTestFreeG16 -nonlinear;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=exPicardNonlinearcoupledSystemTestFreeG32 -nonlinear;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=exPicardNonlinearcoupledSystemTestFreeG64 -nonlinear;

runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=exPicardLinearcoupledSystemTestFreeG8;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=exPicardLinearcoupledSystemTestFreeG16;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=exPicardLinearcoupledSystemTestFreeG32;
runShell -case=coupledSystem -solver=exPicard -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=exPicardLinearcoupledSystemTestFreeG64;



% newton
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=newtonNonlinearcoupledSystemTestFreeG8 -nonlinear;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=newtonNonlinearcoupledSystemTestFreeG16 -nonlinear;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=newtonNonlinearcoupledSystemTestFreeG32 -nonlinear;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=newtonNonlinearcoupledSystemTestFreeG64 -nonlinear;

runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=newtonLinearcoupledSystemTestFreeG8;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=newtonLinearcoupledSystemTestFreeG16;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=newtonLinearcoupledSystemTestFreeG32;
runShell -case=coupledSystem -solver=newton -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=newtonLinearcoupledSystemTestFreeG64;



% fsolve
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=fsolveNonlinearcoupledSystemTestFreeG8 -nonlinear;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=fsolveNonlinearcoupledSystemTestFreeG16 -nonlinear;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=fsolveNonlinearcoupledSystemTestFreeG32 -nonlinear;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=fsolveNonlinearcoupledSystemTestFreeG64 -nonlinear;

runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=80 -ny=80 -bcType=3 -noplot -f=fsolveLinearcoupledSystemTestFreeG8;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=160 -ny=160 -bcType=3 -noplot -f=fsolveLinearcoupledSystemTestFreeG16;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=320 -ny=320 -bcType=3 -noplot -f=fsolveLinearcoupledSystemTestFreeG32;
runShell -case=coupledSystem -solver=fsolve -funcDefFile=coupledSystemTestFuncDef -nx=640 -ny=640 -bcType=3 -noplot -f=fsolveLinearcoupledSystemTestFreeG64;
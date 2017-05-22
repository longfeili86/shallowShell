#!/usr/bin/perl
#
# -- Longfei Li
use threads;
use Getopt::Long qw(GetOptions);


# default values
$rStart=1;
$rEnd=6;
$run="convRate"; #convRate or bifurcation

GetOptions('rStart=i'=>\$rStart,'rEnd=i'=>\$rEnd,'run=s'=>\$run);
print("$run\n");

# convRate commands
@convCmds = (
    #biharmonic conv test
    "convRate -case=biharmonic -test=biharmTrigTest -run -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=biharmonic -test=biharmPolyTest -run -noplot -rStart=$rStart -rEnd=$rEnd",
    #nonlinear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    #"convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    #"convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    #linear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    #"convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    #"convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd"
    );

# bifurcation commands
$nx=320; 
$ny=320;
@bifCmds= (
    "bifurcationRun -bcType=1 -xiMin=-2000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExSupported",
    "bifurcationRun -bcType=2 -xiMin=-8000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExClamped",
    "bifurcationRun -bcType=3 -xiMin=-3000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExFree",
    "bifurcationRun -bcType=4 -xiMin=-3500 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCS",
    "bifurcationRun -bcType=5 -xiMin=-4000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCF"
);

# non-uniform thermal loading
$nx=320;
$ny=320;
@nonuniformTLCmds=(
"runShell -case=coupledSystem -nonlinear -saveIC -solver=fsolve -tol=1.000000e-05 -bcType=1 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFsolveSupported",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=fsolve -tol=1.000000e-05 -bcType=2 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFsolveClamped",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=fsolve -tol=1.000000e-05 -bcType=3 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFsolveFree",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=fsolve -tol=1.000000e-05 -bcType=4 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFsolveCS",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=fsolve -tol=1.000000e-05 -bcType=5 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFsolveCF",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=newton -tol=1.000000e-05 -bcType=1 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLNewtonSupported",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=newton -tol=1.000000e-05 -bcType=2 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLNewtonClamped",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=newton -tol=1.000000e-05 -bcType=3 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLNewtonFree",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=newton -tol=1.000000e-05 -bcType=4 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLNewtonCS",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=newton -tol=1.000000e-05 -bcType=5 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLNewtonCF",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=exPicard -tol=1.000000e-05 -bcType=1 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLExPicardSupported",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=exPicard -tol=1.000000e-05 -bcType=2 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLExPicardClamped",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=exPicard -tol=1.000000e-05 -bcType=3 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLExPicardFree",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=exPicard -tol=1.000000e-05 -bcType=4 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLExPicardCS",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=exPicard -tol=1.000000e-05 -bcType=5 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLExPicardCF",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=imPicard -tol=1.000000e-05 -bcType=1 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLImPicardSupported",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=imPicard -tol=1.000000e-05 -bcType=2 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLImPicardClamped",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=imPicard -tol=1.000000e-05 -bcType=3 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLImPicardFree",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=imPicard -tol=1.000000e-05 -bcType=4 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLImPicardCS",
"runShell -case=coupledSystem -nonlinear -saveIC -solver=imPicard -tol=1.000000e-05 -bcType=5 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLImPicardCF"
);


@thr;

$counter=0;
@cmds;
if($run eq "convRate")
{
    @cmds=@convCmds;
    
}
elsif($run eq "bifurcation")
{
    @cmds=@bifCmds;
}
elsif($run eq "nonuniformTL")
{
    @cmds=@nonuniformTLCmds;
}
else
{
    print("Error: unknown run! use either convRate or bifurcation\n");
    exit
}

foreach(@cmds)
{
    $counter++;
    $cmd="matlab -nodisplay -r \'tic;".$_.";toc;exit;\'>$run$counter.log";
    print("$cmd\n");
    push @thr,threads->create('msc', $cmd);
}



foreach(@thr)
{
    $_->join();
}

print("batchRun -run=$run exited normally!\n");

sub msc{ ## make system call
  system( @_ );
}





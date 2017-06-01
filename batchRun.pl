#!/usr/bin/perl
#
# -- Longfei Li
use threads;
use Getopt::Long qw(GetOptions);


# default values
$rStart=1;
$rEnd=6;
$run="convRate"; #convRate or bifurcation
$solver="";

GetOptions('rStart=i'=>\$rStart,'rEnd=i'=>\$rEnd,'run=s'=>\$run,'solver=s'=>\$solver);
print("$run\n");

# convRate commands
@convCmds = (
    #biharmonic conv test
    "convRate -case=biharmonic -test=biharmTrigTest -run -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=biharmonic -test=biharmPolyTest -run -noplot -rStart=$rStart -rEnd=$rEnd",
    #nonlinear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    #linear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd"
    );

# bifurcation commands
$nx=320; 
$ny=320;
@bifCmds= (
    "bifurcationRun -noplot -bcType=1 -xiMin=-2000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExSupported",
    "bifurcationRun -noplot -bcType=2 -xiMin=-6000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExClamped",
    "bifurcationRun -noplot -bcType=3 -xiMin=-3000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExFree",
    "bifurcationRun -noplot -bcType=4 -xiMin=-3500 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCS",
    "bifurcationRun -noplot -bcType=5 -xiMin=-4000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCF"
);

# bifurcation PAC commands
$nx=320; 
$ny=320;
@bifPACCmds= (
    "bifurcationPACRun -noplot -bcType=1 -xiMin=-2000 -dxi=50 -maxIter=10 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonSupported",
    "bifurcationPACRun -noplot -bcType=2 -xiMin=-6000 -dxi=50 -maxIter=10 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonClamped",
    "bifurcationPACRun -noplot -bcType=3 -xiMin=-3000 -dxi=50 -maxIter=10 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonFree",
    "bifurcationPACRun -noplot -bcType=4 -xiMin=-3500 -dxi=50 -maxIter=10 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonCS",
    "bifurcationPACRun -noplot -bcType=5 -xiMin=-4000 -dxi=50 -maxIter=10 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonCF"
);

# non-uniform thermal loading
$nx=320;
$ny=320;
@nonuniformTLCmds=(
"runShell -noplot -case=coupledSystem -nonlinear -saveIC -solver=$solver -tol=1e-6 -bcType=1 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLSupported_$solver",
"runShell -noplot -case=coupledSystem -nonlinear -saveIC -solver=$solver -tol=1e-6 -bcType=2 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLClamped_$solver",
"runShell -noplot -case=coupledSystem -nonlinear -saveIC -solver=$solver -tol=1e-6 -bcType=3 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLFree_$solver",
"runShell -noplot -case=coupledSystem -nonlinear -saveIC -solver=$solver -tol=1e-6 -bcType=4 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLCS_$solver",
"runShell -noplot -case=coupledSystem -nonlinear -saveIC -solver=$solver -tol=1e-6 -bcType=5 -nx=$nx -ny=$ny -maxIter=100 -funcDefFile=coupledSystemNonUniformTLFuncDef -f=nonuniformTLCF_$solver"
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
elsif($run eq "bifurcationPAC")
{
    @cmds=@bifPACCmds;
}
elsif($run eq "nonuniformTL")
{
    @cmds=@nonuniformTLCmds;
    $run.=$solver;
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





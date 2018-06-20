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
$nx=320; 
$ny=320;
$serial=0;

GetOptions('rStart=i'=>\$rStart,'rEnd=i'=>\$rEnd,'run=s'=>\$run,'solver=s'=>\$solver,'nx=i'=>\$nx,'ny=i'=>\$ny,'serial=i'=>\$serial);
print("$run\n");

# convRate commands
@convBiharmCmds = (
    #biharmonic conv test
    "convRate -case=biharmonic -test=biharmTrigTest -run -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=biharmonic -test=biharmPolyTest -run -noplot -rStart=$rStart -rEnd=$rEnd"
    );

@convNonlinearCmds=(
    #nonlinear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=5 -nonlinear",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=5 -nonlinear"
    );

@convLinearCmds=(
    #linear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=5",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=$rStart -rEnd=5"
    );

# bifurcation commands

@bifCmds= (
    "bifurcationRun -noplot -bcType=1 -xiMin=-2000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExSupported",
    "bifurcationRun -noplot -bcType=2 -xiMin=-6000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExClamped",
    "bifurcationRun -noplot -bcType=3 -xiMin=-3000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExFree",
    "bifurcationRun -noplot -bcType=4 -xiMin=-3500 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCS",
    "bifurcationRun -noplot -bcType=5 -xiMin=-4000 -dxi=50 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCF"
);

# bifurcation PAC commands
$gn=$nx/10;
@bifPACCmds= (
    "bifurcationPACRun -noplot -bcType=1 -xiMin=-2000 -dxi=50 -maxIter=30 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonG${gn}Supported",
    "bifurcationPACRun -noplot -bcType=2 -xiMin=-6000 -dxi=50 -maxIter=30 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonG${gn}Clamped",
    "bifurcationPACRun -noplot -bcType=3 -xiMin=-3000 -dxi=50 -maxIter=30 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonG${gn}Free",
    "bifurcationPACRun -noplot -bcType=4 -xiMin=-3500 -dxi=50 -maxIter=30 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonG${gn}CS",
    "bifurcationPACRun -noplot -bcType=5 -xiMin=-4000 -dxi=50 -maxIter=30 -tol=1e-6 -nx=$nx -ny=$ny -solver=newton -results=bifurcationPACnewtonG${gn}CF"
);

# non-uniform thermal loading

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
if($run eq "convRateBiharm")
{
    @cmds=@convBiharmCmds;
    
}
elsif($run eq "convRateNonlinear")
{
    @cmds=@convNonlinearCmds;
    
}
elsif($run eq "convRateLinear")
{
    @cmds=@convLinearCmds;
    
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
    if($serial==0)
    {# run in multithreads
	#print("multithreads\n");
	push @thr,threads->create('msc', $cmd);
    }
    else 
    {# run in serial
	#print("serial\n");
	system($cmd);
    }
}



foreach(@thr)
{
    $_->join();
}

print("batchRun -run=$run exited normally!\n");

sub msc{ ## make system call
  system( @_ );
}





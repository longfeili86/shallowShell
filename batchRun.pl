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
    "bifurcationRun -bcType=1 -xiMin=-2000 -dxi=10 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExSupported",
    "bifurcationRun -bcType=2 -xiMin=-8000 -dxi=10 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExClamped",
    "bifurcationRun -bcType=3 -xiMin=-3000 -dxi=10 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExFree",
    "bifurcationRun -bcType=4 -xiMin=-3500 -dxi=10 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCS",
    "bifurcationRun -bcType=5 -xiMin=-4000 -dxi=10 -maxIter=100 -tol=1e-6 -nx=$nx -ny=$ny -solver1=exPicard -solver2=imPicard -solver3=exPicard -results=bifurcationExImExCF"
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
else
{
    print("Error: unknown run! use either convRate or bifurcation\n");
    exit
}

foreach(@cmds)
{
    $counter++;
    $cmd="matlab -nodisplay -r \'".$_.";exit;\'>$run$counter.log";
    print("$cmd\n");
    push @thr,threads->create('msc', $cmd);
}



foreach(@thr)
{
    $_->join();
}


sub msc{ ## make system call
  system( @_ );
}





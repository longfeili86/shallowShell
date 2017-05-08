#!/usr/bin/perl
#
# -- Longfei Li
use threads;

$rEnd=6;
@cmds = (
    #biharmonic conv test
    "convRate -case=biharmonic -test=biharmTrigTest -run -noplot -rStart=1 -rEnd=$rEnd",
    "convRate -case=biharmonic -test=biharmPolyTest -run -noplot -rStart=1 -rEnd=$rEnd",
    #nonlinear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd -nonlinear",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd -nonlinear",
    #linear coupled system conv test
    "convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd",
    "convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -run  -noplot -rStart=1 -rEnd=$rEnd"
    );


@thr;

$counter=0;
foreach(@cmds)
{
    $counter++;
    $cmd="matlab -nodisplay -r \'".$_.";exit;\'>convRate$counter.log";
    print("$cmd\n");
    ## Do not use multithread. I don't have that much mem to do all the runs at
    ## the same time. Just do it one by one
    #push @thr,threads->create('msc', $cmd);
    msc($cmd);
}


## Do not use multithread. I don't have that much mem to do all the runs at
## the same time. Just do it one by one
# foreach(@thr)
# {
#     $_->join();
# }


sub msc{ ## make system call
  system( @_ );
}





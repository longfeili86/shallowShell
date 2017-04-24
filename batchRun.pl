#!/usr/bin/perl
#
# -- Longfei Li
use threads;


@cmds = (
    "convRate -case=biharmonic -test=biharmTrigTest -run -noplot -rStart=1 -rEnd=6",
    "convRate -case=biharmonic -test=biharmPolyTest -run -noplot -rStart=1 -rEnd=6"
    );


@thr;

$counter=0;
foreach(@cmds)
{
    $counter++;
    $cmd="matlab -nodisplay -r \'".$_.";exit;\'>convRate$counter.log";
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





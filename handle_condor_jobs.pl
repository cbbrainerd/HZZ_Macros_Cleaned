use strict;
use warnings qw/all/;
use feature qw/say/;

$,=' ';
my$files_per_cluster=$ARGV[0];
shift;

foreach my$dn (@ARGV) {
    -d $dn or die "Invalid directory $dn";
    my@fns=glob("${dn}/*.log");
    my$last_stage=0;
    my@successes;
    my@failures;
    my@unfinished;
    say "Analyzing directory $dn";
    foreach my$fn (@fns){
        my($stage,$cluster,$n)=my@finfo=$fn=~m@/([0-9]+).([0-9]+).([0-9]+).log$@ or die "$fn";
        next if $stage<$last_stage;
        $last_stage=$stage if $stage>$last_stage;
        open(my$fh,'<',$fn) or die;
        seek$fh,-1000,2;
        read$fh,$_,1000;
        if(m/return value \K[0-9]+/) {
            my$retval=$&;
            push@finfo,$retval;
            if($retval==0) {
                push@successes,\@finfo;
            } else {
                push@failures,\@finfo;
            }
        } else {
            push@unfinished,\@finfo;
        }
        close($fh);
    }
    my@failed;
    say 'Successful jobs:',scalar @successes;
    say 'Failed jobs:',scalar @failures;
    say 'Unfinished jobs:',scalar @unfinished;
    foreach my$info(@failures) {
        my($stage,$cluster,$n,$retval)=@$info;
        next if $stage<$last_stage;
        push@failed,$n;
    }
    foreach my$info(@unfinished) {
        my($stage,$cluster,$n)=@$info;
        next if $stage<$last_stage;
        push@failed,$n;
    }
    @failed=sort {$a<=>$b} @failed;
    my@inputs=glob("${dn}/inputs*_${last_stage}.txt");
    my$ofn=$inputs[0];
    $ofn=~s@[0-9]+(?=\.txt$)@1+${last_stage}@e or die;
    -f $ofn and die;
    open(OFH,'>',$ofn);
    @inputs==1 or die @inputs;
    open(my$ifh,'<',$inputs[0]) or die;
    #New filehandle to be read over, reset $. in preparation for loop
    $.=0;
    my$next_stage_run=0;
    foreach my$n(@failed) {
        #Find which files should be output
        my$low_bound=$files_per_cluster*$n+1;
        my$high_bound=$files_per_cluster*($n+1);
        #Read, discarding output until reaching the lower bound
        1 while(1+$.<$low_bound and defined(<$ifh>));
        while($.<$high_bound) {
        #Output to file until reaching the higher bound
            defined($_=<$ifh>) or last;
            print OFH;
            ++$next_stage_run;
        }
    }
    say "Rerunning over $next_stage_run files.";
    close(OFH);
    close($ifh);
}

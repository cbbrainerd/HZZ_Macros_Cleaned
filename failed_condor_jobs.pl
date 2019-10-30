use strict;
use warnings qw/all/;

my$job_dir=$ARGV[0] or die;

my @log_files=<${job_dir}/*.log>;
my @condor_files=<${job_dir}/condor_*.cfg>;

for my$fn(@log_files) {
    my$fh=open("$fn","<") or die;
    while<$fh> {

    }
}

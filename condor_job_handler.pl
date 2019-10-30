use strict;
use warnings qw/all/;
use feature qw/say/;

open(my$fh,'<',"$ARGV[0]") or die;
my @files;
while(<$fh>) {
    chomp;
    push @files,$_;
}
close$fh;

open($fh,'<','failed_job_numbers') or die;
while(<$fh>) {
    chomp;
    say $files[$_];
}

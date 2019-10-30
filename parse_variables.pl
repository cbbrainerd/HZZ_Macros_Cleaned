use strict;
use feature qw/say/;
$|=1;
my $fourmuheader='HZZ4LeptonsAnalysis_4mu.h';
my @filenames=<*.C>;
push(@filenames,$fourmuheader);

open(my $fhl,'<','list_of_variables') or die "Failed to open.";
open(my $fhh,'<',$fourmuheader) or die "Failed to open.";
open(FHU,'>','list_of_unused_variables') or die "Failed to open.";
open(FHG,'>','list_of_used_variables') or die "Failed to open.";
my @variables=();

while(<$fhl>){
    chomp;
    my @var_type=split/\s+/;
    $_=$var_type[1];
    s/^\**//;
    s/[[;].*//;
    my $var=$_;
    my $regex=qr/(?<!["A-Za-z0-9_])$var(?!["A-Za-z0-9_])/;
    my $used=0;
    my $comment_block=0;
    for my$fn(@filenames){
        open(my $fh,'<',$fn) or die "Failed to open \"$fn\".";
        while(<$fh>){
            if(!$comment_block) {
                s!//.*!!; #Remove comments
                s!/\*.*?\*/!!; #Remove inline block comments
                $comment_block=s!/\*.*!!; #Detect start of block comment
            } else {
                $comment_block=!s!.*?\*/!!; #Detect end of comment block
                s!.*!! if $comment_block; #Still in comment block
            }
            $used+=()=m/$regex/g;
            #last if $used;
        }
        close($fh);
        #last if $used;
    }
    say "$var: $used" if $used;
    $_=$var;
    say FHU unless ($used-2) or /^b_/;
    say FHG if ($used-2) and !/^b_/;
}

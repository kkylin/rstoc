#!/usr/bin/env perl

foreach $file (@ARGV) {
    $_ = `cat $file`;

    if ( /Copyright \(C\).*by the authors/gs ) {
    	s/(Copyright \(C\).*?) by the authors/$1,2017 by the authors/gs;

    	$outfile = "${file}.out";
    	open(FOUT, ">$outfile");
    	print FOUT "$_";
    	close(FOUT);
    	system("mv -f $outfile $file");
    }

    # if ( /bruce\@math.duke.edu/gs ) {
    # 	s/bruce\@math.duke.edu/bruce.warren.rogers\@gmail.com/gs;

    # 	$outfile = "${file}.out";
    # 	open(FOUT, ">$outfile");
    # 	print FOUT "$_";
    # 	close(FOUT);
    # 	system("mv -f $outfile $file");
    # }
}

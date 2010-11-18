#!/usr/bin/perl
#
# written by Kai Wang (kw2110@columbia.edu)
#
use warnings;
use strict;

my $cvs_dir = '/users/kw2110/cvs/ARACNE/project/ARACNE';
my $work_dir = '/users/kw2110/project/_ARACNE';

print "Testing ARACNE code ... \n";

if ( !(-e $cvs_dir."/aracne") ) {
	print "Please compile the CVS ARACNE code first!\n";
	exit(1);
}

if ( !(-e $work_dir."/aracne") ) {
    print "Please compile the wokring ARACNE code first!\n";
	exit(1);
}

print "==> Testing full matrix construction ... ";
`$cvs_dir/aracne -i $cvs_dir/test/arraydata100x336.exp -n 0.0121 -o $cvs_dir/test/1.adj`;
`$work_dir/aracne -i $work_dir/test/arraydata100x336.exp -n 0.0121 -o $work_dir/test/1.adj`;
my $rslt = `diff -I "^>" $cvs_dir/test/1.adj $work_dir/test/1.adj`;

if (!$rslt) {
    print "Success!\n";
}
else {
	print "Fail!\n";
	exit(1);
}

print "==> Testing sub-matrix construction ... ";
`$cvs_dir/aracne -i $cvs_dir/test/arraydata100x336.exp -h 25 -k 0.15 -o $cvs_dir/test/2.adj`;
`$work_dir/aracne -i $work_dir/test/arraydata100x336.exp -h 25 -k 0.15 -o $work_dir/test/2.adj`;
$rslt = `diff -I "^>" $cvs_dir/test/2.adj $work_dir/test/2.adj`;

if (!$rslt) {
    print "Success!\n";
}
else {
	print "Fail!\n";
	exit(1);
}

print "==> Testing conditional sub-matrix construction ... ";
`$cvs_dir/aracne -i $cvs_dir/test/arraydata100x336.exp -h 50 -c +75 0.35 -k 0.15 -o $cvs_dir/test/3.adj`;
`$work_dir/aracne -i $work_dir/test/arraydata100x336.exp -h 50 -c +75 0.35 -k 0.15 -o $work_dir/test/3.adj`;
$rslt = `diff -I "^>" $cvs_dir/test/3.adj $work_dir/test/3.adj`;

if (!$rslt) {
    print "Success!\n";
}
else {
	print "Fail!\n";
	exit(1);
}

print "==> Testing reconstruction with annotation information ... ";
`printf "G1\nG4\nG7\n" > $work_dir/test/tf_list.dat`;
`$cvs_dir/aracne -i $cvs_dir/test/arraydata10x336.exp -l $work_dir/test/tf_list.dat -k 0.16 -t 0.05 -e 0.1 -o $cvs_dir/test/4.adj`;
`$work_dir/aracne -i $work_dir/test/arraydata10x336.exp -l $work_dir/test/tf_list.dat -k 0.16 -t 0.05 -e 0.1 -o $work_dir/test/4.adj`;
$rslt = `diff -I "^>" $cvs_dir/test/4.adj $work_dir/test/4.adj`;

if (!$rslt) {
    print " Success!\n";
}
else {
	print " Fail!\n";
	exit(1);
}

`rm $work_dir/test/tf_list.dat`;
`rm $work_dir/test/1.adj`;
`rm $work_dir/test/2.adj`;
`rm $work_dir/test/3.adj`;
`rm $work_dir/test/4.adj`;

`rm $cvs_dir/test/1.adj`;
`rm $cvs_dir/test/2.adj`;
`rm $cvs_dir/test/3.adj`;
`rm $cvs_dir/test/4.adj`;


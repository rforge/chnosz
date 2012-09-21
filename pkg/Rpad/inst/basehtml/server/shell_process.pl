#!/usr/bin/perl -w
#!c:/apps/perl/bin/perl.exe 

# The following line is a test script to see if it works.
#http://localhost/Rpad/server/shell_process.pl?&ID=testtesttest&shell_commands=ls
use strict ;
use CGI qw/:standard send_http_header/;
use utf8; # the decoding stuff is from spell-check-logic.cgi in the htmlarea spellcheck plugin
use Encode;

binmode STDIN, ':utf8';
binmode STDOUT, ':utf8';

my $dir = $ENV{SCRIPT_FILENAME};
$dir =~ s/shell_process.pl//; # strip off the perl file name (but leave the trailing slash)
chdir $dir;  

my $Rpad_ID = param('ID');
chomp($Rpad_ID);
length($Rpad_ID) == 12 or die "bad directory $Rpad_ID";
chdir($Rpad_ID);  

my $commands = decode('UTF-8', param('shell_commands')); 
$commands =~ s/\xA0/ /g; #replace non-breaking spaces with regular spaces
open(my $filehandle, '>exec.sh');
print $filehandle $commands;
close $filehandle;
chmod 0777, './exec.sh';
my $ret = `./exec.sh`; # backticks do an exec
# print # send_http_header('text/plain');
print "Content-type: text/plain\n\n";
print $ret;

#!/usr/bin/perl -w
#!c:/apps/perl/bin/perl.exe 

# The following line is a test script to see if it works.
# http://localhost/Rpad/server/R_process.pl?&ID=ddNTlmHSvWZF&command=R_commands&R_commands=print('hello')
use Statistics::Rpad ;
use strict ;
use CGI qw/:standard send_http_header/;
use Cwd ;
use utf8; # the decoding stuff is from spell-check-logic.cgi in the htmlarea spellcheck plugin
use Encode;

CGI::initialize_globals();

binmode STDIN, ':utf8';
binmode STDOUT, ':utf8';

my $dir = $ENV{SCRIPT_FILENAME}; # works for apache 1.3 & 2.0
$dir =~ s/[Rr]_process.pl//; # strip off the perl file name (but leave the trailing slash)
if ($dir eq '') {
  my $dir = $ENV{PATH_TRANSLATED}; # works for IIS
  $dir =~ s/[Rr]_process.pl//; # strip off the perl file name (but leave the trailing slash)
}

my $Rpad_ID = param('ID');
chomp($Rpad_ID);
length($Rpad_ID) == 12 or die "bad directory $Rpad_ID";
my $p_command = param('command');
chomp $p_command ;
my $R = Statistics::Rpad->new(
#                           log_dir => '/var/www/Rpad/server/' . $Rpad_ID
                           log_dir => $dir . $Rpad_ID
                           ) ;

$R->start_sharedR;

my $ret = "";

if ($p_command eq 'logout') {
  $R->stopR();
}
elsif ($p_command eq 'R_commands') {
  my $R_commands = decode('UTF-8', param('R_commands')); 
  $R_commands =~ s/\xA0/ /g; #replace non-breaking spaces with regular spaces
  $R->send($R_commands) ; 
  $ret = $R->read ;
}
#print send_http_header('text/plain');
#print "Content-type: text/plain\n\n";
print header;
print $ret . " ";  


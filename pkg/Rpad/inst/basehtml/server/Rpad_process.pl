#!/usr/bin/perl -w
#!c:/apps/perl/bin/perl.exe 

# The following line is a test script to see if it works.

#http://localhost/Rpad/server/Rpad_process.pl?command=login
#http://localhost/Rpad/server/Rpad_process.pl?command=savefile&ID=testtesttest&filename=test.txt&content=hellothere
use strict ;
use CGI qw/:standard/;
use File::Path;
use utf8; # the decoding stuff is from spell-check-logic.cgi in the htmlarea spellcheck plugin
use Encode;

CGI::initialize_globals();

binmode STDIN, ':utf8';
binmode STDOUT, ':utf8';

my $dir = $ENV{SCRIPT_FILENAME};
$dir =~ s/Rpad_process.pl//; # strip off the perl file name (but leave the trailing slash)
chdir $dir;  

my $p_command = param('command');
chomp($p_command);

print STDERR $p_command;

if ($p_command eq "login") {
  # make a randomly-named directory for the user
  my $max   = 10;
  my @chars = ('A'..'Z', 'a'..'z');
  my $str   = '';
  $str .= $chars[rand @chars] for 1 .. $max;
  my $dirname = "dd" . $str;
  mkdir $dirname;
  # return an ID to the user
  # print send_http_header('text/plain');
  print "Content-type: text/plain\n\n";
  print $dirname;
}
elsif ($p_command eq "logout") {
  # delete the directory
  my $dirname = param('ID');
  chomp $dirname;
  rmtree $dirname;
# print # send_http_header('text/plain');
  print "Content-type: text/plain\n\n";
  print 'XXLogged out'; # for some reason, need to return something here or IE hangs sometimes
}
elsif ($p_command eq "savefile") {
  my $dirname = param('ID');
  chomp $dirname ;
  chdir $dirname ;  
  my $filename = param('filename');
  chomp $filename ;
  my $content = decode('UTF-8', param('content')); 
  $content =~ s/\xA0/ /g; #replace non-breaking spaces with regular spaces
  open(my $f, ">$filename") ;
  print $f $content ;
  close $f ;
# print # send_http_header('text/plain');
  print "Content-type: text/plain\n\n";
  print 'File saved'; # for some reason, need to return something here or IE hangs sometimes
}
elsif ($p_command eq "savefileinbase") {
  my $filename = param('filename');
  chomp $filename ;
  my $content = decode('UTF-8', param('content')); 
  $content =~ s/\xA0/ /g; #replace non-breaking spaces with regular spaces
  open(my $f, ">../$filename") ;
  print $f $content ;
  close $f ;
# print # send_http_header('text/plain');
  print "Content-type: text/plain\n\n";
  print 'File saved'; # for some reason, need to return something here or IE hangs sometimes
}
elsif ($p_command eq "savehtml") {
  my $header = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN"><html><head><meta http-equiv="Content-Type" content="text/html; charset=utf-8" /><script type="text/javascript" src="editor/Rpad_loader.js"></script></head><body>';
  my $footer = '</body></html>';
  my $filename = param('filename');
  chomp $filename ;
  my $content = decode('UTF-8', param('content')); 
  $content =~ s/\xA0/ /g; #replace non-breaking spaces with regular spaces
  chomp $content ;
  open(my $f, ">../$filename")
      or print "Content-type: text/plain\n\nRpad ERROR: couldn't open file on the server: it may be write protected." and exit;  
  print $f $header . $content . $footer  
      or print "Content-type: text/plain\n\nRpad ERROR: file not saved" and exit;
  close $f or print "Content-type: text/plain\n\nRpad ERROR: file not saved" and exit;
# print # send_http_header('text/plain');
  print "Content-type: text/plain\n\n";
  print 'HTML file saved'; # for some reason, need to return something here or IE hangs sometimes
}

#!c:/apps/perl/bin/perl.exe 
#!/usr/bin/perl -w


# This receives and stores a file uploaded from Rpad.
# It is used with the fileupload.R on the browser side.

# The process is fairly convoluted. Xmlhttp can't be used to upload files,
# so we have to use a hidden iframe to "post" the uploaded file.
# When the file has finished uploading, perl sends back a snippet
# of javascript that turns around and executes an R command (normally to
# load or do some manipulation on the uploaded file).
# For an even fancier file upload widget, see
# http://blog.joshuaeichorn.com/archives/2005/05/01/ajax-file-upload-progress/

#use strict ;
use CGI qw/:standard send_http_header/;

my $filename = param('uploadfilename');
my @t = split('\\\\|/', $filename);
my $basefilename =  pop(@t);

my $dir = $ENV{SCRIPT_FILENAME};
$dir =~ s/\/fileupload.pl//; # strip off the perl file name (but leave the trailing slash)
my $subdir = param('activedir');
chdir "$dir/$subdir";  
print "$dir/$subdir";  
print "$filename<br>";

open (OUTFILE,">$basefilename");
binmode OUTFILE;
while ($bytesread=read($filename,$buffer,1024)) {
    print OUTFILE $buffer;
}
close OUTFILE;

my $Rresponse = param('response');
print "Content-type: text/html\n\n";
print "<script type='text/javascript'>top.R_run_commands('RresponseToUpload(\"$dir/$subdir/$basefilename\")')</script>";



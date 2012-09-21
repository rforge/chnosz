
# p5httpd: Tiny HTTP server written in perl, in order to run perl CGI
# scripts on EPOC machines. May be useful for other purposes.  tested
# on EPOC and linux

# POD documentation at end of file

################ Configuration  ##################################

# Change the next line to let p5httpd know where it lives:
# This directory will also be the server root directory

our $version=0.04;
our $basdir = "c:/e/www/rpad/inst/basehtml/";

# The port on which p5httpd will listen:

our $port = 80;

# Add your own MIME types here; text/plain is the default.
our %mime_types = (
    '\.txt'  =>  'text/plain',
    '\.htm'  =>  'text/html',
    '\.html' =>  'text/html',
    '\.gif'  =>  'image/gif',
    '\.jpg'  =>  'image/jpeg',
    '\.png'  =>  'image/png',
    '\.xbm'  =>  'image/x-xbitmap',
    '\.css'  =>  'text/css',
    '\.js'   =>  'application/x-javascript',
    '\.htc'  =>  'text/x-component',
    '\.xml'  =>  'text/xml',
    '\.pdf'  =>  'application/pdf',
    '\.eps'  =>  'application/postscript',
    '\.ps'   =>  'application/postscript',
    '\.Rpad' =>  'text/html',
    '\.pl'   =>  ''
		  );


###################  no real need to edit below ##################

require 5.6.0;
package p5httpd; # keep namespace separate from CGI scripts

use Socket;
use strict;

our($localname, $EPOC);

initialise();

main_loop();

################################## Subroutines ###################

sub logerr($$); sub logmsg($); sub cat($$;$); # forward declarations


sub initialise {
    # Do we run under EPOC?
    $EPOC = ($^O =~ /epoc/i ? 1 : 0);

    #Test the server on some non-EPOC machine:
    #$basdir = "/usr/local/home/hlub/projekten/webplek/rlwrap/p5httpd-$version" unless $EPOC;
    #$basdir = "/dolomieten/var/www/html" unless $EPOC;
    #$port = 4711 unless $EPOC;


   -d $basdir or die "$basdir is not a directory\n you should probably edit p5httpd.\n";
   $p5httpd::basdir= $basdir; # make this variable visible for CGI scripts
}


sub main_loop {
  my $tcp = getprotobyname('tcp');
  socket(Server, PF_INET, SOCK_STREAM, $tcp)      or die "socket: $!";
  unless ($EPOC) { # to prevent "address in use" errors (not in EPOC (and Win?))
    setsockopt(Server, SOL_SOCKET, SO_REUSEADDR, pack("l", 1))
	  or warn "setsockopt: $!";
  }
  bind(Server, sockaddr_in($port, INADDR_ANY))    or die "bind: $!";
  listen(Server,SOMAXCONN)                        or die "listen: $!";
  logmsg "server started on port $port";

 CONNECT:
  for ( ; accept(Client,Server); close Client) {

    *STDIN = *Client;
    *STDOUT = *Client;

    my $remote_sockaddr  = getpeername(STDIN);
    my (undef, $iaddr)   = sockaddr_in($remote_sockaddr);
    my $peername         = gethostbyaddr($iaddr, AF_INET) || "localhost";
    my $peeraddr         = inet_ntoa($iaddr) || "127.0.0.1";

    my $local_sockaddr   = getsockname(STDIN);
    my (undef, $iaddr)   = sockaddr_in($remote_sockaddr);
    $localname           = gethostbyaddr($iaddr, AF_INET) || "localhost";
    my $localaddr        = inet_ntoa($iaddr) || "127.0.0.1";


    chomp($_ = <STDIN>);
    my ($method, $url, $proto, undef) = split;

    $url =~ s#\\#/#g;
    logmsg "<- $peername: $_";
    my ($file, undef, $arglist) = ($url =~ /([^?]*)(\?(.*))?/); # split at ?
    my $file_escaped = $file;
    $file =~ s/%([\dA-Fa-f]{2})/chr(hex($1))/eg; # %20 -> space

    if ( $method !~ /^(GET|POST|HEAD)$/ ) {
      logerr 400, "I don't understand method $method";
      next CONNECT;
    }


    if (-d  "$basdir/$file" ) {
      unless ($file =~ m#/$#) {
	   redirect ("$file/");
	   next CONNECT;
	  }
	my $dir = "$basdir/$file";
    if (-f "$dir/index.html") {
	  $file .= "/index.html";
    } else {
	  directory_listing($file);
	  next CONNECT;
      }
    }

    if ( not -r "$basdir/$file" ) {
      logerr 404, "$file: $!";
      next CONNECT;
    }

    my @env_vars = qw(USER_AGENT CONTENT_LENGTH CONTENT_TYPE);
    foreach my $var (@env_vars) {
		$ENV{$var} = ""; # delete $ENV{$var} will crash perl on Netbook :-(
		}
    while(<STDIN>) {
	  s/[\r\l\n\s]+$//;
	  /^User-Agent: (.+)/i and $ENV{USER_AGENT} = $1;
	  /^Content-length: (\d+)/i and $ENV{CONTENT_LENGTH} = $1;
	  /^Content-type: (.+)/i    and $ENV{CONTENT_TYPE} = $1;
	  last if (/^$/);
	}


    print "HTTP/1.1 200 OK\n"; # probably OK by now



    if ( $file =~ m/\.pl$/i) {
      $ENV{SERVER_PROTOCOL} = $proto;
      $ENV{SERVER_PORT}     = $port;
      $ENV{SERVER_NAME}     = $localname;
      $ENV{SERVER_URL}      = "http://$localname:$port/";
      $ENV{SCRIPT_NAME}	    = $file;
      $ENV{SCRIPT_FILENAME} = "$basdir/$file";
      $ENV{REQUEST_URI}     = $url;
      $ENV{REQUEST_METHOD}  = $method;
      $ENV{REMOTE_ADDR}     = $peeraddr;
      $ENV{REMOTE_HOST}     = $peername;
      $ENV{QUERY_STRING}    = $arglist;
      $ENV{SERVER_SOFTWARE} = "p5httpd/$version";

      if ($method =~ /POST/) {
	    logmsg "<- Content-length: $ENV{CONTENT_LENGTH}, type: $ENV{CONTENT_TYPE}";
      }
      cgi_run ($file,$arglist);
      next CONNECT;
    }

    my $mime_type =  "text/html"; # default
    foreach my $suffix (keys %mime_types) {
    	if ($file =~ /$suffix$/i) {
	    $mime_type = $mime_types{$suffix};
	    last;
	}
    }
    cat $file, $mime_type, $method || logerr 500, "$file: $!";
      next CONNECT;
  }
    die "Fatal error: accept failed: $!\n"; # This should never happen
  }

#################### other subroutines ####################

sub logmsg ($) {
    my $fulltime = localtime();
    my ($hms) = ($fulltime =~ /(\d\d:\d\d:\d\d)/);
    print STDERR  "$$ $hms @_\n";
}

sub logerr ($$) {
  my ($code, $detail) = @_;
  my %codes =
      ( 200  => 'OK',
        400  => 'Bad Request',
        404  => 'Not Found',
        500  => 'Internal Server Error',
        501  => 'Not Implemented',
      );
  my $msg = "$code " . $codes{$code};
  logmsg "-> $msg $detail";
  print  <<EOF;
HTTP/1.0 $msg
Content-type: text/html

<html><body>
<h1>$msg</h1>
<p>$detail</p>
<hr>
<p><I>p5httpd/$version server at $localname port $port</I></p>
</body></html>
EOF
}



sub cat($$;$){   # cat ($file, $mimetype) writes Content-type header and $file to STDOUT
  my ($file, $mimetype, $method) = @_;
  $method = "GET" unless $method;
  my $fullpath = "$basdir/$file";

  my ($x,$x,$x,$x,$x,$x,$x,$length,$x,$mtime) = stat($fullpath);
  $mtime = gmtime $mtime;
  my ($day, $mon, $dm, $tm, $yr) =
	  ($mtime =~ m/(...) (...) (..) (..:..:..) (....)/);

  print "Content-length: $length\n";
  print "Last-Modified: $day, $dm $mon $yr $tm GMT\n";
  print "Content-type: $mimetype\n\n";
  my $sent=0;
  if ($method eq "GET") {
    local $/; undef $/; # gobble whole files
    open IN, "<$fullpath" || return 0;
    binmode IN;
    my $content = <IN>;
    close IN;
    $sent = length($content);
    print $content;
  }
  logmsg "-> 200 OK $file: $sent bytes sent as $mimetype";
  return 1;
}




sub cgi_run {
  my ($file,$arglist) = @_;
  my ($dir) = ($file =~ /^(.*\/)/);
  my $path = "$basdir/$file";
#  chdir "$basdir/$dir" or return logerr 500, "Cannot chdir to $basdir/$dir: $!";
  $path=~s/[A-Z]://;
  logmsg "-> doing $path";

  {
    package main;
sub header {
    return "";
}
print STDERR header();

    do $path;
  } # or return logerr 500, "Cannot do $path: $!";
  $@ and logerr 500, "in $file:<br>  $@";
#  chdir $basdir;
}


sub directory_listing {
  my ($dir) = @_;
  $dir =~ s#//#/#g;
  chdir "$basdir/$dir" or return logerr 500, "Cannot chdir to $basdir/$dir: $!";
  my @files = glob("*");
  print  <<EOF;
HTTP/1.0 200 OK
Content-type: text/html

<html>
<head><title>$dir</title></head>
<body>
<h1>$dir</h1>
EOF
  print "<p><a href=..>..</a></p>\n";
  foreach my $file (sort @files) {
    -d $file and $file .= "/";
    print "<p><a href=./$file>$file</a></p>\n";
  }
  print <<EOF;
<hr>
<p><I>p5httpd/$version server at $localname port $port</I>
</body></html>
EOF
logmsg "-> 200 OK listing $dir";
}


sub redirect {
	my ($redir) = @_;
	print "HTTP/1.0 301 Moved Permanently\nLocation: $redir\n\n";
    logmsg "-> 301 Moved Permanently to $redir"
}



__END__


=head1 NAME

p5httpd - tiny perl http server

=head1 SYNOPSIS

 perl path/to/p5httpd

=head1 DESCRIPTION

Tiny HTTP server written in perl. Written in order to run perl CGI
scripts on EPOC machines. Tested on EPOC and linux

=head1 USING CGI SCRIPTS WITH p5httpd

As EPOC doesn't have the C<fork()> system call, CGI scripts are run by
the same process as the server (to be precise, they are
evaluated by a C<do $filename> command).

This makes writing CGI scripts a little unusual in the following ways:

=over 4

=item 1

Obviously, a CGI script has to be a perl program. B<p5httpd>
cannot run python or OPL scripts.

=item 2


If you call C<exit()> in your script, the whole server will quit
(C<die()> will just print an error message). CGI scripts may hang, or
even crash, the server. Filehandles remain open across invocations.


=item 3

All CGI scripts run in the same namespace (package C<main> by default).
Package globals remain defined across invocations. This may be useful
in some situations, but you should be aware of it.

=item 4

If you use the C<CGI.pm> package, you have to use the (undocumented)
subroutine C<CGI::initialize_globals()> to get it to re-read the script
parameters.  If you don't, the script will only read them the first
time it runs.

=back

The server does a C<chdir> to a script's directory before running it.
The name of the server root directory is given by C<$p5httpd::basdir>.

The server uses a few small subroutines that can be usefully called
from CGI scripts:

C<p5httpd::cat($file, $mimetype)> prints a C<Content-type: $mimetype> header
and dumps C<$file>.

C<p5httpd::logerr($code, $detail)> prints out an HTTP error code (like
C<404 not found>) and the message C<$detail>.



=head1 CREDITS

Based on phttpd (pure perl httpd, (C) Paul Tchistopolskii 1998, 99
(http://www.pault.com, e-mail: paul@qub.com) Adapted and extended for
use under Olaf Flebbe's port of Perl 5.6.0 to EPOC
(http://www.linuxstart.com/~oflebbe/perl/perl5.html) by Hans Lub
(hlub@knoware.nl)


=cut

# Local Variables:
# mode: cperl
# End:














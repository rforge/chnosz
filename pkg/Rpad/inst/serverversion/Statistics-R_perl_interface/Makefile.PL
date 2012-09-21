###############
# MAKEFILE.PL #
###############

use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'          => 'Statistics::Rpad' ,
    'VERSION_FROM'  => 'lib/Statistics/Rpad.pm' ,
    #'PREREQ_PM'     => {'Foo'  => 0.01} ,
    'EXE_FILES'     => [ qw(script/statistics-r.pl) ],
    ($] >= 5.005 ?
      ( ABSTRACT_FROM => 'lib/Statistics/Rpad.pm',
        AUTHOR        => 'Graciliano M. P. <gm@virtuasites.com.br>'
       ) : ()
     ),
);

1;


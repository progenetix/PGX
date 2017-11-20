package PGX::GenomePlots::PlotParameters;

use Data::Dumper;
use YAML::XS qw(LoadFile DumpFile);

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  read_plot_defaults
  args_modify_plot_parameters
  hex2rgb
);

################################################################################

sub read_plot_defaults {

  use File::Basename;
  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $plotPars  =   LoadFile($path_of_this_module.'/../rsrc/config/plotdefaults.yaml');
  return  $plotPars;

}

################################################################################

sub args_modify_plot_parameters {

=pod

Expects:
  - the dash prefixed input args
  - the plot parameter object

Returns:
  - modified plot parameter object

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my ($plotPars, $args) =   @_;

  # local defaults overwrite general values; but command line parameters
  # take precedence later on
  my $locDefaults       =   {};
  my $defaultsDir;

  if ($args->{'-defaultsfile'} =~ /\w+?\.\w+?$/) {
    $defaultsDir        = $args->{'-defaultsfile'};
    $defaultsDir        =~  s/\/[\w\.\,]+?$//;
  }

  if (-f $args->{'-defaultsfile'}) {
    $locDefaults        =   LoadFile($args->{'-defaultsfile'});
    foreach my $par (keys %$plotPars) {
      if ($locDefaults->{$par}) {
        $plotPars->{$par}   =   $locDefaults->{$par};
      }
    }
  }

  # the -plottype | -colorschema specific mappings are processed first & removed
  # thereafter (there are still fallbacks if no parameters given)
  if (grep{ $args->{'-colorschema'} eq $_ } keys %{ $plotPars->{colorschemas} }) {
    my $colorschema     =   $args->{'-colorschema'};
    foreach (keys %{ $plotPars->{colorschemas}->{ $colorschema } }) {
      if ($plotPars->{colorschemas}->{ $colorschema }->{$_} =~ /^(\#\w{6})$/) {
        $plotPars->{$_} =   $1 }
    }
    delete $plotPars->{colorschemas};
    delete $args->{'-colorschema'};
  }

  # adjusting arguments for the selected plot type
  if (grep{ $args->{'-plottype'} eq $_ } keys %{ $plotPars->{plottype_values} }) {
    foreach (keys %{ $plotPars->{plottype_values}->{ $args->{'-plottype'} } }) {
      $plotPars->{$_} =   $plotPars->{plottype_values}->{ $args->{'-plottype'} }->{$_} }
    delete $plotPars->{plottype_values};
    delete $args->{'-plottype'};
  }

  # arguments to parameters 
  foreach my $par (keys %$plotPars) {

    if ($args->{'-'.$par} !~ /\w/) { next }

    # special evaluation: regions
    if ($par eq 'plotregions') {
      foreach my $plotregion (split(',', $args->{'-plotregions'})) {
        if ($plotregion =~ /^(?:chro?)?(\w\d?)\:(\d+?)\-(\d+?)$/) {
          my $plotR =   {
            reference_name  =>  $1,
            start           =>  $2,
            end             =>  $3,
          };
          push(@{ $plotPars->{'plotregions'} }, $plotR);
    }}}
    
    # special evaluation: markers
    elsif ($par eq 'markers') {
      foreach my $plotregion (split(',', $args->{'-markers'})) {
        if ($plotregion =~ /^(?:chro?)?(\w\d?)\:(\d+?)\-(\d+?)(?:\:([\w \-]+?))?(?:\:(\#\w\w\w(?:\w\w\w)?))?$/) {
           my $mark =   {
            reference_name  =>  $1,
            start           =>  $2,
            end             =>  $3,
            label           =>  $4,
            color           =>  $5,
          };
          if ($mark->{color} !~ /^\#\w\w\w(?:\w\w\w)?$/) {
            $mark->{color}  =   random_hexcolor() }
          push(@{ $plotPars->{'markers'} }, $mark);
    }}}
    
    # list style parameters are provided comma concatenated => deparsed
    elsif (grep{ $par eq $_ } qw(chr2plot label_y_m)) {
      $plotPars->{$par}   =   [ split(',', $args->{'-'.$par}) ] }
    else {
      $plotPars->{$par}   =   $args->{'-'.$par} }

    # modified parameters stored in the local yaml defaults are collected here
    if (grep{ $_ eq $par } @{ $plotPars->{local_overrides} }) {
      $locDefaults->{$par}  =   $plotPars->{$par} }

  }

  # derived
  $plotPars->{pixyfactor}   =   1 * $plotPars->{size_plotarea_h_px} / (2 * $plotPars->{value_plot_y_max});
  foreach my $override (keys %$locDefaults) {
    if (! grep{ $_ eq $override } @{ $plotPars->{local_overrides} }) {
      delete $locDefaults->{$override};
    }
  }

  if (-d $defaultsDir) {
    DumpFile($args->{'-defaultsfile'}, $locDefaults) }

  return $plotPars;

}

################################################################################

sub hex2rgb {

    my ($r, $g, $b)     =   $_[0] =~  m/^\#?(\w{2})(\w{2})(\w{2})$/;

    return [ CORE::hex($r), CORE::hex($g), CORE::hex($b) ];

}

################################################################################

sub random_hexcolor {

  use List::Util qw( shuffle);

=pod

random_hexcolor()

=cut

  my @randRGB   =   (
    (shuffle(0..180))[0],
    (shuffle(60..255))[0],
    (shuffle(200..255))[0],
  );
  return('#'.sprintf("%x%x%x", shuffle(@randRGB)));

}


1;

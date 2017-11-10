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

  # the -do_plottype | -colorschema specific mappings are processed first & removed
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

  if (grep{ $args->{'-do_plottype'} eq $_ } keys %{ $plotPars->{plottype_values} }) {
    foreach (keys %{ $plotPars->{plottype_values}->{ $args->{'-do_plottype'} } }) {
      $plotPars->{$_} =   $plotPars->{plottype_values}->{ $args->{'-do_plottype'} }->{$_} }
    delete $plotPars->{plottype_values};
    delete $args->{'-do_plottype'};
  }

  foreach my $par (keys %$plotPars) {
    if ($args->{'-'.$par} =~ /^\-?\#?\w[\. \-\,\(\)\w\:]*?$/) {

    # list style parameters are provided comma concatenated
      if ($par eq 'plotregions') {
        if ($args->{'-'.$par} =~ /\w\:\d+?\-\d+?(?:\,|$)/) {
          foreach my $plotregion (split(',', $args->{'-'.$par})) {
            if ($plotregion =~ /^(?:chro?)?(\w\d?)\:(\d+?)\-(\d+?)$/) {
              push(
                @{ $plotPars->{$par} },
                {
                  reference_name  =>  $1,
                  start           =>  $2,
                  end             =>  $3,
                }
              );

      }}}}
      elsif (grep{ $par eq $_ } qw(chr2plot label_y_m)) {
        $plotPars->{$par}   =   [ split(',', $args->{'-'.$par}) ] }
      else {
        $plotPars->{$par}   =   $args->{'-'.$par} }
      if (grep{ $_ eq $par } @{ $plotPars->{local_overrides} }) {
        $locDefaults->{$par}  =   $plotPars->{$par};
      }
    }
  }

  # derived
  $plotPars->{pixyfactor}   =   1 * $plotPars->{size_plotarea_h_px} / (2 * $plotPars->{value_plot_y_max});
  foreach my $override (keys %$locDefaults) {
    if (! grep{ $_ eq $override } @{ $plotPars->{local_overrides} }) {
      delete $locDefaults->{$override};
    }
  }

  if (-d $defaultsDir) {
    DumpFile($args->{'-defaultsfile'}, $locDefaults);
  }

  return $plotPars;

}

################################################################################

sub hex2rgb {

    my ($r, $g, $b)     =   $_[0] =~  m/^\#?(\w{2})(\w{2})(\w{2})$/;

    return [ CORE::hex($r), CORE::hex($g), CORE::hex($b) ];

}

1;

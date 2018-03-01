package PGX::GenomePlots::PlotParameters;

use Data::Dumper;
use YAML::XS qw(LoadFile DumpFile);

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  read_plot_defaults
  args_modify_plot_parameters
  hex2rgb
  random_hexcolor
  frequencies2rgb
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
  my $defaultsDir       =   '';

  if (defined($args->{'-defaultsfile'})) {
    $defaultsDir        = $args->{'-defaultsfile'};
    $defaultsDir        =~  s/\/[\w\.\,]+?$//;
    if (-f $args->{'-defaultsfile'}) {
      $locDefaults        =   LoadFile($args->{'-defaultsfile'});
      foreach my $par (keys %$plotPars) {
        if ($locDefaults->{$par}) {
          $plotPars->{$par}   =   $locDefaults->{$par};
        }
      }
    }
  }

  # the -plottype | -colorschema specific mappings are processed first & removed
  # thereafter (there are still fallbacks if no parameters given)
  if (defined($args->{'-colorschema'}) && grep{ $args->{'-colorschema'} eq $_ } keys %{ $plotPars->{colorschemas} }) {
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
# print Dumper([$args->{'-plottype'}, $_, $plotPars->{plottype_values}->{ $args->{'-plottype'} }->{$_}]).'<hr/>';
      $plotPars->{$_} =   $plotPars->{plottype_values}->{ $args->{'-plottype'} }->{$_} }
    delete $plotPars->{plottype_values};
    delete $args->{'-plottype'};
  }

  # arguments to parameters
  foreach my $par (keys %$plotPars) {

    if (! defined($args->{'-'.$par}) || $args->{'-'.$par} !~ /\w/) { next }

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
      foreach (split(',', $args->{'-markers'})) {
        my @markervals  =   split(':', $_);
        if (
          $markervals[0]  =~ /^(chro?)?([\dxy]\d?)$/i
          &&
          $markervals[1]  =~ /^\d+?\-\d+?$/
        ) {
          my $mark     =   { reference_name  =>  $markervals[0] };
          $mark->{reference_name} =~  s/[^xy\d]//gi;
          ($mark->{start}, $mark->{end})  =   split('-', $markervals[1]);
          if ($markervals[2] =~ /^\w[\w \-\(\)\[\]]+?$/) {
            $mark->{label}  =  $markervals[2] }
          if ($markervals[3] =~ /^\#\w\w\w(\w\w\w)?$/) {
            $mark->{color}  =  $markervals[3] }
          if ($mark->{color} !~ /^\#\w\w\w(?:\w\w\w)?$/) {
            $mark->{color}  =   random_hexcolor() }
          push(@{ $plotPars->{'markers'} }, $mark);
    }}}

    # list style parameters are provided comma concatenated => deparsed
    elsif (grep{ $par eq $_ } qw(chr2plot label_y_m)) {
      $plotPars->{$par}   =   [ split(',', $args->{'-'.$par}) ] }
    else {
      $plotPars->{$par}   =   $args->{'-'.$par};
    }

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

  use List::Util qw(shuffle);

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

################################################################################

sub frequencies2rgb {

	my (
    $plotPars,
    $dupF,
    $delF,
    $maxF
  )             =   @_;
	if ($maxF < 0.001) {$maxF = 100}

	my $dupRGB		=		hex2rgb($plotPars->{color_var_dup_hex});
	my $delRGB		=		hex2rgb($plotPars->{color_var_del_hex});
  my @RGB;

  for my $i (0..2) {
    $dupRGB->[$i]	=		int($dupRGB->[$i] * $dupF / $maxF);
	  $delRGB->[$i]	=		int($delRGB->[$i] * $delF / $maxF);
    if (($dupRGB->[$i] + $delRGB->[$i]) < 255) { $RGB[$i] = $dupRGB->[$i] + $delRGB->[$i] }
    else { $RGB[$i] = 255 }
  }

	return	join(',', @RGB);

}



1;

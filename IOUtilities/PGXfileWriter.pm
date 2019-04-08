package PGX::IOUtilities::PGXfileWriter;

use Data::Dumper;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  write_svg
  write_probefile
  write_segmentfile
  write_labelfile
  write_status_matrix
  write_value_matrix
  write_frequency_matrix
);

################################################################################
################################################################################
################################################################################

sub write_svg {

  my $pgx       =   shift;
  my $fileName  =   shift;
 
  if (! -d $pgx->{parameters}->{path_loc}) {
    warn "No specification for output path $pgx->{parameters}->{path_loc} $!";
    return $pgx;
  }

  if ($svgName !~ /\.svg$/) {
    $fileName    =   $pgx->{parameters}->{plotid};    
    if ($fileName !~ /.../) {
      $fileName  =   "$^T"."$$" }
    $fileName    .=  '.svg';
  }
      
  $pgx->{svg_path_loc}  =   $pgx->{parameters}->{path_loc}.'/'.$fileName;

  if ($pgx->{parameters}->{path_web} =~ /./) {
    $pgx->{svg_path_web}    =   $pgx->{parameters}->{path_web}.'/'.$fileName;
    $pgx->{svg_html_embed}  =   '<a href="'.$pgx->{svg_path_web}.'" target="_blank"><img src="'.$pgx->{svg_path_web}.'" /></a>';
  }
  open  (FILE, ">", $pgx->{svg_path_loc}) || warn 'output file '.$pgx->{svg_path_loc}.' could not be created.';
  binmode(FILE, ":utf8");
  print FILE  $pgx->{svg};
  close FILE;
  
  return $pgx;

}

################################################################################

sub write_probefile {

=pod

Expects:
  - a list reference of genome position / value objects:
    [
      {
        no              =>  __integer__,          # 1 -> n
        probe_id        =>  __string__,
        reference_name  =>  __string__,
        position        =>  __integer__,
        value           =>  __long__,
      },
      {
      ...
      },
    ]

Returns:
 - a standard Progenetix style probe file

  ID  chro  pos log2
  cnvi0111187 17  35295593  0.0859121900
  cnvi0111188 8 65499402  -0.1438023000
  cnvi0111189 2 177061178 -0.0113166000
  cnvi0111190 5 70255894  0.0463862400
  ...

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No specification for output $file $!" }

  my @columnKs  =   qw(probe_id reference_name position value);
  open  FILE, '>'."$file";
  print FILE join("\t", @columnKs)."\n";
  foreach my $probe (@{ $pgx->{probedata} }) {
    print FILE join("\t", @{$probe}{@columnKs} )."\n";
  }
  close FILE;

}

################################################################################

sub write_segmentfile {

=pod

Expects:
  - a list reference of genome CNV objects:
    [
      {
        no              =>  __integer__,    # 1 -> n
        callset_id      =>  __string__,
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        variant_type    =>  __string__,     # DUP, DEL
        info            =>  {
          value           =>  __long__,
          svlen           =>  __integer__,
          probes          =>  __integer__,
          assembly_id     =>  __string__,     # GRCh36 ...
          experiment_type =>  __string__,     # aCGH ...
        },
      },
      {
      ...
      },
    ]

Returns:

  - a standard Progenetix segments  file

  sample  chro  start stop  mean  probes
  GSM481286 1 742429  7883881 -0.1594 699
  GSM481286 1 115673158 115705254 -0.3829 8
  GSM481286 1 115722621 119771659 0.167 424
  GSM481286 1 119776776 162617092 0.4168  1587
  GSM481286 1 162621657 165278686 0.6508  350
  GSM481286 1 165280711 167221337 0.4056  241
  GSM481286 1 167248788 168289603 0.6784  130
  ...


=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No specification for output $file $!" }

  my @columnKs  =   qw(sample chro start stop mean probes variant_type);

  if (! $pgx->{segmentdata}) {
    $pgx->{segmentdata} =  [ map{ @{ $_->{variants} } } @{ $pgx->{samples} } ] }

  open  FILE, '>'."$file";
  print FILE join("\t", @columnKs)."\n";
  foreach my $seg (@{ $pgx->{segmentdata} }) {
#print Dumper($seg).'<hr/>';
    print FILE join("\t",
      $seg->{callset_id},
      $seg->{reference_name},
      $seg->{start}->[0],
      $seg->{end}->[-1],
      $seg->{info}->{value},
      $seg->{info}->{probes},
      $seg->{info}->{variant_type},
    )."\n";
  }
  close FILE;

}################################################################################

sub write_labelfile {

=pod

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No specification for output $file $!" }

  my @columnKs  =   qw(sample key label color);

  open  FILE, '>'."$file";
  print FILE join("\t", @columnKs)."\n";
  foreach my $sample (@{ $pgx->{samples} }) {
    print FILE join("\t",
      $sample->{id},
      $sample->{sortkey},
      $sample->{sortlabel},
      $sample->{labels}->[0]->{label_color},
    )."\n";
  }
  close FILE;

}

################################################################################

sub write_frequency_matrix {

=pod

Expects:


Returns:


=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No specification for output $file $!" }
  if (! $pgx->{frequencymaps}->[0]) { return $pgx }

  my $i         =   0;
  open  FILE, '>'."$file";
  print FILE join("\t",
    'label',
    ( map{ 'dup_'.($_+1) } @{ $pgx->{matrixindex} } ),
    ( map{ 'del_'.($_+1) } @{ $pgx->{matrixindex} } ),
  )."\n";
  foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {
    $i++;
    if ($frequencymapsSet->{name} =~ /\w\w/) {
      my $label =   $frequencymapsSet->{name};
      $label    =~  s/[^\w\,]/_/g;
      $label    =~  s/_$//g;
      print FILE $label."\t" }
    else {
      print FILE 'subset_'.$i."\t" }
    print FILE join("\t",
      (map{ $frequencymapsSet->{dupfrequencies}->[$_] + 0 } @{ $pgx->{matrixindex} }),
      (map{ $frequencymapsSet->{delfrequencies}->[$_] + 0 } @{ $pgx->{matrixindex} }),
    )."\n";
  }
  close FILE;

}

################################################################################

sub write_status_matrix {

  no warnings 'uninitialized';

=pod

Expects:


Returns:


=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No specification for output $file $!" }
  if (! $pgx->{samples}->[0]->{statusmaps}) { return $pgx }

  my $i         =   0;
  open  FILE, '>'."$file";
  print FILE join("\t",
    'label',
    ( map{ 'dup_'.($_+1) } @{ $pgx->{matrixindex} } ),
    ( map{ 'del_'.($_+1) } @{ $pgx->{matrixindex} } ),
  )."\n";
  foreach my $sample (@{ $pgx->{samples} }) {
    $i++;
    my $label   =   'sample_'.$i;
    if ($sample->{id} =~ /\w\w/) {
      $label    =   $sample->{id} }
    elsif ($sample->{UID} =~ /\w\w/) {
      $label    =   $sample->{UID} }
    $label      =~  s/[^\w\,\-]/_/g;
    $label      =~  s/_$//g;      
    print FILE $label."\t";
    print FILE join("\t",
      (map{ $sample->{statusmaps}->{dupmap}->[$_] eq 'DUP' ? 1 : 0 } @{ $pgx->{matrixindex} }),
      (map{ $sample->{statusmaps}->{delmap}->[$_] eq 'DEL' ? 1 : 0 } @{ $pgx->{matrixindex} })
    )."\n";
  }
  close FILE;

}

################################################################################

sub write_value_matrix {

  no warnings 'uninitialized';

=pod

Expects:


Returns:


=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $file      =   shift;

  if (! $file) { die "No file $file $!" }
  if (! $pgx->{samples}->[0]->{statusmaps}) { return $pgx }

  my $i         =   0;
  open  FILE, '>'."$file";
  print FILE join("\t",
    'label',
    ( map{ 'dup_'.($_+1) } @{ $pgx->{matrixindex} } ),
    ( map{ 'del_'.($_+1) } @{ $pgx->{matrixindex} } ),
  )."\n";
  foreach my $sample (@{ $pgx->{samples} }) {
    $i++;
    if ($sample->{id} =~ /\w\w/) {
      my $label =   $sample->{id};
      $label    =~  s/[^\w\,]/_/g;
      $label    =~  s/_$//g;
      print FILE $label."\t" }
    else {
      print FILE 'sample_'.$i."\t" }
    print FILE join("\t",
      (map{ $sample->{statusmaps}->{dupmax}->[$_] + 0 } @{ $pgx->{matrixindex} }),
      (map{ $sample->{statusmaps}->{delmin}->[$_] + 0 } @{ $pgx->{matrixindex} }),
    )."\n";
  }
  close FILE;

}


1;

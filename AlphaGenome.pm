=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 AlphaGenome

=head1 SYNOPSIS

 mv AlphaGenome.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin AlphaGenome,file=/path/to/alphageome_scores.tsv.gz

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves pre-computed variant effect predictions
 from AlphaGenome (Google DeepMind). AlphaGenome is a deep learning model that
 inputs 1 Mb of DNA sequence and predicts functional genomic tracks at
 single-base-pair resolution across 11 modalities:

   RNA_SEQ            - RNA-seq gene expression
   CAGE               - CAGE expression at transcription start sites
   PROCAP             - PRO-cap transcription initiation
   DNASE              - DNase I chromatin accessibility
   ATAC               - ATAC-seq chromatin accessibility
   CHIP_HISTONE       - ChIP-seq histone modifications
   CHIP_TF            - ChIP-seq transcription factor binding
   SPLICE_SITES       - Splice site donor/acceptor probabilities
   SPLICE_SITE_USAGE  - Splice site usage fractions
   SPLICE_JUNCTIONS   - Splice junction counts (donor-acceptor pairs)
   CONTACT_MAPS       - 3D chromatin interaction maps

 The plugin reads a tabix-indexed TSV file containing scores exported from the
 AlphaGenome Python SDK via the tidy_scores() function.

 By default all modalities are reported. Use the 'modalities' parameter to
 restrict output to specific modalities of interest.

 Please cite the AlphaGenome publication alongside Ensembl VEP if you use
 this resource:
 https://www.nature.com/articles/s41586-025-10014-0

 AlphaGenome SDK: https://github.com/google-deepmind/alphagenome

 Running options:

  file        : (required) Path to tabix-indexed AlphaGenome scores TSV

  cutoff      : Minimum absolute raw score to report (default: 0, report all)

  quantile_cutoff : Minimum absolute quantile score to report (default: 0)

  modalities  : Colon-separated list of modalities to include
                (default: all modalities)
                Example: modalities=SPLICE_JUNCTIONS:SPLICE_SITES:RNA_SEQ

 Output:

  For tab/VCF output the result is reported as a pipe-delimited string:
    GENE_NAME|OUTPUT_TYPE|TRACK_NAME|RAW_SCORE|QUANTILE_SCORE
  When multiple modalities/tracks pass filters, the highest absolute raw
  score result is reported.

  For JSON output, all passing results are returned as:
    "AlphaGenome": [
      {"gene_id": "...", "gene_name": "...", "output_type": "...",
       "track_name": "...", "raw_score": ..., "quantile_score": ...},
      ...
    ]

 Data preparation:

  1. Run AlphaGenome batch variant scoring using the Python SDK
  2. Export scores via tidy_scores() to a pandas DataFrame
  3. Write to TSV with columns:
       #CHROM  POS  REF  ALT  GENE_ID  GENE_NAME  OUTPUT_TYPE  TRACK_NAME  RAW_SCORE  QUANTILE_SCORE
  4. Sort, compress and index:
       sort -k1,1 -k2,2n alphageome_scores.tsv | bgzip -c > alphageome_scores.tsv.gz
       tabix -s 1 -b 2 -e 2 alphageome_scores.tsv.gz

 The tabix utility must be installed in your path to use this plugin.
 Check https://github.com/samtools/htslib.git for instructions.

=cut

package AlphaGenome;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my %VALID_MODALITIES = map { $_ => 1 } qw(
  RNA_SEQ CAGE PROCAP DNASE ATAC CHIP_HISTONE CHIP_TF
  SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS CONTACT_MAPS
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();

  # File parameter is required
  my $file = $params->{file};
  die "ERROR: file parameter is required for AlphaGenome plugin, e.g.:\n" .
    "  --plugin AlphaGenome,file=/path/to/alphageome_scores.tsv.gz\n"
    unless $file;

  $self->add_file($file);

  # Cutoff for raw scores (default: 0 = no filtering)
  $self->{cutoff} = defined($params->{cutoff}) ? $params->{cutoff} : 0;

  # Cutoff for quantile scores (default: 0 = no filtering)
  $self->{quantile_cutoff} = defined($params->{quantile_cutoff}) ? $params->{quantile_cutoff} : 0;

  # Parse modality filter
  if (defined($params->{modalities})) {
    my @mods = split(/:/, $params->{modalities});
    my %mod_filter;
    for my $m (@mods) {
      die "ERROR: Unknown modality '$m'. Valid modalities are: " .
        join(', ', sort keys %VALID_MODALITIES) . "\n"
        unless $VALID_MODALITIES{$m};
      $mod_filter{$m} = 1;
    }
    $self->{modality_filter} = \%mod_filter;
  }

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  return {
    AlphaGenome => 'AlphaGenome variant effect predictions. ' .
      'Format: GENE_NAME|OUTPUT_TYPE|TRACK_NAME|RAW_SCORE|QUANTILE_SCORE. ' .
      'See https://www.nature.com/articles/s41586-025-10014-0'
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;

  my $alt_alleles = $tva->base_variation_feature->alt_alleles;
  my $ref_allele = $vf->ref_allele_string;

  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if ($vf_start > $vf_end);

  my @data = @{
    $self->get_data(
      $vf->{chr},
      $vf_start,
      $vf_end
    )
  };

  return {} unless @data;

  my @passing_results;

  foreach my $row (@data) {
    # Match alleles
    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref_allele,
        alts   => $alt_alleles,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $row->{ref},
        alts => [$row->{alt}],
        pos  => $row->{start},
      }
    );
    next unless @$matches;

    my $result = $row->{result};

    # Apply modality filter
    if ($self->{modality_filter}) {
      next unless $self->{modality_filter}{$result->{output_type}};
    }

    # Apply raw score cutoff
    if ($self->{cutoff} > 0) {
      next unless defined($result->{raw_score}) &&
        $result->{raw_score} ne '.' &&
        abs($result->{raw_score}) >= $self->{cutoff};
    }

    # Apply quantile score cutoff
    if ($self->{quantile_cutoff} > 0) {
      next unless defined($result->{quantile_score}) &&
        $result->{quantile_score} ne '.' &&
        abs($result->{quantile_score}) >= $self->{quantile_cutoff};
    }

    push @passing_results, $result;
  }

  return {} unless @passing_results;

  # JSON output: return all passing results
  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    return {
      AlphaGenome => \@passing_results
    };
  }

  # Tab/VCF output: return the top-scoring result as pipe-delimited string
  my $top = $passing_results[0];
  my $top_abs = abs($top->{raw_score} || 0);
  for my $r (@passing_results[1..$#passing_results]) {
    my $abs = abs($r->{raw_score} || 0);
    if ($abs > $top_abs) {
      $top = $r;
      $top_abs = $abs;
    }
  }

  return {
    AlphaGenome => join('|',
      $top->{gene_name}      // '.',
      $top->{output_type}    // '.',
      $top->{track_name}     // '.',
      $top->{raw_score}      // '.',
      $top->{quantile_score} // '.',
    )
  };
}

sub parse_data {
  my ($self, $line) = @_;
  chomp $line;
  my @f = split /\t/, $line;

  return {
    chr   => $f[0],
    start => $f[1],
    ref   => $f[2],
    alt   => $f[3],
    result => {
      gene_id        => $f[4],
      gene_name      => $f[5],
      output_type    => $f[6],
      track_name     => $f[7],
      raw_score      => $f[8],
      quantile_score => $f[9],
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

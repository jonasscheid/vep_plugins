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
 ./vep -i variations.vcf --plugin AlphaGenome,file=/path/to/alphagenome_scores.tsv.gz

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

 The plugin reads a tabix-indexed TSV file in wide format with one row per
 variant and pre-aggregated top scores per modality. The data preparation step
 aggregates the full AlphaGenome tidy_scores() output (thousands of rows per
 variant across cell types and tracks) into a single row per variant.

 Please cite the AlphaGenome publication alongside Ensembl VEP if you use
 this resource:
 https://www.nature.com/articles/s41586-025-10014-0

 AlphaGenome SDK: https://github.com/google-deepmind/alphagenome

 Running options:

  file        : (required) Path to tabix-indexed AlphaGenome scores TSV

  cutoff      : Minimum absolute raw score to report a modality
                (default: 0, report all)

  quantile_cutoff : Minimum absolute quantile score to report a modality
                    (default: 0)

  modalities  : Plus-separated list of modalities to include
                (default: all modalities)
                Example: modalities=SPLICE_JUNCTIONS+SPLICE_SITES+RNA_SEQ

 Output:

  The plugin reports one VEP field per modality:
    AlphaGenome_GENE_ID       - Ensembl gene ID (top-scoring gene)
    AlphaGenome_GENE_NAME     - Gene symbol (top-scoring gene)
    AlphaGenome_ATAC          - ATAC raw_score/quantile_score
    AlphaGenome_CAGE          - CAGE raw_score/quantile_score
    AlphaGenome_CHIP_HISTONE  - ChIP histone raw_score/quantile_score
    AlphaGenome_CHIP_TF       - ChIP TF raw_score/quantile_score
    AlphaGenome_CONTACT_MAPS  - Contact maps raw_score/quantile_score
    AlphaGenome_DNASE         - DNase raw_score/quantile_score
    AlphaGenome_PROCAP        - PRO-cap raw_score/quantile_score
    AlphaGenome_RNA_SEQ       - RNA-seq raw_score/quantile_score
    AlphaGenome_SPLICE_JUNCTIONS   - Splice junctions raw_score/quantile_score
    AlphaGenome_SPLICE_SITES       - Splice sites raw_score/quantile_score
    AlphaGenome_SPLICE_SITE_USAGE  - Splice site usage raw_score/quantile_score

  Score format: RAW_SCORE/QUANTILE_SCORE (e.g. -0.0234/0.876)

  For JSON output, a structured hash is returned with all fields.

 Data preparation:

  1. Run AlphaGenome batch variant scoring using the Python SDK
  2. Export scores via tidy_scores() to a pandas DataFrame
  3. Aggregate to wide format (one row per variant) with columns:
       #CHROM  POS  REF  ALT  GENE_ID  GENE_NAME  ATAC  CAGE
       CHIP_HISTONE  CHIP_TF  CONTACT_MAPS  DNASE  PROCAP  RNA_SEQ
       SPLICE_JUNCTIONS  SPLICE_SITES  SPLICE_SITE_USAGE
     Each modality column contains: RAW_SCORE/QUANTILE_SCORE
     GENE_ID/GENE_NAME: gene with highest absolute score across
     gene-centric modalities, or "." if none.
     Use "." for modalities without scores.
  4. Sort, compress and index:
       sort -k1,1 -k2,2n alphagenome_scores.tsv | bgzip -c > alphagenome_scores.tsv.gz
       tabix -s 1 -b 2 -e 2 alphagenome_scores.tsv.gz

 The tabix utility must be installed in your path to use this plugin.
 Check https://github.com/samtools/htslib.git for instructions.

=cut

package AlphaGenome;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my @MODALITIES = qw(
  ATAC CAGE CHIP_HISTONE CHIP_TF CONTACT_MAPS DNASE PROCAP
  RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE
);

my %VALID_MODALITIES = map { $_ => 1 } @MODALITIES;

# Column indices in the wide-format TSV
my %COL_IDX;
my $idx = 0;
for my $col (qw(CHROM POS REF ALT GENE_ID GENE_NAME), @MODALITIES) {
  $COL_IDX{$col} = $idx++;
}

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
    "  --plugin AlphaGenome,file=/path/to/alphagenome_scores.tsv.gz\n"
    unless $file;

  $self->add_file($file);

  # Cutoff for raw scores (default: 0 = no filtering)
  $self->{cutoff} = defined($params->{cutoff}) ? $params->{cutoff} : 0;

  # Cutoff for quantile scores (default: 0 = no filtering)
  $self->{quantile_cutoff} = defined($params->{quantile_cutoff}) ? $params->{quantile_cutoff} : 0;

  # Parse modality filter
  if (defined($params->{modalities})) {
    my @mods = split(/\+/, $params->{modalities});
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

  my $cite = 'See https://www.nature.com/articles/s41586-025-10014-0';
  my %header = (
    AlphaGenome_GENE_ID   => "AlphaGenome top-scoring gene Ensembl ID. $cite",
    AlphaGenome_GENE_NAME => "AlphaGenome top-scoring gene symbol. $cite",
  );

  for my $mod (@MODALITIES) {
    $header{"AlphaGenome_$mod"} =
      "AlphaGenome $mod variant effect score (raw_score/quantile_score). $cite";
  }

  return \%header;
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

    # Wide format: one row per variant, so first match is the result
    my $result = $row->{result};
    my %output;

    # Gene info
    my $gene_id   = $result->{gene_id};
    my $gene_name = $result->{gene_name};
    $output{AlphaGenome_GENE_ID}   = $gene_id   if defined($gene_id)   && $gene_id   ne '.';
    $output{AlphaGenome_GENE_NAME} = $gene_name  if defined($gene_name) && $gene_name ne '.';

    # Modality scores
    for my $mod (@MODALITIES) {
      # Apply modality filter
      if ($self->{modality_filter}) {
        next unless $self->{modality_filter}{$mod};
      }

      my $val = $result->{$mod};
      next unless defined($val) && $val ne '.';

      # Parse raw_score/quantile_score
      my ($raw, $quantile) = split(/\//, $val, 2);

      # Apply raw score cutoff
      if ($self->{cutoff} > 0) {
        next unless defined($raw) && $raw ne '' && abs($raw) >= $self->{cutoff};
      }

      # Apply quantile score cutoff
      if ($self->{quantile_cutoff} > 0) {
        next unless defined($quantile) && $quantile ne '' &&
          abs($quantile) >= $self->{quantile_cutoff};
      }

      $output{"AlphaGenome_$mod"} = $val;
    }

    return \%output if %output;
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  chomp $line;
  my @f = split /\t/, $line;

  my %result = (
    gene_id   => $f[$COL_IDX{GENE_ID}],
    gene_name => $f[$COL_IDX{GENE_NAME}],
  );

  for my $mod (@MODALITIES) {
    $result{$mod} = $f[$COL_IDX{$mod}];
  }

  return {
    chr    => $f[$COL_IDX{CHROM}],
    start  => $f[$COL_IDX{POS}],
    ref    => $f[$COL_IDX{REF}],
    alt    => $f[$COL_IDX{ALT}],
    result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

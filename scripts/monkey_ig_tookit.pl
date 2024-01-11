#!/usr/bin/env perl
package MyApp {
    use MooseX::App qw(Color);
    use Log::Any '$log';

    has 'log' => (
        is            => 'ro',
        isa           => 'Log::Any::Proxy',
        required      => 1,
        default       => sub { Log::Any->get_logger },
        documentation => 'Keep Log::Any::App object',
    );

    __PACKAGE__->meta->make_immutable;
}

package MyApp::FixHeavyChain {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use Bio::SeqIO;
    use namespace::autoclean;

    command_short_description q[Correct FASTA insertion sequence containg heavy chaind (it looks for _HC_ and FWD_T4 in sequence name)];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[Input fasta sequence!],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 1,
        documentation => q[Output fasta sequence!],
    );


    sub run {
        my ($self) = @_;
        my $in = Bio::SeqIO->new(
            -file   => $self->input_file->stringify,
            -format => 'fasta'
        );
        my $out = Bio::SeqIO->new(
            -file   => '>' . $self->output_file->stringify,
            -format => 'fasta'
        );
        while ( my $seq = $in->next_seq ) {
            if ( $seq->id =~ /\_HC\_/i && $seq->id =~ /FWD\_T4/i ) {

                # search pattern GGTG-AGC
                my $dna = $seq->seq;
                if ( $dna =~ /GGTGAGC/i ) {
                    if ( $-[0] <= 80 ) {
                        $dna =~ s/GGTGAGC/GGTGAAGC/g;
                        $seq->seq($dna);
                        #$seq->id($seq->id.'_fixed');
                    }
                }
                # Second fix (VH5.7)
                elsif ( $dna =~ /AGCAGAGGTGAAAGGCCC/i ) {
                    if ( $-[0] <= 80 ) {
                        $dna =~ s/AGCAGAGGTGAAAGGCCC/AGCAGAGGTGAAAAGGCCC/g;
                        $seq->seq($dna);
                        #$seq->id($seq->id.'_fixed2');
                    }
                }
                # 3rd fix (VH3.15)
                elsif ( $dna =~ /GCTTGGCAAGCCTG/ ){
                    if ( $-[0] <= 80 ) {
                        $dna =~ s/GCTTGGCAAGCCTG/GCTTGGCAAAGCCTG/g;
                        $seq->seq($dna);
                        #$seq->id($seq->id.'_fixed3');
                    }
                }
            }
            $out->write_seq($seq);
        }
    }

    __PACKAGE__->meta->make_immutable;
}

package MyApp::PairSequenceChains {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use Bio::SeqIO;
    use Data::Printer;
    use namespace::autoclean;

    command_short_description q[Correct FASTA insertion sequence containg heavy chaind (it looks for _HC_ and FWD_T4 in sequence name)];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[Input fasta sequence!],
    );

    has_file 'index_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(d)],
        required      => 1,
        documentation => q[Index File],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 0,
        documentation => q[Output fasta sequence!],
    );


    sub parse_index_file {
        my ($self) = @_;
        open( my $in, '<', $self->index_file ) 
            || die "Cannot open/read file " . $self->index_file . "!";
        my %h;
        while ( my $row = <$in> ){
            chomp $row;
            my @F = split "\t",$row;
            push @{$h{$F[0].'_'.$F[1]}},$F[2];
        }
        close( $in );
        return \%h;
    }

    sub parse_fasta_file {
        my ($self) = @_;
        my %h;
        my $in = Bio::SeqIO->new(
            -file   => $self->input_file->stringify,
            -format => 'fasta'
        );
         while ( my $seq = $in->next_seq ) {
            my $id = $seq->id;
            #$id =~ s/\_fixed$//g;
            $h{$id} = $seq;
        }
        return \%h;
    }
    
    sub get_codon_seq {
        my ($self, $string) = @_;
        my @aux = $string =~ /\w{3}/g;
        my $codon_string;
        foreach my $codon (@aux) {
            $codon_string .= $codon if length($codon) == 3;
        }
        return $codon_string;
    }

    sub run {
        my ($self) = @_;
        my $index = $self->parse_index_file;
        my $seqs = $self->parse_fasta_file;
        my @pairs;
        my @single_hc;
        my @single_lc;

        PAIR_KEYS:
        foreach my $k (sort {$a cmp $b} keys %{$index}) {
            #if (scalar @{$index->{$k}} == 2){
                my @aux = @{$index->{$k}};
                my $c=0;
                my $HC;
                my $LC;
                foreach my $i (@aux) {
                   if ($seqs->{$i}){
                       $c++;
                       $HC = $i if $i =~ /\_HC\_|FWD\_T4/i;
                       $LC = $i if $i =~ /\_LCL\_|LC\_LAMBDA/i;
                   }
                   else{
                       #warn "Cannot find sequence for: $i";
                       next PAIR_KEYS;
                   }
                
                }
                if ($c == 2 && $HC && $LC){
                    say "$k:";
                    say "\t".$seqs->{$HC}->id;
                    #say "\t".$aux[0];
                    say "\t".$seqs->{$LC}->id;
                    #say "\t".$aux[1]
                    my $id     = $seqs->{$HC}->id . '-' . $seqs->{$LC}->id;
                   
                    # Get sequences codon corrected (multiple of 3)
                    my $HC_seq = $self->get_codon_seq($seqs->{$HC}->seq);
                    my $LC_seq = $self->get_codon_seq($seqs->{$LC}->seq);


                    my $seqobj = Bio::Seq->new(
                        -display_id => $id,
                        -seq        => $HC_seq . $LC_seq
                    );
                    push @pairs, $seqobj;
                    delete $seqs->{$HC};
                    delete $seqs->{$LC};
                }
                elsif ($HC){
                    push @single_hc, $seqs->{$HC};
                    delete $seqs->{$HC};
                }
                elsif ($LC){
                    push @single_lc, $seqs->{$LC};
                    delete $seqs->{$LC};
                }
                else{
                    die "Cannot define chain type: ". join ", ", @aux;
                }
            #}
        }

        # Get sequences without pairing information if any:
        foreach my $i (sort {$a cmp $b} keys %{$seqs}) {
            my ($HC,$LC);
            $HC = $i if $i =~ /\_HC\_|FWD\_T4/i;
            $LC = $i if $i =~ /\_LCL\_|LC\_LAMBDA/i;
            if ($HC) {
                push @single_hc, $seqs->{$HC};
                delete $seqs->{$HC};
            }
            elsif ($LC) {
                push @single_lc, $seqs->{$LC};
                delete $seqs->{$LC};
            }
            else{
                die "Cannot define chain type: $i";
            }
        }

        my $out_pairs = Bio::SeqIO->new(
            -file   => '>' . 'pairs.fasta',
            -format => 'fasta'
        );

        foreach my $seq (@pairs) {
            $out_pairs->write_seq($seq);
        }

        my $out_hc = Bio::SeqIO->new(
            -file   => '>' . 'single_hc.fasta',
            -format => 'fasta'
        );

        foreach my $seq (@single_hc) {
            $out_hc->write_seq($seq);
        }

        my $out_lc = Bio::SeqIO->new(
            -file   => '>' . 'single_lc.fasta',
            -format => 'fasta'
        );

        foreach my $seq (@single_lc) {
            $out_lc->write_seq($seq);
        }

        say join "\n", keys %{$seqs};

    }

    __PACKAGE__->meta->make_immutable;
}


package MyApp::AirrToGenbank {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use Bio::SeqIO;
    use Bio::Seq;
    use Bio::SeqFeature::Generic;
    use Data::Printer;
    use namespace::autoclean;

    command_short_description q[Parses a Airr igblast output and create genbank annotated sequences.];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[Input fasta sequence!],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 0,
        documentation => q[Output fasta sequence!],
    );


    sub parse_airr_file {
        my ($self) = @_;
        open( my $in, '<', $self->input_file ) 
            || die "Cannot open/read file " . $self->input_file . "!";
        my @sequences;
        my @header;
        while ( my $row = <$in> ){
            chomp $row;
            my @F = split "\t", $row;
            if ($F[0] eq 'sequence_id' && ! @header){
                @header = @F;
            }
            else{
                my %h;
                @h{@header} = @F;
                push @sequences, \%h;
            }
            
        }
        close( $in );
        return \@sequences;
    }

    sub run {
        my ($self) = @_;
        my $sequences = $self->parse_airr_file;

        my $out = Bio::SeqIO->new( -file => ">test.gb", -format => 'genbank' );

        # Annotate all columns with start and end
        my @regions;
        foreach my $f (sort {$a cmp $b } keys %{$sequences->[0]} ) {
            next if $f =~ /germline|alignment/;
            push @regions, $1 if $f =~ /(\S+)\_end/;
        }

        foreach my $seq ( @{$sequences} ) {
            my $seqobj = Bio::Seq->new(
                -id  => $seq->{sequence_id},
                -seq => $seq->{sequence}
            );
            #say $seq->{sequence_id};
            my @features;
            foreach my $r (@regions) {
                next unless defined($seq->{ $r . '_start'  }) && defined($seq->{ $r . '_end'  });
                next unless $seq->{ $r . '_start'  } =~ /\d+/ && $seq->{ $r . '_end'  } =~ /\d+/;
                #say $r . "\t". $seq->{ $r . '_start'  }.".." . $seq->{ $r . '_end'  };
                my $feat = new Bio::SeqFeature::Generic(
                    -start       => $seq->{ $r . '_start' } + 1,
                    -end         => $seq->{ $r . '_end' },
                    -strand      => 1,
                    -primary_tag => uc($r),
                    -tag         => {
                        evidence => 'predicted',
                    }
                );
                push @features, $feat;
            }
            
            if ( $seq->{'productive'} ) {

                # Define as productive or not
                my $status = 'PRODUCTIVE';
                $status = 'NOT_PRODUCTIVE' if $seq->{'productive'} eq 'F';

                my $feat = new Bio::SeqFeature::Generic(
                    -start       => 1,
                    -end         => length $seq->{'sequence'},
                    -strand      => 1,
                    -primary_tag => uc($status),
                    -tag         => {
                        evidence => 'predicted',
                    }
                );

                push @features, $feat;
            }
            if ( $seq->{'locus'} ) {

                my $feat = new Bio::SeqFeature::Generic(
                    -start       => 1,
                    -end         => length $seq->{'sequence'},
                    -strand      => 1,
                    -primary_tag => uc($seq->{'locus'}),
                    -tag         => {
                        evidence => 'predicted',
                    }
                );

                push @features, $feat;
            }
          
            $seqobj->add_SeqFeature(@features);
            $out->write_seq($seqobj);
        }
    }

    __PACKAGE__->meta->make_immutable;
}

package MyApp::CreatePairedClones {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use Bio::SeqIO;
    use Bio::Seq;
    use Bio::SeqFeature::Generic;
    use Data::Printer;
    use namespace::autoclean;

    command_short_description q[Parses ChangeO output for heavy and light chainsand create clones.];

    has_file 'input_HC_clones' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[ChangeO clone pass file for Heavy Chain sequences.],
    );

    has_file 'input_LC_clones' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(l)],
        required      => 1,
        documentation => q[ChangeO clone pass file for Light Chain sequences.],
    );

    has_file 'index_pairs' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(p)],
        required      => 1,
        documentation => q[TSV file with sequence pairs.],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 0,
        documentation => q[Output fasta sequence!],
    );


    sub parse_changeo_clone_file {
        my ( $self, $file ) = @_;
        open( my $in, '<', $file )
          || die "Cannot open/read file " . $file . "!";
        my %sequences;
        my @header;
        while ( my $row = <$in> ) {
            chomp $row;
            $row =~ s/\s+$//g;
            $row =~ s/^\s+$//g;
            my @F = split "\t", $row;
            if ( $F[0] =~ /sequence_id/i && !@header ) {
                @header = @F;
            }
            else {
                my %h;
                @h{@header} = @F;
                my $id = $h{SEQUENCE_ID};

                # remove fixed suffix
                #$id =~ s/_fixed.*$//g;
                $sequences{$id} = \%h;
            }
        }
        close($in);
        return \%sequences;
    }

    sub parse_index_pairs {
        my ($self, $file) = @_;;
        my (%h,%l);
        open( my $in, '<', $file )
            || die "Cannot open/read file " . $file . "!";

        while ( my $row = <$in> ){
            chomp $row;
            next if $row =~ /^$/;
            my @F = split "\t", $row;
            die "$file doesnt'contain two columns" unless scalar @F == 2;
            # remove sequence extension
            $F[0] =~ s/\.ab1//ig;
            $F[1] =~ s/\.ab1//ig;
            # remove fixed suffix
            $F[0] =~ s/_fixed.*$//g;
            $F[1] =~ s/_fixed.*$//g;
            
            $h{$F[0]} = $F[1];
            $l{$F[1]} = $F[0];
        }
        close( $in );
        return (\%h,\%l);
    }

    sub check_index {
        my ($self, $index, $seqs);
        foreach my $id (sort {$a cmp $b} keys %{$seqs}) {
            die "Cannot find $id" unless $index->{$id};
        }
    }

    sub run {
        my ($self) = @_;
        my $hc = $self->parse_changeo_clone_file($self->input_HC_clones->stringify);
        my $lc = $self->parse_changeo_clone_file($self->input_LC_clones->stringify);
        my ($h_index, $l_index ) = $self->parse_index_pairs($self->index_pairs->stringify);

        $self->check_index($h_index, $hc);
        $self->check_index($l_index, $lc);

        # use heavy chain to go over pairs:
        my %clones;
        foreach my $h_seq (sort {$a cmp $b} keys %{$h_index}) {
            my $l_seq = $h_index->{$h_seq};
            if ($hc->{$h_seq} && $lc->{$l_seq}){
                my $h_id = $hc->{$h_seq}->{CLONE};
                my $h_func = $hc->{$h_seq}->{FUNCTIONAL};
                my $h_v = $hc->{$h_seq}->{V_CALL};
                my $h_d = $hc->{$h_seq}->{D_CALL};
                my $h_j = $hc->{$h_seq}->{J_CALL};
                my $h_cdr3 = $hc->{$h_seq}->{CDR3_IMGT};
                my $h_bioseq = Bio::Seq->new(-seq => $h_cdr3);
                $h_cdr3 = $h_bioseq->translate->seq;
                my $l_id = $lc->{$l_seq}->{CLONE};
                my $l_func = $lc->{$l_seq}->{FUNCTIONAL};
                my $l_v = $lc->{$l_seq}->{V_CALL};
                my $l_j = $lc->{$l_seq}->{J_CALL};
                my $l_cdr3 = $lc->{$l_seq}->{CDR3_IMGT};
                my $l_bioseq = Bio::Seq->new(-seq => $l_cdr3);
                $l_cdr3 = $l_bioseq->translate->seq;
 
                say join "\t", ("$h_id,$l_id", $h_seq, $h_func, $h_v, $h_d, $h_j,$h_cdr3, $l_seq, $l_func, $l_v, $l_j, $l_cdr3)  ;
            }
            elsif ($hc->{$h_seq}){
                my $h_id = $hc->{$h_seq}->{CLONE};
                my $h_func = $hc->{$h_seq}->{FUNCTIONAL};
                my $h_v = $hc->{$h_seq}->{V_CALL};
                my $h_d = $hc->{$h_seq}->{D_CALL};
                my $h_j = $hc->{$h_seq}->{J_CALL};
                my $h_cdr3 = $hc->{$h_seq}->{CDR3_IMGT};
                my $h_bioseq = Bio::Seq->new(-seq => $h_cdr3);
                $h_cdr3 = $h_bioseq->translate->seq;
                say join "\t", ("$h_id,", $h_seq, $h_func, $h_v, $h_d, $h_j,$h_cdr3, $l_seq, '', '', '', '')  ;
            }
            elsif ($lc->{$l_seq}){
                my $l_id = $lc->{$l_seq}->{CLONE};
                my $l_func = $lc->{$l_seq}->{FUNCTIONAL};
                my $l_v = $lc->{$l_seq}->{V_CALL};
                my $l_j = $lc->{$l_seq}->{J_CALL};
                my $l_cdr3 = $lc->{$l_seq}->{CDR3_IMGT};
                my $l_bioseq = Bio::Seq->new(-seq => $l_cdr3);
                $l_cdr3 = $l_bioseq->translate->seq;
                say join "\t", (",$l_id", $h_seq, '', '', '', '','', $l_seq, $l_func, $l_v, $l_j, $l_cdr3)  ;
            }
            else {
                say join "\t", (",", $h_seq, '', '', '', '','', $l_seq, '', '', '', '')  ;
            }
        }
   }

    __PACKAGE__->meta->make_immutable;
}


use MyApp;
#use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();


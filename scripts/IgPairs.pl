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

package MyApp::ParseExcel {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use Spreadsheet::ParseXLSX;
    use Spreadsheet::Read;
    use Excel::Writer::XLSX;
    use Spreadsheet::WriteExcel::Styler;
    use Data::Printer;
    use Bio::Seq;
    use namespace::autoclean;

    command_short_description q[Parse XLSX from IgInterface];
    command_long_description q[Parse XLSX from IgInterface];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        must_exist    => 1,
        documentation => q[Very important option!],
    );

    has_file 'pairs_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(p)],
        required      => 1,
        must_exist      => 1,
        documentation => q[File containing pair information!],
    );

    option 'output_prefix' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases   => [qw(o)],
        required      => '1',
        documentation => q[Output prefix.],
    );


    our @TABLE_HEADER;

    sub parse_clone_summary {
        my ($self) = @_;
        my %structure;
        my %indexed_sequences;

        my $parser = Spreadsheet::ParseXLSX->new();
        my $workbook = $parser->parse( $self->input_file->stringify );

        if ( !defined $workbook ) {
            die $parser->error(), ".\n";
        }
        for my $worksheet ( $workbook->worksheets() ) {
            my $worksheet_name = $worksheet->get_name;

            # Filter some  worksheets
            next if $worksheet_name =~ /chart/i;
            next unless $worksheet_name =~ /filtered/i;

            my $chain = $worksheet_name;
            $chain =~ s/^(\S+).*/$1/g;
            $chain = uc $chain;

            my ( $row_min, $row_max ) = $worksheet->row_range();
            my ( $col_min, $col_max ) = $worksheet->col_range();

            my @header;
            my $i;
            my ($clone_id, $n_seqs); # hold clone_id and n_seqs because those those are merged cells
            for my $row ( $row_min .. $row_max ) {
                my @values;
                for my $col ( $col_min .. $col_max ) {

                    my $cell = $worksheet->get_cell( $row, $col );
                    next unless $cell;

                    # Get header if it is the first row
                    unless ($i){
                        push @header, $cell->unformatted();
                    }
                    # or get values
                    else {
                        push @values, $cell->unformatted();
                    }
                }
                if (@values && $values[2]){
                    # Update clone_id variable using valeu from the first row
                    # with the clone_id
                    if ($values[0]){
                        $clone_id = $values[0];
                    }
                    # Update values with clone information if is empty
                    else {
                        die "Something is wrong! I cannot find clone id for row: $row" unless $clone_id;
                        $values[0] = $clone_id;
                    }

                    my %row_hash;
                    @row_hash{@header} = @values;

                    # add to the structure 
                    unless ($structure{$chain}{$row_hash{'SEQUENCE_ID'}}){
                        $structure{$chain}{$row_hash{'SEQUENCE_ID'}} = \%row_hash;
                    }
                    else{
                        die "Something is wrong. $row_hash{SEQUENCE_ID} is repeated in the same spreadsheet!";
                    }

                    # Add to index hash
                    $indexed_sequences{$row_hash{'SEQUENCE_ID'}} = {chain => $chain, clone_id => $clone_id };
                }
                $i++;
            }
            @TABLE_HEADER = @header;
        }
        return \%structure;
    }


    sub parse_pairs_files {
        my ($self) = @_;
        my %index;
        my @pairs;
        my $workbook = Spreadsheet::Read->new ($self->pairs_file->stringify);
        my $info = $workbook->[0];
        my $sheet = $workbook->sheet(1);

        # check if there is 3 columns
        if ($sheet->maxcol != 3){
            die $self->pairs_file . " doesn't look a valid pair files. It must have only 3 columns named: heavy, kappa, lambda. Please, check your file!";
        }
        my @header;
        foreach my $row (1.. $sheet->maxrow) {
            if ($row == 1){
                foreach my $value ($sheet->row($row)){
                    push @header, uc $value;
                }
                next;
            }
            my %hash;
            @hash{@header} = $sheet->row($row);
            push @pairs, \%hash;
            #create index;
            foreach my $seq_id ($sheet->row($row)) {
                next unless $seq_id;
                $index{$seq_id} =\%hash;
            }
        }
        return (\%index, \@pairs);
    }


    sub check_not_used_sequences_by_make_pairs {
        my ( $self, $seq_classification ) = @_;
        foreach my $chain ( sort { $a cmp $b } keys %{$seq_classification} ) {
            foreach my $seq_id ( sort { $a cmp $b }
                keys %{ $seq_classification->{$chain} } )
            {
                die "Something is wrong! Sequence $seq_id was not used for pairing. Please check sequence names or if all the sequence pairs are in the file '"
                  . $self->pairs_file . "'";
            }
        }
    }


    sub make_pairs {
        my ( $self, $seq_classification, $pairs ) = @_;
        my @pairs_classification;
        my $i = 0;
        foreach my $pair ( @{$pairs} ) {
            foreach my $chain ( sort { $a cmp $b } keys %{$pair} ) {
                my $seq_id = $pair->{$chain};
                next unless $seq_id;
                if ( $seq_classification->{$chain} ) {
                    if ( $seq_classification->{$chain}{$seq_id} ) {
                        $pairs_classification[$i]->{$chain} = $seq_classification->{$chain}{$seq_id};
                        delete $seq_classification->{$chain}{$seq_id};
                        delete $pairs->[$i]->{$chain};
                    }
                    else {
                        #push @not_found_seqs, $seq_id;
                    }
                }
                else {
                    die "Something is wrong! Cannot find chain: '$chain' !";
                }
            }
            $i++;
        }

        # Check if there is any sequence left
        $self->check_not_used_sequences_by_make_pairs($seq_classification);
        return \@pairs_classification;
    }


    sub generate_groups {
        my ($self, $pairs_classification) = @_;
        my %pair_types;
        my $i = 0;
        foreach my $pair (@{$pairs_classification}) {
            my @chains = keys %{$pair};
            if ($pair->{'HEAVY'} and scalar @chains == 2){
                foreach my $chain (@chains){
                    next if $chain =~ /heavy/i;
                    push @{$pair_types{'PROPER'}{$chain}}, $pair;
                    delete $pairs_classification->[$i];
                }
            }
            # tripplet
            elsif (scalar @chains == 3){
                    push @{$pair_types{'TRIPLETS'}}, $pair;
                    delete $pairs_classification->[$i];
            }
            # dual light chaing
            elsif (scalar @chains == 2){
                    push @{$pair_types{'LIGHT_ONLY'}}, $pair;
                    delete $pairs_classification->[$i];
            }
            # single chain
            elsif (scalar @chains == 1){
                    push @{$pair_types{'SINGLES'}{$chains[0]}}, $pair;
                    delete $pairs_classification->[$i];
            }
            $i++;
        }
        return \%pair_types;
    }


    sub find_clusters {
        my ($self, $pair_types) = @_;
        my %clusters;
        foreach my $type (sort {$a cmp $b} keys %{$pair_types}) {
            if ($type =~ /PROPER|SINGLES/){
               foreach my $class_chain (keys %{$pair_types->{$type}} ) {
                    my $info_list = $pair_types->{$type}{$class_chain};
                    foreach my $info (@{$info_list}) {
                        my @V_regions;
                        my @VJ_regions;
                        foreach my $chain (sort {$a cmp $b} keys %{$info}) {
                            push @V_regions, $info->{$chain}{V_CALL};
                            push @VJ_regions, join ",", @{$info->{$chain}}{('V_CALL','J_CALL')};
                        }
                        my $strict_cluster_key =  join "|", @VJ_regions;
                        my $relaxed_cluster_key =  join "|", @V_regions;
                        push @{ $clusters{'strict'}{$type}{$strict_cluster_key} }, $info;
                        push @{ $clusters{'relaxed'}{$type}{$relaxed_cluster_key} }, $info;
                    }
                }
            }
            else {
                my $info_list = $pair_types->{$type};
                foreach my $info ( @{$info_list} ) {
                    my @V_regions;
                    my @VJ_regions;
                    foreach my $chain ( sort { $a cmp $b } keys %{$info} ) {
                        push @V_regions,  $info->{$chain}{V_CALL};
                        push @VJ_regions, join ",",
                          @{ $info->{$chain} }{ ( 'V_CALL', 'J_CALL' ) };
                    }
                    my $strict_cluster_key  = join "|", @VJ_regions;
                    my $relaxed_cluster_key = join "|", @V_regions;
                    push @{ $clusters{'strict'}{$type}{$strict_cluster_key} },
                      $info;
                    push @{ $clusters{'relaxed'}{$type}{$relaxed_cluster_key} },
                      $info;
                }
            }
        }
        return \%clusters;
    }


    sub create_douhgnut_chart {
        my ($self, $workbook, $sheet_name, $data, $chart_colors, $styler ) = @_;
        my $worksheet = $workbook->add_worksheet($sheet_name);
        my $chart     = $workbook->add_chart( type => 'doughnut', embedded => 1 );
        my $total = 0;
        $total += $_ foreach @{$data->[1]}[ 1 .. $#{ $data->[1] }];
        $worksheet->write( 'A1', $data );

        # Configure the chart.
        $chart->add_series(
            categories => [ $sheet_name, 1, $#{ $data->[0] }, 0, 0 ],
            values     => [ $sheet_name, 1, $#{ $data->[1] }, 1, 1 ],
            points     => $chart_colors,
            border     => { color => 'black' },
        );
        my $title = $sheet_name;
        $title =~ s/(\S+?)\_.*/$1/g;
        $chart->set_title( name => $title);
        # Insert the chart into the worksheet.
        $worksheet->insert_chart( 'C2', $chart, {'x_offset' => 15, 'y_offset' => 7} );
        my $text = $workbook->add_shape(
            type   => 'rect',
            text   => $total,
            width  => 70,
            height => 60,
            line        => 'FFFFFF',
            format => $styler->('floating_text'),

        );
        $worksheet->insert_shape(
            'F9', $text
        );
    }


    sub create_xlsx {
        my ( $self, $clusters, $header, $xls_name ) = @_;

        my @header_first = qw/cluster_id num_seqs/;
        # Create a new Excel workbook
        my $workbook = Excel::Writer::XLSX->new($xls_name);

        #  Add and define a format
        # Create a styler with some styles
        my $styler = Spreadsheet::WriteExcel::Styler->new($workbook);
        my @colors = (
            qw [FF0000
              FF8000
              FFFF00
              80FF00
              00FF00
              00FF80
              00FFFF
              0080FF
              0000FF
              7F00FF
              FF00FF
              FF007F
              808080
              ]
        );
        my %color_hash = map { $_ => {bg_color => "#$_"} } @colors;
        $styler->add_styles(
            title => {
                align    => "center",
                bold     => 1,
                border => 1,
            },
            HEAVY => {
                bg_color => 'silver'
            },
            KAPPA => {
                bg_color => 'orange'
            },
            LAMBDA => {
                bg_color => 'yellow'
            },

            %color_hash,
            floating_text => {
                align    => "center",
                bold     => 1,
                size => 15
            },
           
        );
        my @types = qw/PROPER TRIPLETS LIGHT_ONLY SINGLES/;

        my $cluster_count = 1;
        foreach my $type (@types) {
            my @chart_colors;
            my @chart_cluster_id = ('cluster_id');
            my @chart_cluster_size = ('num_seqs');
            my $single_size = 0;
            my @local_colors = @colors;
            my $c_type = $clusters->{$type};
            my $i = 0;
            # Add a worksheet
            my $worksheet = $workbook->add_worksheet($type);
            # add header
            my @chains = qw/HEAVY KAPPA LAMBDA/;
            # add merged columns
            my $col_size = scalar @{$header};
            my $col_end = $col_size + 1;
            my $col_start = 2; # always start on 2
            foreach my $chain (@chains) {
                $worksheet->merge_range($i,$col_start,$i,$col_end,$chain, $styler->('title',$chain));
                $col_start += $col_size;
                $col_end += $col_size;
            }

            $i++;
            $worksheet->write_row($i,0,[ @header_first, @{$header}, @{$header}, @{$header}], $styler->(qw/title/));
            $i++;
            foreach my $cluster_key ( sort { $#{ $c_type->{$b} } <=> $#{ $c_type->{$a} } || $a cmp $b } keys %{$c_type} ) {
                my $cluster_size = scalar @{ $c_type->{$cluster_key} };

                my $color;
                if ( $cluster_size > 1 ) {
                    $color = shift @local_colors;
                    push @local_colors, $color;
                    push @chart_cluster_id, $cluster_count;
                    push @chart_cluster_size, $cluster_size;
                    push @chart_colors, { fill => { color => '#'.$color }, border => {color => '#000000'} }
                }
                else {
                    $single_size++;
                }

                foreach my $info ( @{ $c_type->{$cluster_key} } ) {
                    my @columns = ( $cluster_count, $cluster_size );
                    #foreach my $chain ( sort { $a cmp $b } keys %{$info} ) {
                    foreach my $chain (@chains) {
                        if ( $info->{$chain} ) {
                            push @columns, @{ $info->{$chain} }{ @{$header} };
                            # THIS IS A HACK TO TRANSLATE THE CORRECTED INPUT
                            # FROM START
                            my @translate = grep { $_ =~ /translation_of_/i } @{$header};
                            foreach my $column (@translate) {
                                my $from = $column;
                                $from =~ s/translation_of_//gi;
                                if ($info->{$chain}{$from}) {
                                    my $seq = Bio::Seq->new(-id => 'test', -seq => $info->{$chain}{$from});
                                    my @position = grep { $header->[$_] =~ /$column/ } 0 .. $#{$header};
                                    die "Error, ". scalar @position . " positions found in header!" if scalar @position != 1;
                                    $columns[($#columns - $#{$header}) + $position[0]] = $seq->translate->seq;
                                }
                                else{
                                    die "Cannot find $from required by $column";
                                }

                            }

                        }
                        else {
                            my @empty;
                            push @empty,   '' for 0 .. $#{$header};
                            push @columns, @empty;
                        }
                    }

                    #say "\t". join ",", @columns;
                    # Add color to cluster > 1
                    if ($color) {
                        $worksheet->write_row( $i, 0, \@columns,
                            $styler->($color) );
                    }
                    else {
                        $worksheet->write_row( $i, 0, \@columns );
                    }
                    $i++;
                }

                # skip line between clusters
                $worksheet->write_row( $i, 0, [] );
                $i++;
                $cluster_count++;
            }
            # Add chart
            if ($single_size) {
                push @chart_cluster_id,   'single';
                push @chart_cluster_size, $single_size;
                push @chart_colors, { fill => { color => '#FFFFFF' }, border => {color => '#000000'}  }
            }
            my @data = (\@chart_cluster_id,\@chart_cluster_size);
            my $chart_name = $type."_CHART";
            $self->create_douhgnut_chart($workbook,$chart_name,\@data, \@chart_colors, $styler);
        }
    }


    sub run {
        my ($self) = @_;
        my ( $index, $pairs ) = $self->parse_pairs_files;
        my $seq_classification = $self->parse_clone_summary;
        my $pairs_classification = $self->make_pairs( $seq_classification, $pairs );
        my $pair_types       = $self->generate_groups($pairs_classification);
        my $clusters = $self->find_clusters($pair_types);
        my @expected_header = @TABLE_HEADER[ 2 .. $#TABLE_HEADER ];
        $self->create_xlsx($clusters->{strict}, \@expected_header, $self->output_prefix . '_strict.xlsx');
        my @selected = qw/SEQUENCE_ID V_CALL D_CALL J_CALL corrected_input_from_start translation_of_corrected_input_from_start v_insertions v_deletions nt_mismatches_V_region v_region_aa_mismatches  cdr3_aa cdr3_aa_length/;
        $self->create_xlsx($clusters->{strict}, \@selected,  $self->output_prefix . '_strict_selected_columns.xlsx');
    }


    __PACKAGE__->meta->make_immutable;
}


use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

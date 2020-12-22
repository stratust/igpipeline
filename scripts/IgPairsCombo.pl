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
    use List::MoreUtils qw(uniq);
    use  Algorithm::Combinatorics qw(combinations);
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

    has 'chain_pairing_classification' => (
        is      => 'rw',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_pairing_classification',
        documentation => 'Hashref of chain pairing classification'
    );

    our @TABLE_HEADER;

    sub parse_clone_summary {
        my ($self, $index) = @_;
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
                    if ($index->{$row_hash{SEQUENCE_ID}}){
                        $row_hash{SAMPLEID} = $index->{$row_hash{SEQUENCE_ID}}->{SAMPLEID}
                    }
                    else {
                        p $index;
                        die "Cannot find SAMPLEID for $row_hash{SEQUENCE_ID}";
                    }

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
        if ($sheet->maxcol != 6){
            die $self->pairs_file . " doesn't look a valid pair files. It must have only 6 columns named: sampleid, igg, iga, igm, kappa, lambda. Please, check your file!";
        }
        my @header;
        foreach my $row (1.. $sheet->maxrow) {
            if ($row == 1){
                foreach my $value ($sheet->row($row)){
                    push @header, uc $value;
                }
                # Be sure that header is correct
                foreach my $value (@header) {
                    if ($value !~ /SAMPLEID|IGG|IGA|IGM|KAPPA|LAMBDA/){
                        die "Column name: $value is not acceptable! Valid column names are: sampleid, heavy, kappa, lambda";
                    }
                }
                next;
            }
            my %hash;
            @hash{@header} = $sheet->row($row);
            push @pairs, \%hash;
            #create index;
            foreach my $seq_id ($sheet->row($row)) {
                next unless $seq_id;
                next if $seq_id eq $hash{SAMPLEID};
                $index{$seq_id} = \%hash;
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
                next unless $chain =~ /IGG|IGM|IGA|LAMBDA|KAPPA/;
                my $igblast_chain = $chain;
                my $seq_id = $pair->{$chain};
                next unless $seq_id;
                #igchains
                if ($chain =~ /IGG|IGM|IGA/){
                    $igblast_chain = 'HEAVY';
                }
                if ( $seq_classification->{$igblast_chain} ) {
                    if ( $seq_classification->{$igblast_chain}{$seq_id} ) {
                        $pairs_classification[$i]->{$chain} = $seq_classification->{$igblast_chain}{$seq_id};
                        delete $seq_classification->{$igblast_chain}{$seq_id};
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

    sub _is_proper_paired {
        my ($self, $chains) = @_;
        my $result;
        my $chain_pair = join ",", sort @{$chains};
        foreach my $heavy ('IGA','IGG','IGM') {
            foreach my $light ('KAPPA','LAMBDA') {
                my $proper_pair =  join ",", sort ($heavy,$light);
                $result = 1 if $chain_pair eq $proper_pair;
            }
        }
        return $result;
    }


    sub _build_pairing_classification {
        my ($self) = @_;
        my %classif;
        my @heavy = qw/IGA IGG IGM/;
        my @light = qw/KAPPA LAMBDA/;
        $classif{join ',', sort @heavy } = 'THREE_HEAVY_ONLY';
        $classif{join ',', sort @light } = 'TWO_LIGHT_ONLY';
        $classif{join ',', sort (@heavy,@light) } = 'THREE_HEAVY_TWO_LIGHT';
         foreach my $l (@light) {
            my $key = join ",", sort (@heavy, $l);
            $classif{$key} = 'THREE_HEAVY_ONE_LIGHT';
        }
        foreach my $h (@heavy) {
            my $key = join ",", sort (@light,$h);
            $classif{$key} = 'ONE_HEAVY_TWO_LIGHT';
        }
        my $iter = combinations(\@heavy, 2);
        while (my $c = $iter->next) {
            my $key = join ",", sort @{$c};
            $classif{$key} = 'TWO_HEAVY_ONLY';
            # Added to TWO_HEAVY TWO_LIGHT
            $key = join ",", sort (@{$c},@light);
            $classif{$key} = 'TWO_HEAVY_TWO_LIGHT';
            foreach my $l (@light) {
                my $key = join ",", sort (@{$c},$l);
                $classif{$key} = 'TWO_HEAVY_ONE_LIGHT';
            }
        }
        return \%classif;
    }


    sub _get_not_proper_classification {
        my ($self, $chains) = @_;
        my $chain_pair = join ",", sort @{$chains};
        if ($self->chain_pairing_classification->{$chain_pair}){
            return $self->chain_pairing_classification->{$chain_pair};
        }
        else {
            die "Cannot find Pairing Classification for $chain_pair !";
        }
    }


    sub generate_groups {
        my ($self, $pairs_classification) = @_;
        my %pair_types;
        my $i = 0;
        
        foreach my $pair (@{$pairs_classification}) {
            if ($pair) {
                my @chains = keys %{$pair};

                #if ($pair->{'HEAVY'} and scalar @chains == 2){
                if ( $self->_is_proper_paired( \@chains ) ) {
                    foreach my $chain (@chains) {
                        next if $chain =~ /^IG/i;
                        push @{ $pair_types{'PROPER'}{$chain} }, $pair;
                        delete $pairs_classification->[$i];
                    }
                }

                # single chain
                elsif ( scalar @chains == 1 ) {
                    push @{ $pair_types{'SINGLES'}{ $chains[0] } }, $pair;
                    delete $pairs_classification->[$i];
                }
                else {
                    my $classif =
                      $self->_get_not_proper_classification( \@chains );
                    push @{ $pair_types{$classif} }, $pair;
                    delete $pairs_classification->[$i];
                }
            }
            $i++;
        }
        return \%pair_types;
    }


    sub _get_cluster_info {
        my ($self, $clusters, $clusters_by_sample, $type, $info_list) = @_;
        foreach my $info (@{$info_list}) {
            my @V_regions;
            my @VJ_regions;
            my @clone_ids;
            my @sample_ids;
            foreach my $chain (sort {$a cmp $b} keys %{$info}) {
                push @V_regions, $info->{$chain}{V_CALL};
                push @VJ_regions, join ",", @{$info->{$chain}}{('V_CALL','J_CALL')};
                push @clone_ids, $info->{$chain}{clone_id};
                push @sample_ids, $info->{$chain}{SAMPLEID};
            }
            my $strict_cluster_key =  join "|", @VJ_regions;
            my $relaxed_cluster_key =  join "|", @V_regions;
            my $clonal_cluster_key =  join "|", @clone_ids;
            my @SAMPLEID = uniq @sample_ids ;
            die "Something is wrong, chains have different sampleids!" unless scalar @SAMPLEID == 1;
            push @{ $clusters->{'strict'}{$type}{$strict_cluster_key} }, $info;
            push @{ $clusters_by_sample->{$SAMPLEID[0]}{'strict'}{$type}{$strict_cluster_key} }, $info;
            push @{ $clusters->{'relaxed'}{$type}{$relaxed_cluster_key} }, $info;
            push @{ $clusters_by_sample->{$SAMPLEID[0]}{'relaxed'}{$type}{$relaxed_cluster_key} }, $info;
            push @{ $clusters->{'clonal'}{$type}{$clonal_cluster_key} }, $info;
            push @{ $clusters_by_sample->{$SAMPLEID[0]}{'clonal'}{$type}{$clonal_cluster_key} }, $info;
        }
    }

    sub find_clusters {
        my ($self, $pair_types) = @_;
        my %clusters;
        my %clusters_by_sample;
        foreach my $type (sort {$a cmp $b} keys %{$pair_types}) {
            if ($type =~ /PROPER|SINGLES/){
               foreach my $class_chain (keys %{$pair_types->{$type}} ) {
                    my $info_list = $pair_types->{$type}{$class_chain};
                    $self->_get_cluster_info(\%clusters, \%clusters_by_sample, $type, $info_list);
                }
            }
            else {
                my $info_list = $pair_types->{$type};
                $self->_get_cluster_info(\%clusters, \%clusters_by_sample, $type, $info_list);
            }
        }
        return \%clusters, \%clusters_by_sample;
    }


    sub create_douhgnut_chart {
        my ($self, $workbook, $worksheet, $data, $indexed_colors, $styler, $second_chart ) = @_;
        my @chart_colors;
        foreach my $cluster_id (@{$data->[0]}) {
            if ($indexed_colors->{$cluster_id}){
                push @chart_colors, $indexed_colors->{$cluster_id};
            }
        }
        #my $worksheet = $workbook->add_worksheet($worksheet->get_name);
        my $chart     = $workbook->add_chart( type => 'doughnut', embedded => 1 );
        my $total = 0;
        $total += $_ foreach @{$data->[1]}[ 1 .. $#{ $data->[1] }];
        $worksheet->write( 'A1', $data );

        # Configure the chart.
        $chart->add_series(
            categories => [ $worksheet->get_name, 1, $#{ $data->[0] }, 0, 0 ],
            values     => [ $worksheet->get_name, 1, $#{ $data->[1] }, 1, 1 ],
            points     => \@chart_colors,
            border     => { color => 'black' },
        );
        my $title = $worksheet->get_name;
        $title =~ s/\_chart//ig;
       # Insert the chart into the worksheet.
        if ($second_chart){
            $chart->set_title( name => $title . '_SHARED');
            $worksheet->insert_chart( 'K2', $chart, {'x_offset' => 15, 'y_offset' => 7} );
            my $text = $workbook->add_shape(
                type   => 'rect',
                text   => $total,
                width  => 70,
                height => 60,
                line        => 'FFFFFF',
                format => $styler->('floating_text'),
            );
            $worksheet->insert_shape(
                'N9', $text
            );
        }
        else {
            $chart->set_title( name => $title );
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
    }


    sub create_xlsx {
        my ( $self, $clusters, $header, $xls_name ) = @_;

        my @header_first = qw/sampleid cluster_id num_seqs/;
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
            ]
              #808080
              #]
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
            IGA => {
                bg_color => 'silver'
            },
             IGG => {
                bg_color => 'silver'
            },
            IGM => {
                bg_color => 'silver'
            },
 
            KAPPA => {
                bg_color => 'orange'
            },
            LAMBDA => {
                bg_color => 'yellow'
            },
            right_border => {right       => 6,         # double line
                   right_color => 'blue'},

            %color_hash,
            floating_text => {
                align    => "center",
                bold     => 1,
                size => 15
            },

        );

        #my @others = uniq values %{$self->chain_pairing_classification};
        #@others= sort @others;
        my @others = sort( grep(! /PROPER|SINGLES/ , keys %{$clusters}) );
        my @types = ('PROPER', 'SINGLES', @others);

        my $cluster_count = 1;
        foreach my $type (@types) {
            my %ig_proper_isotype;
            my @chart_cluster_id = ('cluster_id');
            my @chart_cluster_size = ('num_seqs');
            my $single_size = 0;
            my @local_colors = @colors;
            my %indexed_colors;
            my %shared_clusters;
            my $c_type = $clusters->{$type};
            my $i = 0;
            # Add a worksheet
            my $worksheet = $workbook->add_worksheet($type);
            my $worksheet_chart = $workbook->add_worksheet($type.'_CHART');
            if ($type =~ /PROPER/){
                foreach my $isotype (qw/IGA IGM IGG/) {
                    $ig_proper_isotype{$isotype}{worksheet} = $workbook->add_worksheet($type.'_'.$isotype);
                    $ig_proper_isotype{$isotype}{chart} = $workbook->add_worksheet($type.'_'.$isotype.'_CHART');
                }
            }
            # add header
            #my @chains = qw/HEAVY KAPPA LAMBDA/;
            my @chains = qw/IGA IGM IGG KAPPA LAMBDA/;
            # add merged columns
            my $col_size = scalar @{$header};
            my $col_end = $col_size + 1;
            my $col_start = 2; # always start on 2
            foreach my $chain (@chains) {
                $worksheet->merge_range($i,$col_start,$i,$col_end,$chain, $styler->('title',$chain));
                # add to the isotype worksshets if they exist
                foreach my $isotype (sort {$a cmp $b}keys %ig_proper_isotype) {
                    $ig_proper_isotype{$isotype}{worksheet}->merge_range($i,$col_start,$i,$col_end,$chain, $styler->('title',$chain));
                }
                $col_start += $col_size;
                $col_end += $col_size;
            }

            $i++;
            $worksheet->write_row($i,0,[ @header_first, @{$header}, @{$header}, @{$header},@{$header},@{$header}], $styler->(qw/title/));
            # add to the isotype worksshets if they exist
            foreach my $isotype (sort {$a cmp $b} keys %ig_proper_isotype) {
                $ig_proper_isotype{$isotype}{worksheet}->write_row($i,0,[ @header_first, @{$header}, @{$header}, @{$header},@{$header},@{$header}], $styler->(qw/title/));
                $ig_proper_isotype{$isotype}{row} = 2;
            }
            $i++;

            # Recalculate clones by isotype
            foreach my $cluster_key ( sort { $#{ $c_type->{$b} } <=> $#{ $c_type->{$a} } || $a cmp $b } keys %{$c_type} ) {
                my $cluster_size = scalar @{ $c_type->{$cluster_key} };
                foreach my $info ( @{ $c_type->{$cluster_key} } ) {
                    foreach my $chain (@chains) {
                        if ( $info->{$chain} ) {
                            if (%ig_proper_isotype and $chain =~ /^IG/){
                                $ig_proper_isotype{$chain}{clusters}{$cluster_key}{size}++;
                            }
                        }
                    }
                }
            }

            foreach my $cluster_key ( sort { $#{ $c_type->{$b} } <=> $#{ $c_type->{$a} } || $a cmp $b } keys %{$c_type} ) {
                my $cluster_size = scalar @{ $c_type->{$cluster_key} };

                my $color;
                if ( $cluster_size > 1 ) {
                    $color = shift @local_colors;
                    push @local_colors, $color;
                    push @chart_cluster_id, $cluster_count;
                    push @chart_cluster_size, $cluster_size;
                    $indexed_colors{$cluster_count} = { fill => { color => '#'.$color }, border => {color => '#000000'} };
                }
                else {
                    $single_size++;
                    $indexed_colors{single} = { fill => { color => '#FFFFFF' }, border => {color => '#000000'}  }
                }

                my $current_isotype; # keep current isotype

                foreach my $info ( @{ $c_type->{$cluster_key} } ) {
                    my $SAMPLEID;
                    $current_isotype = undef; # keep current isotype
                    # Assume that all chains have the same SAMPLEID
                    foreach my $chain (@chains) {
                        $SAMPLEID = $info->{$chain}->{SAMPLEID} if $info->{$chain};
                    }
                    my @columns = ( $SAMPLEID, $cluster_count, $cluster_size );
                    #foreach my $chain ( sort { $a cmp $b } keys %{$info} ) {
                    foreach my $chain (@chains) {
                        if ( $info->{$chain} ) {
                            if (%ig_proper_isotype and $chain =~ /^IG/){
                                $current_isotype = $chain;
                           }
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
                        $worksheet->write_row( $i, 0, \@columns, $styler->($color) );
                        if ($current_isotype){
                            # Change cluster number
                            $columns[2] = $ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{size};
                            $ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{count}++;
                            $ig_proper_isotype{$current_isotype}{worksheet}->write_row( $ig_proper_isotype{$current_isotype}{row}, 0, \@columns, $styler->($color) );
                        }
                        #$worksheet->write( $i, 0, $columns[0], $styler->($color,'right_border'));
                    }
                    else {
                        $worksheet->write_row( $i, 0, \@columns );
                        if ($current_isotype){
                            $ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{count}++;
                            $ig_proper_isotype{$current_isotype}{worksheet}->write_row($ig_proper_isotype{$current_isotype}{row}, 0, \@columns);
                        }
                    }

                    if ($current_isotype){
                        $ig_proper_isotype{$current_isotype}{row}++;
                        $ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{cluster_id} = $cluster_count;
                        if ($ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{count} == $ig_proper_isotype{$current_isotype}{clusters}{$cluster_key}{size}){
                            $ig_proper_isotype{$current_isotype}{worksheet}->write_row($ig_proper_isotype{$current_isotype}{row}, 0, []);
                            $ig_proper_isotype{$current_isotype}{row}++;
                        }
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
            }
            my @data = (\@chart_cluster_id,\@chart_cluster_size);
            $self->create_douhgnut_chart($workbook, $worksheet_chart, \@data, \%indexed_colors, $styler);

            # Calculate Shared clusters
            foreach my $isotype (sort {$a cmp $b} keys %ig_proper_isotype ) {
                my $clusters = $ig_proper_isotype{$isotype}{clusters};
                foreach my $cluster_key ( sort { $clusters->{$b}->{size} <=> $clusters->{$a}->{size} || $a cmp $b } keys %{$clusters} ) {
                    if ($clusters->{$cluster_key}{size}){
                        $shared_clusters{$clusters->{$cluster_key}{cluster_id}}{$isotype}++;
                    }
                }
            }

            my %shared_colors = %indexed_colors;
            foreach my $cluster_id (keys %shared_clusters) {
                my $n_shared_isotypes = scalar keys %{$shared_clusters{$cluster_id}};
                if ($n_shared_isotypes == 1){
                    $shared_colors{$cluster_id} = { fill => { color => '#CACACA' }, border => {color => '#000000'}  } 
                }
            }

            # Add isotype charts
            foreach my $isotype (sort {$a cmp $b} keys %ig_proper_isotype ) {
                my $clusters = $ig_proper_isotype{$isotype}{clusters};
                my @isotype_data_cluster_id = ('cluster_id');
                my @isotype_data_size = ('num_seqs');
                my $isotype_single_size;
                foreach my $cluster_key ( sort { $clusters->{$b}->{size} <=> $clusters->{$a}->{size} || $a cmp $b } keys %{$clusters} ) {

                    if ($clusters->{$cluster_key}{size} > 1){
                        push @isotype_data_cluster_id, $clusters->{$cluster_key}{cluster_id};
                        push @isotype_data_size, $clusters->{$cluster_key}{size};
                    }
                    else{
                        $isotype_single_size++;
                    }
                }

                if ($isotype_single_size) {
                    push @isotype_data_cluster_id,   'single';
                    push @isotype_data_size, $isotype_single_size;
                }

                my @isotype_data = (\@isotype_data_cluster_id, \@isotype_data_size);
                $self->create_douhgnut_chart($workbook, $ig_proper_isotype{$isotype}{chart}, \@isotype_data, \%indexed_colors, $styler);
                $self->create_douhgnut_chart($workbook, $ig_proper_isotype{$isotype}{chart}, \@isotype_data, \%shared_colors, $styler, 'TRUE');
            }

        }
    }


    sub run {
        my ($self) = @_;
        my ( $index, $pairs ) = $self->parse_pairs_files;
        my $seq_classification = $self->parse_clone_summary($index);
        my $pairs_classification = $self->make_pairs( $seq_classification, $pairs );
        my $pair_types       = $self->generate_groups($pairs_classification);
        my ($clusters, $clusters_by_sample) = $self->find_clusters($pair_types);
        my @expected_header = @TABLE_HEADER[ 2 .. $#TABLE_HEADER ];

        my @selected = qw/SEQUENCE_ID V_CALL D_CALL J_CALL corrected_input_from_start translation_of_corrected_input_from_start v_insertions v_deletions nt_mismatches_V_region v_region_aa_mismatches  cdr3_aa cdr3_aa_length/;
        foreach my $sampleid ( sort { $a cmp $b } keys %{$clusters_by_sample} ) {
            foreach my $type (qw/strict clonal/) {
                $self->create_xlsx( $clusters_by_sample->{$sampleid}{$type}, \@expected_header,
                    $self->output_prefix . '_'.$sampleid.'_' . $type . '.xlsx' );
                $self->create_xlsx( $clusters_by_sample->{$sampleid}{$type}, \@selected,
                        $self->output_prefix
                      . '_'. $sampleid . '_'
                      . $type
                      . '_selected_columns.xlsx' );
            }

        }

        #create combo
        foreach my $type (qw/strict clonal/) {
            $self->create_xlsx($clusters->{$type}, \@expected_header, $self->output_prefix . '_all_samples_' . $type . '.xlsx');
            $self->create_xlsx($clusters->{$type}, \@selected,  $self->output_prefix . '_all_samples_' . $type . '_selected_columns.xlsx');
        }
   }


    __PACKAGE__->meta->make_immutable;
}


use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

#!/usr/bin/perl

# perl script to parse schematic reaction files.

# Initialise all variables.  The large block of text below is 
# documentation printed when requested.

sub InitAll {
    $codename = "0";
    $reaction_file = "$codename.reactions";
    $component_file = "$codename.components";
    @name_list = ();
    @init_vals = ();
    %components = ();
    $num_eq_blocks = 100;
    $num_prod_blocks = 0;
    $num_steps = 10000;
    $freq_anal = 100;
    $help = <<EOF;

Parses schematic reaction files and generates an input file for Gillespie
on stdout.  Schematic reaction files are parsed line by line, each line is
one of the following forms.

1. Blank lines and lines not containing ":" are ignored.

2. All text after '#' is ignored.

3. A line of the type
 <value> : <identifier> 
resets default values.  The <identifier> is from the following table

-------------------------------------------------------------------------
 <identifier>          minimum components in <identifier>  default
-------------------------------------------------------------------------
 number of equilibrium block   "num" "eq"                   50
 number of production blocks   "num" "prod"                 0
 number of steps               "num" "steps"                500
 frequency analysis            "freq"                       10 
 code name                     "name"                       0
 reaction file                 "react" "file"               0.reactions
 component file                "comp" "file"                0.components
-------------------------------------------------------------------------

Changing the code name also has the effect of changing the file names.

4. A line of the type
 <space_separated_list_of_components> : components
declares components, and is used to determine a particular order of
appearance, otherwise components appear in the order in which they
appear in reactions.

5. A line of the type
 <component> <equals> <val> : init 
causes <component> name to have a initial value <val> (default 0).  
The separator <equals> can be anything.

6. Lstly a line of the type 
 <reactants> <transform> <products> : <k_values>
declares a reaction with the following grammar.

Annihilation reactions (can write "null" or "NULL"),
 A --> null : k = <val>
 A + B --> null : k = <val>

Simple reactions,
 A + B --> C + D + E : k = <val>

Production reactions, in which the reactants are not destroyed,
 A + B +--> C + D : k = <val>
is equivalent to
 A + B --> A + B + C + D : k = <val>

Equilibriums, which are converted to reaction pairs,
 A + B <--> C + D : k = <val_forward> : k = <val_backward>
is equivalent to
 A + B --> C + D : k = <val_forward>
 C + D --> A + B : k = <val_backward>

Finally simple and production reactions can be written with reactants
and products as lists of alternatives.  In such cases all possible
pairs of alternatives are implemented, with the same k, thus:
 A | B | C --> null : k = <val_1>
 A + B +--> C | D : k = <val_2>
is equivalent to
 A --> null : k = <val_1>
 B --> null : k = <val_1>
 C --> null : k = <val_1>
 A + B --> A + B + C : k = <val_2>
 A + B --> A + B + D : k = <val_2>

At the end, reaction and component files suitable for direct reading
by Gillespie code are generated.  Standard output can be redirected to
Gillespie.inp (the default avoids accidental overwriting).  Various
comments are generated on standard error.

EOF
}

# Returns id number of a component of a given name, either one
# previously used or a new one which is initialised appropriately.

sub component_id {
    local ($name, $id);
    $name = $_[0];
    $name =~ s/\s//g;
    if (exists $components{$name}) {
        $id = $components{$name};
    }
    else {
        $id = scalar(@name_list);
        $components{$name} = $id;
        @init_vals[$id] = 0;
        @name_list[$id] = $name;
    }
    return $id;
}

# Add a list of components, in the order in which they appear.

sub add_components {
    local ($name);
    @names = split(' ', $_[0]);
    foreach $name (@names) {
        component_id($name);
    }
}

# Set or reset the initial concentration for a given component.

sub set_initial {
    local ($name, $expr, $val, $id);
    ($name, $expr, $val) = split(' ', $_[0]);
    printf STDERR "Initialising component %s to %i\n", $name, $val;
    $id = component_id($name);
    @init_vals[$id] = $val;
}

# Add a reaction, specifically compute the output lines for later
# writing.

sub add_reaction {
    local ($reacts, $prods, $kval, $component, $id, $nreacts, $nprods);
    printf STDERR "Adding reaction %s --> %s, k = %f\n", $_[0], $_[1], $_[2];
    ($reacts, $prods, $kval) = ($_[0], $_[1], $_[2]);
    $kval =~ s/\s//g;
    $line2 = "";
    @reacts = split(/\+/, $reacts);
    $nreacts = scalar(@reacts);
    $addmore = 0;
    foreach $component (@reacts) {
        if ($addmore) {$line2 .= "+ ";}
        $id = component_id($component);
        $bit = sprintf("X %i ", $id);
        $line2 .= $bit;
        $addmore = 1;
    }
    $line2 .= "->";
    if ($prods =~ /null/ || $prods =~ /NULL/) {
        $line2 .= " 0";
        $nprods = 0;
    }
    else {
        @prods = split(/\+/, $prods);
        $nprods = scalar(@prods);
        $addmore = 0;
        foreach $component (@prods) {
            if ($addmore) {$line2 .= " +";}
            $id = component_id($component);
            $bit = sprintf(" 1 X %i", $id);
            $line2 .= $bit;
            $addmore = 1;
        }
    }
    $line1 = sprintf("%f\t%i\t%i\tRateConstant_k_Nreactants_Nproducts",
        $kval, $nreacts, $nprods);
    push(@react1lines, $line1);
    push(@react2lines, $line2);
}

# Parse a single line from a schematic reaction file.

sub parseline {
    local ($line, $front, $back, $reacts, $prods);
    $line = $_[0];
    $line =~ s/#.*//;
    if ($line =~ /:/) {
        ($front, $back) = split(":", $line, 2);
        if ($back =~ /name/) {
            $codename = $front;
            $codename =~ s/\s//g;
            $reaction_file = "$codename.reactions";
            $component_file = "$codename.components";
            print STDERR "Setting codename to $codename\n";
        }
        elsif ($back =~ /file/) {
            if ($back =~ /react/) {
                $reaction_file = $front;
                $reaction_file =~ s/\s//g;
                print STDERR "Setting reaction file to $reaction_file\n";
            }
            elsif ($back =~ /comp/) {
                $component_file = $front;
                $component_file =~ s/\s//g;
                print STDERR "Setting component file to $component_file\n";
            }
        }
        elsif ($back =~ /num/) {
            if ($back =~ /eq/) {
                $num_eq_blocks = $front;
                $num_eq_blocks =~ s/\s//g;
                printf STDERR "Setting # equilibration blocks to %i\n",
                    $num_eq_blocks;
            }
            elsif ($back =~ /prod/) {
                $num_prod_blocks = $front;
                $num_prod_blocks =~ s/\s//g;
                printf STDERR "Setting # production blocks to %i\n",
                    $num_prod_blocks;
            }
            elsif ($back =~ /step/) {
                $num_steps = $front;
                $num_steps =~ s/\s//g;
                printf STDERR "Setting # steps to %i\n", $num_steps;
            }
        }
        elsif ($back =~ /freq/) {
            $freq_anal = $front;
            $freq_anal =~ s/\s//g;
            printf STDERR "Setting frequency analysis to %i\n", $freq_anal;
        }
        elsif ($back =~ /init/) {
            set_initial($front);
        }
        elsif ($back =~ /comp/) {
            add_components($front);
        }
        elsif ($front =~ /<-->/) {
            ($optreacts, $optprods) = split('<-->', $front);
            @optreacts = split(/\|/, $optreacts);
            @optprods = split(/\|/, $optprods);
            if ($back =~ /:/) {
                ($kforward, $kback) = split(/:/, $back, 2);
            }
            else {
                $kforward = $back;
                $kback = $back;
            }
            $kforward =~ s/.*=\s*//;
            $kback =~ s/.*=\s*//;
            foreach $reacts (@optreacts) {
                foreach $prods (@optprods) {
                    add_reaction($reacts, $prods, $kforward);
                    add_reaction($prods, $reacts, $kback);
                }
            }
        }
        elsif ($front =~ /\+-->/) {
            ($optreacts, $optprods) = split('\+-->', $front);
            @optreacts = split(/\|/, $optreacts);
            @optprods = split(/\|/, $optprods);
            $kval = $back;
            $kval =~ s/.*=\s*//;
            foreach $reacts (@optreacts) {
                foreach $prods (@optprods) {
                    add_reaction($reacts, "$reacts + $prods", $kval);
                }
            }
        }
        elsif ($front =~ /-->/) {
            ($optreacts, $optprods) = split('-->', $front);
            @optreacts = split(/\|/, $optreacts);
            @optprods = split(/\|/, $optprods);
            $kval = $back;
            $kval =~ s/.*=\s*//;
            foreach $reacts (@optreacts) {
                foreach $prods (@optprods) {
                    add_reaction($reacts, $prods, $kval);
                }
            }
        }
        else {
            printf STDERR "Didn't recognise %s\n", $line;
        }
    }
}

# Parse a schematic reaction file line by line.

sub ParseFile {
    $in_file = $_[0];
    open(FIN, $in_file) || die "Cannot open $in_file: $!";
    printf STDERR "Reading from $in_file\n";
    while ($line = <FIN>) {
        chop($line);
        parseline($line);
    }
    close(FIN);
}

# Write a list of current reactions to a file using PRs format.

sub WriteReactionFile {
    $out_file = $_[0];
    open(FOUT, "> $out_file") || die "Cannot open $out_file: $!";
    $nreact = scalar(@react1lines);
    printf FOUT "%i\t\tNumber_of_reaction_channels\n\n", $nreact;
    for ($i = 0; $i < $nreact; $i++) {
        printf FOUT "%s\n%s\n", @react1lines[$i], @react2lines[$i];
    }
    print FOUT "\n";
    close(FOUT);
    printf STDERR "Written %i reactions to $out_file\n", $nreact;
}

# Write a list of current components to a file using PRs format.

sub WriteComponentFile {
    local ($ncomp, $id);
    $out_file = $_[0];
    open(FOUT, "> $out_file") || die "Cannot open $out_file: $!";
    $ncomp = scalar(@name_list);
    printf FOUT "%i\n", $ncomp;
    for ($id = 0; $id < $ncomp; $id++) {
        printf FOUT "%i\t\t%s\n", @init_vals[$id], @name_list[$id];
    }
    print FOUT "\n";
    close(FOUT);
    printf STDERR "Written %i components to $out_file\n", $ncomp;
}

# Write a control file to stdot using current information, suitable
# for PRs Gillespie code.  It is intended this should be re-directed
# as appropriate.

sub WriteGillespieInp {
    printf "%i\t\t\tName\n", $codename;
    printf "%i\t\t\tNumber_of_components\n", scalar(@name_list);
    printf "%i\t\t\tNumber_of_reactions\n", scalar(@react1lines);
    printf "%i\t%i\t%i\tNumber_eq_blocks;_Number_prod_blocks;_Number_steps\n",
        $num_eq_blocks, $num_prod_blocks, $num_steps;
    printf "%i\t\t\tFrequency_analysis\n\n", $freq_anal;
}

# This is where the action starts.

InitAll;

if (scalar(@ARGV) == 0) {
    printf STDERR "Usage: ./parse.pl <input_files> > Gillespie.inp\n";
    printf STDERR "./parse.pl --help for more info\n";
    exit(0);
}

if (@ARGV[0] =~ /--help/) {
    print "Usage: ./parse.pl <input_files> > Gillespie.inp\n";
    print $help;
    exit(0);
}

foreach $file (@ARGV) {
    ParseFile($file);
}

WriteReactionFile($reaction_file);
WriteComponentFile($component_file);
WriteGillespieInp;

# End of file.

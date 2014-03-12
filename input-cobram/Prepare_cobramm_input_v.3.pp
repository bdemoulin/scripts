#!/usr/bin/perl -w

# routine needed to check if atom names the same in real.pdb and resid.mol2 files
# make the script search for all solvent molecules in a radius and include them in the high layer

$ENV{AMBERHOME}="/usr/local/src/amber11";
system('export PATH=$AMBERHOME/exe:$PATH');
use strict;
use warnings;
no warnings 'recursion';
use POSIX qw/floor/;
use POSIX qw/ceil/;
use Cwd;
use Env;
use feature 'switch';

my ($BASEDIR,$input,@PDB,@RES,@high_layer,@medium_layer,@medium_layer1,@low_layer,$type,@hlinks,@mlinks,@llinks,$new_resid);
my ($resid,@UNIQUE_RES,$solvent,$solvent_name,@solvent,$dist,@link_number,$AMBERHOME,$BABELHOME,@alphabet);
my (@nonstd_resid,$nonstd_resid,@std_resid,@ff,@predef_resid,@model_resid,%PTE,$solshell,@REASSIGNED_RES);

$BASEDIR    = &getcwd();
$input      = "$BASEDIR/cobram.parm";
$AMBERHOME  = $ENV{'AMBERHOME'};
$BABELHOME  = $ENV{'BABELHOME'};
@high_layer = ();
@medium_layer = ();
@low_layer = ();
$solvent_name = 'none';
$solshell = 'CoMH';
@alphabet = qw/A B C D E F G H I J K L M N O P Q R S T U V W X Y Z/;
%PTE = (H  => [1 ,  1.007940], He => [2 ,  4.002602], Li => [3 ,  6.941000], Be => [4 ,  9.012182], B  => [5 , 10.811000], C  => [6 , 12.010700], N  => [7 , 14.006700], O  => [8 , 15.999400], F  => [9 , 18.998403], Ne => [10, 20.179700], Na => [11, 22.989770], Mg => [12, 24.305000], Al => [13, 26.981539], Si => [14, 28.085500], P  => [15, 30.973762], S  => [16, 32.065000], Cl => [17, 35.453000], Ar => [18, 39.948000], K => [19, 39.0983], Ca => [20, 40.078], Se => [34, 78.96 ], Br => [35, 79.904], I => [53, 126.90447 ]);

############### Main ################
#####################################

read_coord();
get_high_layer(); 
get_medium_layer(); 
get_low_layer();
sol_specs();
get_links();
prepare_real_layers();
get_connection();
my $ff = read_input("force fields");
@ff = split(/\s+/,$ff);
if ($nonstd_resid ne '0'){
	prepare_nonstd();
}
prepare_model();
prepare_real();
check_atom_types();
keep_files();

#####################################
########### Subroutines #############
#####################################

sub read_input{
my $key = $_[0];
open (INP, "<".$input) or die;
while (<INP>){
	if (/^$key/){
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/=/,$_);
		if ($_[1]){
                	$_[1] =~ s/^\s*//;
			close(INP);
                	return $_[1];
		} else {
			close(INP);
			return 0;
		}
	}
}
}

#####################################

sub sol_specs{
$solvent = read_input("solvent");
if ($solvent =~ /(\w+)\s+(\d+[\.\d+]*)\s+(\w+)/) {
        $solvent_name = $1;
        $solvent = $2;
        $solshell = $3;
        include_solvent($solvent);
} elsif ($solvent =~ /(\w+)\s+(\d+[\.\d+]*)/) {
        $solvent_name = $1;
        $solvent = $2;
        include_solvent($solvent);
} elsif ($solvent =~ /(\w+)\s+(\w)/){
        $solvent_name = $1;
        $solshell = $2;
        get_H_bonds();
} elsif ($solvent =~ /(\w+)/){
        $solshell = 'M';
        $solvent_name = $1;
        get_H_bonds();
}
}

#####################################

sub read_coord{
my $old_type;
my $coord = read_input("coordinate file");
open (PDB, "<".$coord) or die "Input file $coord is missing!";
print "Reading $coord ... ";
$resid = 0;
$old_type = '';
while (<PDB>){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
	$_ =~ s/-/ -/g;
        @_ = split(/\s+/,$_);
        if ($_[0] eq "ATOM"){
		if ($_[1] == 1 || $old_type eq 'TER'){$old_type = $_[4];}
                push(@{$PDB[0]},$_[3]);
                push(@{$PDB[1]},$_[5]);
                push(@{$PDB[2]},$_[6]);
                push(@{$PDB[3]},$_[7]);
		$_[2] =~ s/\d+//g;
		until (grep {$_ eq $_[2]} keys %PTE){
	                $_[2] =~ s/.$//g;
		}
		push(@{$PDB[4]},$_[2]);
		if ($_[4] != $old_type){
			$resid++;
		}
	        push(@{$RES[$resid]},$_[1]);
		$old_type = $_[4];
        }
        if ($_[0] eq "TER"){
                $resid++;
		$old_type = 'TER';
        }
}
#for my $i (0..$#RES){
#	print "Residue: $i:";
#	for my $j (0..$#{$RES[$i]}){
#		print " ${$RES[$i]}[$j] ";
#	}
#	print "\n";
#}
close (PDB);
print "done!\n";
}

#####################################

sub get_high_layer{
my @range;
$_ = read_input("high layer atoms");
if ($_ eq '0'){
	return;
}
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        @range = split(/-/,$_[$i]);
	if (!$range[1]){
		push(@high_layer,$range[0]);	
	} else {
        	for my $j ($range[0]..$range[1]){
                	push(@high_layer,$j);
		}
        }
	undef @range;
}
}

######################################

sub get_medium_layer{
$_ = read_input("medium layer atoms");
if ($_ eq '0'){
        return;
}
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        my @range = split(/-/,$_[$i]);
        for my $j ($range[0]..$range[1]){
                push(@medium_layer,$j);
        }
}
@medium_layer1 = @medium_layer;
}

#####################################

sub get_low_layer{
for my $j (1..$#{$PDB[0]}+1){
        if ((grep {$_ == $j} @medium_layer) || (grep {$_ == $j} @high_layer)){
        } else {
                push(@low_layer,$j);
        }
}
}

#####################################

sub get_H_bonds{
my @heteroatoms = ("N", "O", "F", "P", "S", "Cl", "Br");
my $dist;
print "Looking for H-bonds ... \n";
foreach my $i (@high_layer,@medium_layer){
	if (grep {$_ eq ${$PDB[4]}[$i-1]} @heteroatoms){
		foreach my $j (@low_layer){
			if (${$PDB[4]}[$j-1] eq "H" && ${$PDB[0]}[$j-1] eq $solvent_name){
				$dist = calculate_dist($i,$j);
				if ($dist <= 2.5){
					foreach my $k (@low_layer){
						if ($k == $j){
						} else {
							$dist = calculate_dist($k,$j);
							if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
								print "Heteroatom: $i(${$PDB[4]}[$i-1])\n";
								push(@solvent,complete_solvent($j));
							}
						}
					}
				}
			}
		}
	}
}
print "done for heteroatoms in the high/medium layer (Het(H/M)---H-Het(L))\n\n";
print "Looking for H-bonds ... \n";
foreach my $i (@high_layer){
        if (${$PDB[4]}[$i-1] eq "H"){
		foreach my $k (@high_layer){
			if ($k == $i){
			} else {
				$dist = calculate_dist($i,$k);
				if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
					foreach my $j (@low_layer){
                 			        if (grep {$_ eq ${$PDB[4]}[$j-1]} @heteroatoms){
                                			$dist = calculate_dist($i,$j);
                                			if ($dist <= 2.5){
								print "Hydrogen: $i(${$PDB[4]}[$i-1]) ... $k(${$PDB[4]}[$k-1])\n";
								push(@solvent,complete_solvent($j));
                                			}
                        			}
                			}
				}
			}
                }
        }
}
foreach my $i (@medium_layer){
        if (${$PDB[4]}[$i-1] eq "H"){
                foreach my $k (@medium_layer){
                        if ($k == $i){
                        } else {
                                $dist = calculate_dist($i,$k);
                                if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
                                        foreach my $j (@low_layer){
                                                if (grep {$_ eq ${$PDB[4]}[$j-1]} @heteroatoms){
                                                        $dist = calculate_dist($i,$j);
                                                        if ($dist <= 2.5){
								print "Hydrogen: $i(${$PDB[4]}[$i-1]) ... $k(${$PDB[4]}[$k-1])\n";
								push(@solvent,complete_solvent($j));
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}
print "done for H's attached to heteroatoms in the high/medium layer (Het(H/M)-H(H/M)---Het(L)\n\n";
if ($solshell eq 'H'){
	move_to_high_layer(@solvent);
} elsif ($solshell eq 'M') {
	move_to_medium_layer(@solvent);
}
print "High layer atoms: ";
for my $i (@high_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
print "Medium layer atoms: ";
for my $i (@medium_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
#print "Low layer atoms: ";
##for my $i (@low_layer){
##        print "$i(${$PDB[0]}[$i-1]) ";
##}
##print "\n\n";
}

#####################################

sub include_solvent {
my $r = $_[0];
if (($solshell eq 'CoMH') || ($solshell eq 'CoMHM')){
	my @R = center_of_mass($solshell);
	print "Center of mass: @R\n";
	push(@{$PDB[1]},$R[0]);
	push(@{$PDB[2]},$R[1]);
	push(@{$PDB[3]},$R[2]);
	foreach my $j (@low_layer){
		if (${$PDB[0]}[$j-1] eq $solvent_name){
			$dist = calculate_dist($#{$PDB[1]}+1,$j);
			if ($dist <= $r){
				print "Include low layer atom: $j(${$PDB[4]}[$j-1])\n";
				push(@solvent,complete_solvent($j));
			}
		}
	}
} elsif ($solshell eq 'DtH'){
	foreach my $i (@high_layer){
	        foreach my $j (@low_layer){
		        $dist = calculate_dist($i,$j);
		        if ($dist <= $r){
			        print "High layer atom: $i(${$PDB[4]}[$i-1]) ... Low layer atom: $j(${$PDB[4]}[$j-1])\n";
			        push(@solvent,complete_solvent($j));
	                }
	        }
	}
} elsif ($solshell eq 'DtHM'){
        foreach my $i (@high_layer,@medium_layer){
                foreach my $j (@low_layer){
                        $dist = calculate_dist($i,$j);
                        if ($dist <= $r){
                                print "High/Medium layer atom: $i(${$PDB[4]}[$i-1]) ... Low layer atom: $j(${$PDB[4]}[$j-1])\n";
                                push(@solvent,complete_solvent($j));
                        }
                }
        }
}
move_to_medium_layer(@solvent);
print "High layer atoms: ";
for my $i (@high_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
print "Medium layer atoms: ";
for my $i (@medium_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
}

#####################################

sub center_of_mass{
my @R = 0;
my $M = 0;
if ($_[0] eq 'CoMH'){
	@_ = @high_layer;
} elsif ($_[0] eq 'CoMHM'){
	@_ = (@high_layer,@medium_layer);
}
foreach my $i (@_){
	$M += ${$PTE{${$PDB[4]}[$i-1]}}[1];
}
foreach my $i (@_){
	$R[0] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[1]}[$i-1]/$M;
	$R[1] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[2]}[$i-1]/$M;
	$R[2] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[3]}[$i-1]/$M;
}
return @R;
}

#####################################

sub calculate_dist{
my $a = $_[0];
my $b = $_[1];
$dist = sqrt((${$PDB[1]}[$a-1] - ${$PDB[1]}[$b-1])**2 + (${$PDB[2]}[$a-1] - ${$PDB[2]}[$b-1])**2 + (${$PDB[3]}[$a-1] - ${$PDB[3]}[$b-1])**2);
return $dist;
}

#####################################

sub complete_solvent{
my $atom = $_[0];
my @molecule = ();
print "Solvent atom: $atom(${$PDB[4]}[$atom-1])\n";
for my $i (0..$#RES){
	if (grep {$_ == $atom} @{$RES[$i]}){
		for my $j (@{$RES[$i]}){
			push(@molecule,$j);
		}
		last;
	}
}
print "Solvent molecule: ";
foreach my $i (@molecule){
	print "$i(${$PDB[4]}[$i-1]) ";
}
print "\n";
return @molecule;
}

#####################################

sub move_to_high_layer{
my %dummy = ();
my @solvent = grep { ! $dummy{ $_ }++ } @_; # removes duplicated solvent molecules
printf "Number of solvent molecules that will be added to the high layer:%4d\n\n", ($#solvent+1)/($#{$RES[-1]}+1);
foreach my $i (@solvent){
        push(@high_layer,$i);
        my $count = 0;
        foreach my $j (@low_layer){
                if ($j == $i){
                        last;
                }
                $count++;
        }
        splice(@low_layer,$count,1);
}
@high_layer = sort {$a <=> $b} @high_layer;
@low_layer = sort {$a <=> $b} @low_layer;
}

#####################################

sub move_to_medium_layer{
my %dummy = ();
my @solvent = grep { ! $dummy{ $_ }++ } @_; # removes duplicated solvent molecules
printf "Number of solvent molecules that will be added to the medium layer:%4d\n\n", ($#solvent+1)/($#{$RES[-1]}+1);
foreach my $i (@solvent){
	push(@medium_layer,$i);
	my $count = 0;
	foreach my $j (@low_layer){
		if ($j == $i){
			last;
		}
		$count++;
	}
	splice(@low_layer,$count,1);
}
@medium_layer = sort {$a <=> $b} @medium_layer;
@low_layer = sort {$a <=> $b} @low_layer;
}

#####################################

sub get_links{
my $link_atoms = read_input("links");
if ($link_atoms eq '0'){
        return;
}
@_ = split(/\s+/,$link_atoms);
for my $i (0..$#_){
        my @b_atoms = split(/-/,$_[$i]);
        if (grep {$_ == $b_atoms[0]} @high_layer){
		if (grep {$_ == $b_atoms[1]} @medium_layer){
	                push(@hlinks,$b_atoms[0]);
                	push(@mlinks,$b_atoms[1]);
		} elsif (grep {$_ == $b_atoms[1]} @low_layer){
			push(@hlinks,$b_atoms[0]);
                        push(@llinks,$b_atoms[1]);
		} else {
			print "Link atom $b_atoms[0] was found in the high layer but link atom $b_atoms[1] was found neither in the medium layer nor in the low layer!\n"; exit;
		}
        } elsif (grep {$_ == $b_atoms[0]} @medium_layer){
		if (grep {$_ == $b_atoms[1]} @high_layer){
                        push(@hlinks,$b_atoms[1]);
                        push(@mlinks,$b_atoms[0]);
		} else {
			print "Link atom $b_atoms[0] was found in the medium layer but link atom $b_atoms[0] was not found in the high layer!\n"; exit;
		}
        } elsif (grep {$_ == $b_atoms[0]} @low_layer){
                if (grep {$_ == $b_atoms[1]} @high_layer){
                        push(@hlinks,$b_atoms[1]);
                        push(@llinks,$b_atoms[0]);
                } else {
                        print "Link atom $b_atoms[0] was found in the low layer but link atom $b_atoms[0] was not found in the high layer!\n"; exit;
                }
	} else {
                print "Link atom $b_atoms[0] wasn't found in any layers!\n"; exit;
        }
}
if (@mlinks){
	print "Link atoms:    H        M\n";
	foreach my $i (0..$#mlinks){
		print "              $hlinks[$i](${$PDB[4]}[$hlinks[$i]-1])    $mlinks[$i](${$PDB[4]}[$mlinks[$i]-1])\n";
	}
}
if (@llinks){
	print "Link atoms:    H        L\n";
	foreach my $i (0..$#llinks){
	        print "              $hlinks[$i](${$PDB[4]}[$hlinks[$i]-1])    $llinks[$i](${$PDB[4]}[$llinks[$i]-1])\n";
	}
}
print "\n";
}

########################################

sub prepare_real_layers{
my (@dist,$r,@H_pos,@resid_size);
my $default;
$default->{'C'} = 1.089;
$default->{'N'} = 1.008;
$default->{'O'} = 0.947;
$nonstd_resid = read_input("non-standard residues");
if ($nonstd_resid eq '0'){
       print "Only standard residues used throughout\n";
} else {
	my $nonstd_param = read_input("non-standard residue parameters");
	my $nonstd_crg = read_input("non-standard residue charges");
	@nonstd_resid = split(/\s+/,$nonstd_resid);
	my @nonstd_param = split(/\s+/,$nonstd_param);
	my @nonstd_crg = split(/\s+/,$nonstd_crg);
	printf "%2d non-standard residues will be generated\n",scalar(@nonstd_resid);
	print "Name  Charge  Model\n";
	for my $i (0..$#nonstd_resid){
		if ($nonstd_resid[$i] =~ /\w{4,}/){
	                print "Non-standard residue names must have a 3-letter code!\n";
	                exit;
		}
		print "  $nonstd_resid[$i]    $nonstd_crg[$i]    $nonstd_param[$i]\n";
		push(@resid_size,"0");
	}
}
print "Generating the real_layers.xyz and performing H/M/L layer sorting!\n";
open (MAP, ">pdb_2_HL.map") or die;
open (OUT, ">real_layers.xyz") or die;
open (XYZ, ">model-H.xyz") or die;
printf XYZ "%6d\n\n",scalar(@high_layer)+scalar(@hlinks);
for my $i (1..$#{$PDB[0]}+1){
	for my $j (0..$#nonstd_resid){
		if (${$PDB[0]}[$i-1] eq $nonstd_resid[$j]){
			$resid_size[$j]++;
		}
	}
	if (grep {$_ == $i} @high_layer){
		 printf MAP "%10d\n", $i;
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  H\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
		 printf XYZ "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1]; 
	} elsif ((grep {$_ == $i} @medium_layer) && (grep {$_ == $i} @mlinks)){
                 my $k = 0;
                 while ($k <= $#mlinks){
 	               if ($i == $mlinks[$k]){
	                       last;
                       }
                       $k++;
                 }
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  M H  %10d\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1], $hlinks[$k];
		 push(@dist,${$PDB[1]}[$hlinks[$k]-1] - ${$PDB[1]}[$i-1]);
		 push(@dist,${$PDB[2]}[$hlinks[$k]-1] - ${$PDB[2]}[$i-1]);
                 push(@dist,${$PDB[3]}[$hlinks[$k]-1] - ${$PDB[3]}[$i-1]);
                 $r = sqrt($dist[0]**2 + $dist[1]**2 + $dist[2]**2);
		 if (grep {$_ eq ${$PDB[4]}[$hlinks[$k]-1]} (keys %$default)){
		 	push(@H_pos,${$PDB[1]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[0]);
                 	push(@H_pos,${$PDB[2]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[1]);
                 	push(@H_pos,${$PDB[3]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[2]);
		 	@dist = ();
	 	 } else {
			print "One of the boundery atoms is not C,N or O. Currently, this situation is not implemented.\n";
			exit;
		 }
	} elsif (grep {$_ == $i} @medium_layer){
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  M\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] ne $solvent_name) && (grep {$_ == $i} @llinks)){
                 my $k = 0;
                 while ($k <= $#llinks){
                       if ($i == $llinks[$k]){
                               last;
                       }
                       $k++;
                 }
                 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L H  %10d\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1], $hlinks[$k];
                 push(@dist,${$PDB[1]}[$hlinks[$k]-1] - ${$PDB[1]}[$i-1]);
                 push(@dist,${$PDB[2]}[$hlinks[$k]-1] - ${$PDB[2]}[$i-1]);
                 push(@dist,${$PDB[3]}[$hlinks[$k]-1] - ${$PDB[3]}[$i-1]);
                 $r = sqrt($dist[0]**2 + $dist[1]**2 + $dist[2]**2);
                 if (grep {$_ eq ${$PDB[4]}[$hlinks[$k]-1]} (keys %$default)){
                        push(@H_pos,${$PDB[1]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[0]);
                        push(@H_pos,${$PDB[2]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[1]);
                        push(@H_pos,${$PDB[3]}[$hlinks[$k]-1] - ($default->{${$PDB[4]}[$hlinks[$k]-1]}/$r)*$dist[2]);
                        @dist = ();
                 } else {
                        print "One of the boundery atoms is not C,N or O. Currently, this situation is not implemented.\n";
                        exit;
                 } 
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] ne $solvent_name)){
                 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] eq $solvent_name)){
		 next;
        } else {
                 print "Atom $i${$PDB[4]}[$i-1] was not found in any layer!\n";
		 exit;
        }
}
close(MAP);
for my $i (1..$#{$PDB[0]}+1){
	if ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] eq $solvent_name)){
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	}
}
close(OUT);
open (INP, "<real_layers.xyz") or die;
my $count = 0;
while (<INP>){
	$_ =~ s/^\s*//;
	$_ =~ s/\s*$//;
	@_ = split(/\s+/,$_);
	if ($_[0] ne ${$PDB[4]}[$count]){
		print "Something went wrong during sorting! Please compare atom sequence in real_layer.xyz and the .pdb file!\n";
		exit;
	}
	$count++;
}
print "Sorting was successful!\n";
close(INP);
if ($nonstd_resid ne '0'){
	for my $i (0..$#nonstd_resid){
		open (XYZ2, ">".$nonstd_resid[$i].".xyz") or die;
		printf XYZ2 "%6d\n\n",$resid_size[$i];
		for my $j (0..$#{$PDB[0]}){
			if (${$PDB[0]}[$j] eq $nonstd_resid[$i]){
				printf XYZ2 "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$j], ${$PDB[1]}[$j], ${$PDB[2]}[$j], ${$PDB[3]}[$j];
			}
		}
		close(XYZ2);
	}
}
for my $i (0..(scalar(@hlinks)-1)){
	printf XYZ "%s    %12.6f  %12.6f  %12.6f\n", 'H', $H_pos[3*$i], $H_pos[3*$i+1], $H_pos[3*$i+2];
}
close(XYZ);
}

#####################################

sub get_connection{
my (@zmat,$dist);
foreach my $i (@hlinks){
	foreach my $j (@high_layer){
		if ($j == $i){
			next;
		} else {
			$dist = calculate_dist($i,$j);
			if ($dist <= 2.0){
				push(@{$zmat[0]},$j);
				push(@{$zmat[1]},$dist);
			}
		}
	}
	given (scalar(@{$zmat[0]})){
		when (0) {print "There are no atoms in the high layer within 2.0 Angstrom from the boundery atom $i.\n"; exit}
		when (1) {print_data($i,\@{$zmat[0]},\@{$zmat[1]}); 
			 	foreach my $k (@high_layer){ 
					if ($k == $i || $k == ${$zmat[0]}[0]){
						next;
					} else {
						$dist = calculate_dist(${$zmat[0]}[0],$k);
						if ($dist <= 2.0){
						        push(@{$zmat[0]},$k);
							push(@{$zmat[1]},$dist);
							last;
						}
					}
			 	}
			 }
		when (2) {print_data($i,\@{$zmat[0]},\@{$zmat[1]})}
		when (3) {print_data($i,\@{$zmat[0]},\@{$zmat[1]})}
		when ($_ > 3) {print "There are $_ atoms in the high layer within 2.0 Angstrom from the boundery atom.\n"; exit}
	}
	@zmat = ();
}
}

#####################################

sub print_data{
my $boundary = $_[0];
my @atoms = @{$_[1]};
my @dist = @{$_[2]};
print "Boundery atom $boundary(${$PDB[4]}[$boundary-1]) is connected to\n";
foreach my $i (0..$#atoms){
	print "atom $atoms[$i](${$PDB[4]}[$atoms[$i]-1]) by $dist[$i] Angstrom\n";
}
}

#####################################

sub prepare_nonstd{
my (@gaff_parm,@amber_parm,@parm);
my $nonstd_param = read_input("non-standard residue parameters");
my $nonstd_crg = read_input("non-standard residue charges");
my $parmdat = read_input("force field parameters");
my $predef_resid = read_input("pre-defined non-standard residues");
my @nonstd_param = split(/\s+/,$nonstd_param);
my @nonstd_crg = split(/\s+/,$nonstd_crg);
my @parmdat = split(/\s+/,$parmdat);
foreach my $i (@parmdat){
	my ($dummy1,$dummy2) = split(/\./,$i);
	if ($dummy1 eq "gaff"){
		@gaff_parm = ($dummy1, $dummy2);
	} else {
		@amber_parm =($dummy1, $dummy2);
	}
}
for my $i (0..$#nonstd_resid){
	`babel -ixyz $nonstd_resid[$i].xyz -omol2 $nonstd_resid[$i].mol2`;
	`$AMBERHOME/exe/antechamber -fi mol2 -i $nonstd_resid[$i].mol2 -fo mol2 -o $nonstd_resid[$i]\_out.mol2 -at $nonstd_param[$i] -nc $nonstd_crg[$i] -c bcc`;
	if ($nonstd_param[$i] eq "gaff"){
		@parm = @gaff_parm;
	} else {
		@parm = @amber_parm;
	}
	`$AMBERHOME/exe/parmchk -i $nonstd_resid[$i]\_out.mol2 -f mol2 -p $AMBERHOME/dat/leap/parm/$parm[0].$parm[1] -o $nonstd_resid[$i]\_$parm[0].parm`;
	open (IN, ">leap_".$nonstd_resid[$i].".inp") or die;
	for my $j (@ff){
		print IN "source $j\n";
	}
	print IN "loadamberparams $nonstd_resid[$i]\_$parm[0].parm\n$nonstd_resid[$i] = loadmol2 $nonstd_resid[$i]\_out.mol2\n";
	print IN "check $nonstd_resid[$i]\nsaveoff $nonstd_resid[$i] $nonstd_resid[$i].lib\nquit\n";
	close(IN);
	`tleap -f leap_$nonstd_resid[$i].inp`;
}
if ($predef_resid ne '0'){
	@predef_resid = split(/\s+/,$predef_resid);
	for my $i (0..($#predef_resid-2)/3){
		open (IN, ">leap_".$predef_resid[3*$i].".inp") or die;
	        for my $j (@ff){
			print IN "source $j\n";
 	        }
	        print IN "loadamberparams $predef_resid[3*$i+2]\n";
		$_ = (split(/\./,$predef_resid[3*$i+1]))[1];
		given ($_){
			when ("in"){print IN "loadamberprep $predef_resid[3*$i+1]\n"}
			when ("mol2"){print IN "$predef_resid[3*$i] = loadmol2 $predef_resid[3*$i+1]\n"}
			when ("pdb"){print IN "$predef_resid[3*$i] = loadpdb $predef_resid[3*$i+1]\n"}
		}
	        print IN "check $predef_resid[3*$i]\nsaveoff $predef_resid[3*$i] $predef_resid[3*$i].lib\nsavemol2 $predef_resid[3*$i] $predef_resid[3*$i]\_out.mol2 1\nquit\n";
	        close(IN);
	        `tleap -f leap_$predef_resid[3*$i].inp`;	
	}
}
}

#####################################

sub prepare_model{
my ($correct_type,$correct_charge);
for my $i (0..$#{$PDB[0]}){
	if (grep {$_ eq ${$PDB[0]}[$i]} (@nonstd_resid,@predef_resid,@std_resid)){
	} else {
		push(@std_resid,${$PDB[0]}[$i]);
		open (OUT, ">leap_".${$PDB[0]}[$i].".inp") or die;
		for my $j (@ff){
                	print OUT "source $j\n";
        	}
	        print OUT "saveoff ${$PDB[0]}[$i] ${$PDB[0]}[$i].lib\nsavemol2 ${$PDB[0]}[$i] ${$PDB[0]}[$i]_out.mol2 1\nquit\n";
	        close(OUT);
	        `tleap -f leap_${$PDB[0]}[$i].inp`;
	}
}
`babel -ixyz model-H.xyz -omol2 model-H.mol2`;
open (INP, "<model-H.mol2") or die;
open (OUT, ">tmp") or die;
while(<INP>){
	$_ =~ s/^\s*//;
        @_ = split(/\s+/,$_);
	if ($_[0] && $_[0] =~/\d+/ && $_[0] <= scalar(@high_layer) && $_[5] && $_[5] =~ /\w+/){
		my $model_num = $_[0];
		#print "Number of atom $_[1] in the model-H: $model_num\n";
		my $real_num = $high_layer[$_[0]-1];
		#print "Number of atom $_[1] in real_layers: $real_num\n";
		my $curr_resid = ${$PDB[0]}[$real_num-1];
		push(@model_resid,$curr_resid);
		my $res_num = 0;
		for my $i (0..$real_num-1){
			if (${$PDB[0]}[$i] eq $curr_resid){
				$res_num++;
			}
		}
		#print "Number of atom $_[1] in residue $curr_resid: $res_num\n";
		open (IN, "<".$curr_resid."_out.mol2") or die;
		while (my $tmp = <IN>){
			$tmp =~ s/^\s*//;
		        my @tmp = split(/\s+/,$tmp);
			if ($tmp[0] && $tmp[0] !~ /^\d+$/ ){
				next;
			}
			if ($tmp[0] && $tmp[0] =~ /^\d+$/ && $tmp[1] =~ /^\d+$/ && $tmp[2] =~ /^\d+$/ && $tmp[3] =~ /^\d+$/ && $tmp[4] =~ /^\d+$/){
				print "$tmp\n";
				$res_num = recursive_substraction($res_num,\@tmp);
			#print "Recursive Res_num $res_num\n";
			}
			if ($tmp[0] && $tmp[0] == $res_num && $tmp[5] && $tmp[5] =~ /^\w+/){
				$correct_type = $tmp[5];
				$correct_charge = $tmp[8];
				#print "Correct_type $correct_type\n Correct_charge $correct_charge\n\n";
				last;
			}
		}
		close(IN);
		$_ = join('  ',$_[0],$_[1],$_[2],$_[3],$_[4],$correct_type,$_[6],$_[7],$correct_charge,"\n");
	}
	if ($_[0] && $_[0] =~/\d+/ && $_[0] > scalar(@high_layer) && (split(//,$_[1]))[0] eq 'H'){
		$_[5] = 'AL';
		$_ = join('  ',$_[0],$_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],$_[8],"\n");
	} 
	print OUT $_;
}
close(INP);
close(OUT);
`mv tmp model-H_out.mol2`;
my %dummy = ();
@_ = grep { ! $dummy{ $_ }++ } @model_resid;
@model_resid = @_;
open (IN, ">leap_model.inp") or die;
for my $j (@ff){
	print IN "source $j\n";
}
for my $i (@model_resid){
	print IN "loadoff $i.lib\n";
	if (grep {$_ eq $i} @predef_resid){
		for my $j (0..($#predef_resid-2)/3){
			if ($i eq $predef_resid[3*$j]){
				print IN "loadamberparams $predef_resid[3*$j+2]\n";
			}
		}
	} elsif (grep {$_ eq $i} @std_resid){
	} else {
		$_ = `echo -n $i\*.parm`;
		if ($_ eq $i."*.parm"){
		      	print "There is no .parm file for residue $i\n";
		        exit;	
		} elsif (scalar(split(/\s+/,$_)) > 1){
			print "There is more than one .parm file for residue $i\n";
			exit;
		} else {
			print IN "loadamberparams $_\n";
		}
	}
}
print IN "model = loadmol2 model-H_out.mol2\n";
print IN "check model\nsaveamberparm model model-H.top model-H.crd\nsavepdb model model-H.pdb\nquit\n";
close(IN);
`rm leap.log; tleap -f leap_model.inp`;
my $link_atoms = read_input("links");
if ($link_atoms eq '0'){
        return;
} else {
	open (OUT, ">model-H.parm") or die;
	print OUT "remark goes here\nMASS\n";
	@_ = ();
	print OUT "AL 0.000         0.000               ATTN, need revision\n\n";
	open (IN, "<leap.log") or die;
	print OUT "BOND\n";
	open (IN, "<leap.log") or die;
	while (<IN>){
		if (/^\s*Could not find bond parameter for:\s+(\w.*)/){
			if (grep {$_ eq $1} @_){
			} else {
				print OUT "$1    0.00   0.000       ATTN, need revision\n";
				push(@_,$1);
			}
		}
	}
	close(IN);
	open (IN, "<leap.log") or die;
	print OUT "\nANGLE\n";
	while (<IN>){
	        if (/^\s*Could not find angle parameter:\s+(\w.*)/){
			if (grep {$_ eq $1} @_){
			} else {
	                	print OUT "$1   0.00   0.000       ATTN, need revision\n";
				push(@_,$1);
			}
	        }
	}
	close(IN);
	print OUT "\nDIHE\n\n";
	print OUT "IMPROPER\n\n";
	print OUT "NONBON\n  AL          0.0000  0.0000             ATTN, need revision\n"; 
	close(OUT);
	open (IN, "<leap_model.inp") or die;
	open (LEAP, ">tmp") or die;
	while (<IN>){
		if (/check model/){
			print LEAP "loadamberparams model-H.parm\n";
			print LEAP "check model\n";
		} else {
			print LEAP "$_";
		}
	}
	close(LEAP);
	`mv tmp leap_model.inp; tleap -f leap_model.inp`;
}
}

######################################

sub recursive_substraction{

my $res_num = $_[0];	
my $tmp = $_[1];
my @tmp = @$tmp;
if ($res_num > $tmp[0]){
	$res_num = $res_num-$tmp[0];
	recursive_substraction($res_num,\@tmp);
} elsif ($res_num <= $tmp[0]){
	return $res_num;
}
}

######################################

sub prepare_real{
my $count = 0;
for my $i (@hlinks){
	if (grep {$_ eq ${$PDB[0]}[$i-1]} @UNIQUE_RES){
		for my $j (0..$#RES){
			if (grep {$_ eq $i} @{$RES[$j]}){
				push (@REASSIGNED_RES,@{$RES[$j]})
			}
		}
        } else {
		push(@UNIQUE_RES,${$PDB[0]}[$i-1]);
                $new_resid->{${$PDB[0]}[$i-1]} = $alphabet[$count].$alphabet[$count].$alphabet[$count];
                $count++;
		for my $j (0..$#RES){
	                if (grep {$_ eq $i} @{$RES[$j]}){
		                push (@REASSIGNED_RES,@{$RES[$j]})
		        }       
		}
        }
}
my %dummy = ();
my @REASSIGNED_RES = grep { ! $dummy{ $_ }++ } @REASSIGNED_RES;
print "ATOMS IN RESIDUES TO BE REASSIGNED: @REASSIGNED_RES\n";
print "UNIQUE RESIDUES FOR CHARGE REASSIGNMENT: @UNIQUE_RES\n";
for my $i (keys %$new_resid){
        print "Old residue name: $i, New residue name: $new_resid->{$i}\n";
}
foreach my $i (@UNIQUE_RES){
        reassign_charges($i);
}
my $coord = read_input("coordinate file");
open (INP, "<".$coord) or die;
open (OUT, ">real.pdb") or die;
my @tmp;
while (my $tmp = <INP>){
	if ($tmp =~ /^\s*ATOM\s+(\d+)/){
	       if ((grep {$_ eq $1} @REASSIGNED_RES)){
			$tmp =~ s/^\s*/$tmp/;
			$tmp =~ s/\s*$/$tmp/;
			@_ = split(/\s+/,$tmp);
			@tmp = split(//,$_[2]);
			if (scalar(@tmp) == 4){
				$_[2] = $tmp[0].$tmp[1].$tmp[2].$tmp[3];
			} elsif ((scalar(@tmp) == 3 && $tmp[2] eq '*') || (scalar(@tmp) == 3 && $tmp[0] =~ /[A-Z]/)){
				$_[2] = " ".$tmp[0].$tmp[1].$tmp[2];
			} elsif (scalar(@tmp) == 3 && $tmp[0] =~ /[0-9]/){
				$_[2] = $tmp[0].$tmp[1].$tmp[2]." ";
			} elsif (scalar(@tmp) == 2){
				$_[2] = " ".$tmp[0].$tmp[1]." ";
			} elsif (scalar(@tmp) == 1){
				$_[2] = " ".$tmp[0]."  ";
			}
			printf OUT "ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n",$_[1],$_[2],$new_resid->{$_[3]},$_[4],$_[5],$_[6],$_[7];
	       } elsif ((grep {$_ == $1} @high_layer,@medium_layer) || ((grep {$_ == $1} @low_layer) && (${$PDB[0]}[$1-1] ne $solvent_name))) {
			print OUT "$tmp";
	       } else {
			next;
	       }
	}
}
close(INP);
open (INP, "<".$coord) or die;
while (my $tmp = <INP>){
	if ($tmp =~ /^\s*ATOM\s+(\d+)/){
		if ((grep {$_ == $1} @low_layer) && (${$PDB[0]}[$1-1] eq $solvent_name)){
			print OUT "$tmp";
		}
	}
}
close(INP);
close(OUT);
open (IN, ">leap_real.inp") or die;
for my $j (@ff){
        print IN "source $j\n";
}
foreach my $i (@nonstd_resid){
	print IN "loadoff $i.lib\n";
	$_ = `echo -n $i\*.parm`;
	if ($_ eq $i."*.parm"){
	        print "There is no .parm file for residue $i\n";
		exit;
	} elsif (scalar(split(/\s+/,$_)) > 1){
		print "There is more than one .parm file for residue $i\n";
		exit;
	} else {
		print IN "loadamberparams $_\n";
	}
}
foreach my $i (@std_resid){
	print IN "loadoff $i.lib\n";
}
foreach my $i (0..($#predef_resid-2)/3){
	print IN "loadoff $predef_resid[3*$i].lib\n";
	print IN "loadamberparams $predef_resid[3*$i+2]\n";
}
foreach my $i (@UNIQUE_RES){
        print IN "loadoff $new_resid->{$i}.lib\n";
}
print IN "real = loadpdb real.pdb\n";
if (read_input("additional bonds")){
	$_ = read_input("additional bonds");
	my @addbonds = split(/\s+/,$_);
	for my $i (@addbonds){
		@_ = split(/-/,$i);
		print IN "bond real.$_[0] real.$_[1]\n";
	}
}
print IN "saveamberparm real real.top real.cord\n";
print IN "savepdb real real_libs.pdb\nquit\n";
close(IN);
`tleap -f leap_real.inp`;
}

######################################

sub reassign_charges{
my $resid = $_[0];
my @resid_high_layer = ();
my $resid_charge;
open (IN, ">leap_crg.inp") or die;
for my $j (@ff){
        print IN "source $j\n";
}
print IN "loadoff $resid.lib\n";
print IN "charge $resid\nquit\n";
close(IN);
`if [ -f leap.log ]; then rm leap.log; fi; tleap -f leap_crg.inp`;
open (LEAP, "<leap.log") or die;
while (<LEAP>){
	if (/^\s*Total unperturbed charge:\s+(-{0,1}\d+\.\d+)/){
		$resid_charge = $1;
	}
}
print "Re-assigning $resid\n";
foreach my $i (@high_layer){
	if (${$PDB[0]}[$i-1] eq $resid){
		my $count = 0;
		for my $j (0..$i-2){
			if (${$PDB[0]}[$j] eq $resid){
				$count++;
			}
		}
		push(@resid_high_layer,$count);
	}
}
print "Resid_high_layer: @resid_high_layer\n";
print "Collecting residual charge.\n";
open (LIB, "<".$resid.".lib") or die;
my $model_charge = 0;
my $QM_model_charge = read_input("high layer charge");
my $charge = 0;
my $abs_charge = 0;
while (<LIB>){
	if (/!entry\.$resid\.unit\.atoms table/){
		my $count = 0;
		while (my $tmp = <LIB>){
			if ($tmp =~ /!entry\.$resid\.unit\.atomspertinfo/){
                                last;
                        }
                        @_ = split(/\s+/,$tmp);
			if (grep {$_ == $count} @resid_high_layer){
				$model_charge += $_[8];
			} elsif ($_[1] eq '"P"') {
			        $charge += $_[8];
			} else {
				$charge += $_[8];
				$abs_charge += abs($_[8]);
			}
			$count++;
		}
		last;
	}
}
$model_charge = $model_charge - $QM_model_charge;
printf "Total charge the residue should have: %9.6f\n", $resid_charge;
printf "Total residual charge to redistribute: %9.6f\n", $model_charge;
printf "Total charge the residue has before redistribution of the residual charge: %9.6f\n", $charge;
$charge = 0;
close(LIB);
#$abs_charge -= abs($resid_charge);
open (INP, "<".$resid.".lib") or die;
open (LIB, ">".$new_resid->{$resid}.".lib") or die;
while (<INP>){
	if (/!entry\.$resid\.unit\.atoms table/){
		my $count = 0;
		print LIB "!entry.$new_resid->{$resid}.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n";
		while (my $tmp = <INP>){
			if ($tmp =~ /!entry\.$resid\.unit\.atomspertinfo/){
				print LIB "!entry.$new_resid->{$resid}.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n";
                                last;
                        }
			@_ = split(/\s+/,$tmp);
			if (grep {$_ == $count} @resid_high_layer){
                                print LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] 0.000000\n";
                        } elsif ($_[1] eq '"P"') {
				printf LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] $_[8]\n";
				$charge += $_[8];
			} else {
				printf LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] %9.6f\n", $_[8]+abs($_[8])*$model_charge/$abs_charge;
				$charge += $_[8]+abs($_[8])*$model_charge/$abs_charge;
                        }
                        $count++;
		}
	} elsif (/(.*)$resid(.*)/){
		print LIB "$1$new_resid->{$resid}$2\n";
	} else {
		print LIB $_;
	}
}
close(INP);
close(LIB);
printf "Total charge the residue has after redistribution of the residual charge: %9.6f\n", $charge;
}

########################################

sub check_atom_types{
my (@model_types,$model_types,@real_types,@pos);
my ($line,$field,$sed,@tmp);
foreach my $i (@high_layer){
	$line = floor(($i-1)/20);
        $field = ($i-1)%20;
	open (INP, "<real.top") or die;
	while (<INP>){
	        if (/^\s*%FLAG AMBER_ATOM_TYPE/){
			$_ = <INP>;
			for my $j (0..$line){
				$_ = <INP>;
			}
			chomp;
			$_ =~ s/^\s*//;
			$_ =~ s/\s*$//;
			@_ = split(/\s+/,$_);
			if (scalar(split(//,$_[$field])) == 1){
				$_[$field] = $_[$field]." ";
			} elsif ((split(//,$_[$field]))[1] eq '*'){
				@tmp = split(//,$_[$field]);
	                        $_[$field] = $tmp[0]."\\".$tmp[1];
			}
			push(@real_types,$_[$field]);
			last;
		}
	}
	close(INP);
}
print "Atom types in real.top   : @real_types\n";
$model_types = types_from_model();
print "Atom types in model-H.top: @$model_types\n";
($sed,@pos) = compare_real_model(\@real_types,$model_types);
if ($sed eq 'yes'){
print "Atoms types of model-H_out.mol2 will be corrected.\n";
open (MOL, "<model-H_out.mol2") or die;
open (MOL2, ">tmp") or die;
while (<MOL>){
	if (/\@<TRIPOS>ATOM/){
		print MOL2 "@<TRIPOS>ATOM\n";
		for my $i (0..$#$model_types){
			my $_ = <MOL>;
			if (grep {$_ == $i} @pos){
				$_ =~ s/^\s*//;
				$_ =~ s/\s*$//;
				@_ = split(/\s+/,$_);
				printf MOL2 "$_[0] $_[1]        $_[2]  $_[3]  $_[4]  $real_types[$i]         $_[6]  $_[7]  $_[8]\n";
			} else {
				print MOL2 $_;
			}
		}
	} else {
		print MOL2 $_;
	}
}
close (MOL2);
close(MOL);
open (IN, ">leap_model.inp") or die;
for my $j (@ff){
        print IN "source $j\n";
}
for my $i (@model_resid){
        print IN "loadoff $i.lib\n";
}
print IN "\nmodel = loadmol2 model-H_out.mol2\n";
print IN "check model\nsaveamberparm model model-H.top model-H.crd\nsavepdb model model-H.pdb\nquit\n";
close(IN);
`tleap -f leap_model.inp`;
$model_types = types_from_model();
($sed,@pos) = compare_real_model(\@real_types,$model_types);
}
}	

#####################################

sub compare_real_model{
my $sed = 'no';
my @real_types = @{$_[0]};
my @model_types = @{$_[1]};
my @pos = ();
for my $i (0..$#model_types-scalar(@hlinks)){
	if ($real_types[$i] ne $model_types[$i]){
		printf "real.top and model-H.top differ in atom type for atom $high_layer[$i](real)/%d(high): $real_types[$i] vs. $model_types[$i]\n", $i+1;
		push(@pos,$i);
		$sed = 'yes';
	}
}
if ($sed eq 'no'){
	print "Atom types of real.top and model-H.top identical.\n";
}
return ($sed,@pos);
}

#####################################

sub types_from_model{
my (@model_types,@tmp);
open (INP, "<model-H.top") or die;
while (<INP>){
	if (/^\s*%FLAG AMBER_ATOM_TYPE/){
		$_ = <INP>;
		for my $i (0..floor(scalar(@high_layer)/20)){
			$_ = <INP>;
			chomp;
			$_ =~ s/^\s*//;
			$_ =~ s/\s*$//;
			@_ = split(/\s+/,$_);
			for my $j (0..$#_){
				if (scalar(split(//,$_[$j])) == 1){
                                	$_[$j] = $_[$j]." ";
                        	} elsif ((split(//,$_[$j]))[1] eq '*'){
                                	@tmp = split(//,$_[$j]);
                                	$_[$j] = $tmp[0]."\\".$tmp[1];
                        	}
			}
			push(@model_types,@_);
		}
		last;
	}
}
close(INP);
return(\@model_types);
}

######################################

sub keep_files{
	my $coord = read_input("coordinate file");
        system("if [ ! -d \$PWD/input_files ]; then mkdir input_files; fi; cp -f real_layers.xyz real.top model-H.top cobram.parm $coord pdb_2_HL.map input_files/.");
}




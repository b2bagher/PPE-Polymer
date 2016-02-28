#!/usr/bin/perl

# UNKs is how atom is marked in a PDB file (HETATM or ATOM)
$record = "HETATM";

# labeling of the atoms
# old atom index, new atom index, new atom name, residue name, residue number  
@map = (
[1, 1, HA5, RNG,  1],

[2,2 , CA1, RNG, 1],
[3,3 , CA2, RNG, 1],
[4,4 , CA3, RNG, 1],
[5,5 , CA4, RNG, 1],
[21,6, HA1, RNG, 1],
[6,7 , HA2, RNG, 1],
[22,8, HA3, RNG, 1],
[7,9, HA4, RNG, 1],
[8,10 , CA5, RNG,1],
[9,11 , CAL, RNG,1],

[10,12, CL1, ETH, 2],
[11,13, CL2, ETH, 2],

[18,14, CAL, RNG, 3],
[12,15, CA2, RNG, 3],
[13,16, CA1, RNG, 3],
[14,17, CA4, RNG, 3],
[15,18, CA3, RNG, 3],
[23,19, HA2, RNG, 3],
[16,20, HA1, RNG, 3],
[24,21, HA4, RNG,3],
[17,22, HA3, RNG, 3],
[19,23,CA5, RNG, 3],
[20,24, HA5, RNG, 3],

);

print "MAPPING USED:\n";

for $number (0 .. $#map ) {
    for $label ( 0 .. $#{ $map[$number] } ) {
     print "$map[$number][$label] ";
    }
    print "\n";
}

my @pdb = read_pdb("trajectory.pdb");
$n = 1; $conf = 1; $a = 1;
my $pdb_out = "trajectory_reordered.pdb";


print "Opening file $pdb_out\n";
open(OUT, ">$pdb_out") || die "cannot open $pdb_out for reading: $!  ";

# loop over the possible configurations 
for (my $i=0; $i< @pdb; $i++ ) {

  my $record_name    =  trim( substr $pdb[$i], 0, 6 );
  if ( $record_name eq $record ) {
    my $atom_number    =  trim( substr $pdb[$i], 6, 5);
    my $atom_name      =  trim( substr $pdb[$i], 12, 4);
    my $residue_name   =  trim( substr $pdb[$i], 17, 3);
    my $residue_number =  trim( substr $pdb[$i], 22, 4);
    my $x              =  trim( substr $pdb[$i], 30, 8);
    my $y              =  trim( substr $pdb[$i], 38, 8);
    my $z              =  trim( substr $pdb[$i], 46, 8);
    my $element        =  trim( substr $pdb[$i], 76, 2);
    my $charge         =  substr $pdb[$i], 78, 2;
    chomp $charge;     
    #print "$record_name $atom_number $atom_name $residue_name $residue_number $x $y $z $element $charge\n";
    #$old_line = chomp $pdb[$i]; print "OLD: $old_line\n"; 



    
    for ($a=1; $a<=$#map+1; $a++ ){		
	
    if ($map[$a-1][0] eq $atom_number) {
      $new_number = $map[$a-1][1];

      # find the corresponding label
      for $number (0 .. $#map) {
	if ($map[$number][0] eq $atom_number) {
	  $new_atom_name = $map[$number][2];
	  $new_residue_name = $map[$number][3];
	  $new_residue_number = $map[$number][4];
          #print "$number $atom_name $new_atom_name\n";	
	}
      }

      $new_line = sprintf '%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f                     %3s%2s', 
      $record_name, $new_number, $new_atom_name, $new_residue_name, $new_residue_number, $x, $y, $z, $element, $charge;
      $final_line[$new_number]=$new_line; 
      #print "FIN: $final_line[$new_number], $new_number\n";
      #$a++;
    } 
}
    $n++;
  } else {
    #print "$a,$n\n";
    $n = $n-1; $a = $a-1;
    #print "$pdb[$i], configuration: $conf, atoms processed: $n \n";
    if ($n > 1) {
      for $number(1 .. $#map+1) {
	print OUT "$final_line[$number]\n";
	#print "$final_line[$number]\n";
      }
      $conf++;
    };
    if ($a = 1) {
      print OUT "$pdb[$i]";
      #print "$pdb[$i],$a,$n\n";
    }
    $n = 1; $a = 1;
  }
}

close(OUT);


sub read_pdb {
  my( $pdb_file ) = @_;
  print "Opening file $pdb_file for reading.\n";
  open(IN, "$pdb_file") || die "cannot open $pdb_file for reading: $!  ";
  @pdb = <IN>;
  my $natoms = @pdb;
  close(IN);
  print "Number of atoms in a molecule: $natoms\n";
  return @pdb;
}

sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}


# pdb file specification (field - column)
# 1 -  6        Record name     "ATOM  "                                            
# 7 - 11        Integer         Atom serial number.                   
#13 - 16        Atom            Atom name.                            
#17             Character       Alternate location indicator.         
#18 - 20        Residue name    Residue name.                         
#22             Character       Chain identifier.                     
#23 - 26        Integer         Residue sequence number.              
#27             AChar           Code for insertion of residues.       
#31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
#39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
#47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
#55 - 60        Real(6.2)       Occupancy.                            
#61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
#73 - 76        LString(4)      Segment identifier, left-justified.   
#77 - 78        LString(2)      Element symbol, right-justified.      
#79 - 80        LString(2)      Charge on the atom.       


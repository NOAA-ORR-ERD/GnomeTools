import re


WHITESPACE_PATTERN = re.compile( "\s+" )


def load_mass( mass_filename ):
    mass_file = open( mass_filename, "rU" )

    mass_file.readline() # Skip the header lines.
    mass_file.readline()

    data = WHITESPACE_PATTERN.split( mass_file.readline().strip() )
    mass_file.close()

    volume = float( data[ 3 ].strip() )
    density = float( data[ 6 ].strip() )
    mass = volume * density

    return mass

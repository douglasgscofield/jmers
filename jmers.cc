/*  This file is part of jmers.

    jmers is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jmers is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jmers.  If not, see <http://www.gnu.org/licenses/>.
*/


/*
    This file is heavily based on the Jellyfish "query_per_sequence.cc" example.
    If it breaks the GPL to use quite a bit of it as a template, please let me
    know!

    / Martin Norling
*/

//#include "Seq.h"
#include "Input.h"
#include "FosmidEndFragment.h"
#include "KmerBoundary.h"
#include "jmerfish.hpp"
#include <zlib.h>
#include <iostream>

#include <vector>

int main(int argc, char *argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "USAGE: " << argv[0] << " jellyfish_db.jf seq_file [seq_file [...]]" << std::endl;
        return 0;
    }
    
    // Load jellyfish database
    jmers::JellyfishDatabase jf_db( argv[1] );
    jmers::KmerBoundarySimple kmer_boundary( &jf_db );
    
    // fastq-parsing variables
    jmers::Seq s;
    
    for (int i = 2; i < argc; i++)
    {
        try
        {
            jmers::Input seq_file = jmers::Input(argv[i]);
            while ( seq_file.read(s) )
            {
                jmers::FosmidEndFragment fosmid_end_fragment = kmer_boundary.detect_boundary(s);
                // fosmid_end_fragment.infer_fragment_structure( db );
                // fosmid_end_fragment.split_fragment();
                // fosmid_end_fragment.write_pair_fasta();
                fosmid_end_fragment.dump();
                
                break;
            }
        }
        catch (std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
    
    return 0;
}

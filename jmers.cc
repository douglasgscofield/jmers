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

// from https://code.google.com/archive/p/simpleopt/
#include "SimpleOpt.h"
#include "SimpleGlob.h"

/*
    CURRENT TODO:
 
  o figure out if a canonical database works
    o should at least have different workings for canonical/non-canonical (prob. prefer/enforce non-canonical)
  o is a bloom-counter a smaller database?
  o calculate some probabilities!
        o do I need poisson?
        o is GC-content biased * qual ok?
  o do DFS on paths, and basically follow forever if needed
  o evaluate ranges between kmer-stretches instead of at edges to try to "meet" kmer-walking
  o make something simple to visualize paths in OpenGL
  o use SimpleOpt.h for arg parsing (find at https://github.com/douglasgscofield/samla)
  o 
 
 
 
 
 
 
 
 
 
 */

// SimpleOpt command line argument handling
enum { OPT_HELP, OPT_ERROR, OPT_KMER, OPT_LIMIT };
CSimpleOpt::SOption g_rgOptions[] = {
    { OPT_ERROR, "-e",     SO_NONE    }, // "-e"        default: true
    { OPT_KMER,  "-k",     SO_NONE    }, // "-k"        default: false
    { OPT_LIMIT, "-l",     SO_REQ_SEP }, // "-l"        read limiter
    { OPT_HELP,  "-?",     SO_NONE    }, // "-?"
    { OPT_HELP,  "-h",     SO_NONE    }, // "-h"
    { OPT_HELP,  "--help", SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};

void Help() {
    printf( "USAGE: ./jmers [-e] [-k] [-?|-h|--help] <Jellyfish database> <fastq file>\n" );
}


int main(int argc, char *argv[])
{
    // Process arguments with SimpleOpts
    bool error_correct = true;
    bool print_kmer_content = false;
    int read_limit = -1;
    std::string database;
    CSimpleOpt args(argc, argv, g_rgOptions);

    while (args.Next())
    {
        if (args.LastError() == SO_SUCCESS) {
            if (args.OptionId() == OPT_HELP) {
                Help();
                return 0;
            }
            else if (args.OptionId() == OPT_ERROR)
                error_correct = false;
            else if (args.OptionId() == OPT_KMER)
                print_kmer_content = true;
            else if (args.OptionId() == OPT_LIMIT)
                read_limit = atoi(args.OptionArg());
            else
                printf("Option, ID: %d, Text: '%s', Argument: '%s'\n",
                    args.OptionId(), args.OptionText(),
                    args.OptionArg() ? args.OptionArg() : "");
        }
        else {
            printf("Invalid argument: %s\n", args.OptionText());
            return 1;
        }
    }
    
    CSimpleGlob glob(SG_GLOB_NODOT|SG_GLOB_NOCHECK);
    if (SG_SUCCESS != glob.Add(args.FileCount(), args.Files())) {
        printf("Error while globbing files\n");
        return 1;
    }
    
    if ( args.FileCount() < 2 )
    {
        Help();
        return 0;
    }
    
    database = args.Files()[0];
    
    // Load jellyfish database
    jmers::JellyfishDatabase jf_db( database );
    
    //std::cout << "Loaded Jellyfish database, with kmer " << jf_db.kmer << std::endl;
    
    // fastq-parsing variables
    jmers::Seq s;
    
    int c = 0;
    
    for (int i = 1; i < args.FileCount(); i++)
    {
        try
        {
            jmers::Input seq_file = jmers::Input( args.Files()[i] );
            while ( seq_file.read(s) )
            {
                jmers::KmerBoundarySimple kmer_boundary( &jf_db );
                jmers::FosmidEndFragment fosmid_end_fragment = kmer_boundary.detect_boundary(s, error_correct, print_kmer_content);
                //fosmid_end_fragment.split_fragment();
                //fosmid_end_fragment.write_pair_fasta();
                //fosmid_end_fragment.dump();
                
                if ( read_limit > 0 and ++c > read_limit )
                    break;
            }
        }
        catch (std::exception &e)
        {
            std::cerr << "Something went wrong!" << std::endl;
            std::cerr << e.what() << std::endl;
        }
    }
    
    return 0;
}

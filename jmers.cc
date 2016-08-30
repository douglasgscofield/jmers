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

// SimpleOpt command line argument handling
enum { OPT_HELP, OPT_ERROR };
CSimpleOpt::SOption g_rgOptions[] = {
    { OPT_ERROR, "-e",     SO_NONE    }, // "-e"
    { OPT_HELP,  "-?",     SO_NONE    }, // "-?"
    { OPT_HELP,  "-h",     SO_NONE    }, // "-?"
    { OPT_HELP,  "--help", SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};

void Help() {
    printf( "USAGE: ./jmers [-e] [-?|-h|--help] <Jellyfish database> <fastq file>\n" );
}


int main(int argc, char *argv[])
{
    // Process arguments with SimpleOpts
    bool error_correct = true;
    std::string database;
    CSimpleOpt args(argc, argv, g_rgOptions);

    while (args.Next()) {
        if (args.LastError() == SO_SUCCESS) {
            if (args.OptionId() == OPT_HELP) {
                Help();
                return 0;
            }
            else if (args.OptionId() == OPT_ERROR)
                error_correct = false;
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
    jmers::KmerBoundarySimple kmer_boundary( &jf_db );
    
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
                jmers::FosmidEndFragment fosmid_end_fragment = kmer_boundary.detect_boundary(s);
                //fosmid_end_fragment.split_fragment();
                //fosmid_end_fragment.write_pair_fasta();
                //fosmid_end_fragment.dump();
                
                if ( ++c >= 100)
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

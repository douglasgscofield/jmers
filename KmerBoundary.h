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

#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include "jmerfish.hpp"

namespace jmers {

    /*
    TODO: use KmerDB class instead of jellyfish db once things are working.
    
class KmerDB {
    bool        is_empty;
    std::string description;
  public:
    KmerDB(...)  // initialise with jellyfish db?
    { }
    KmerDB()
        : is_empty(true)
    { }
};  // class KmerDB

KmerDB empty_db = KmerDB();
*/
class KmerBoundary {
    
  public:
    KmerBoundary() {}
    ~KmerBoundary() {}

};   // class KmerBoundary

class KmerBoundarySimple : KmerBoundary {
    JellyfishDatabase *A;
    enum sense_t { sense_none = 0, sense_AB, sense_BA };
    
    std::string extend_path(Seq& seq, int max_depth = 4, int i = 0, bool expand = false)
    {
        /* TODO: This function should be improved to try changes up until max_depth.
           Right now it just corrects single errors...
        */
        if ( i == seq.sequence.length() )
            return std::string();
        
        if ( A->query( seq.sequence.substr(i, A->kmer) ) == 0 )
        {
            if (expand)
            {
                for ( auto c : "ACTG" )
                {
                    if ( c == seq.sequence[i+A->kmer] )
                        continue;
                    if ( A->query(seq.sequence.substr(i, A->kmer-1) + c) )
                        seq.sequence[i+ A->kmer -1] = c;
                }
            }
        }
        else
        {
            /* Once we're out of the first zero region, we can start trying to 
               error correct things!
             */
            expand = true;
        }
        return seq.sequence[i] + extend_path(seq, max_depth, i+1, expand);
    }
    
    Seq find_path(Seq& seq, int max_depth = 4, int i = 0, bool expand=false )
    {
        Seq out;
        out.comment     = seq.comment + " error corrected";
        out.has_quality = seq.has_quality;
        out.name        = seq.name;
        /*
         *  I've made this so that the first base is never changed. The reasoning
         *  is that if it's a zero, it's considered to be part of the "junk" region
         *  and if it's not, it shouldn't be changed...
         */
        out.extend( seq.sequence[i] + extend_path(seq, max_depth, i+1, expand) );
        return out;
    }
    
  public:
    KmerBoundarySimple(JellyfishDatabase *a) { A = a; }
    ~KmerBoundarySimple() {}
    std::string describe() const {
        std::stringstream s;
        s << "A:" << A->describe();
        return s.str();
    }
    jmers::FosmidEndFragment detect_boundary(Seq& seq)
    {
        Seq corrected_seq = find_path(seq);
        seq.dump();
        corrected_seq.dump();
        jmers::FosmidEndFragment fosmid_end(corrected_seq);
        
        int i;
        std::vector<int> kmer_content;
        
        for (int i = 0; i < seq.sequence.length()-A->kmer; i++)
        {
            kmer_content.push_back( A->query( seq.sequence.substr(i, A->kmer) ) );
            std::cout << A->query( seq.sequence.substr(i, A->kmer) ) << " ";
        }
        std::cout << std::endl;
        
        return fosmid_end;
    }
};   // class KmerBoundarySimple



// perhaps a better design would be KmerBoundaryChimeric/KmerBoundarySimilar,
// where either side of the boundary is from the same kmer database, and
// KmerBoundaryDissimilar/KmerBoundaryDistinct, where the two sides are from
// different databases, one of them possibly being no database

// KmerBoundarySimilar would not have a sense
// KmerboundaryDissimilar would
/*
class KmerBoundaryDissimilar {
    KmerDB A;
    KmerDB B;
    enum sense_t { sense_none = 0, sense_AB, sense_BA };
    KmerBoundaryDissimilar(KmerDB& a = empty_db, KmerDB& b = empty_db)
        : A(a), B(b)
    { }
    ~KmerBoundaryDissimilar()
    { }
    std::string describe() const {
        std::stringstream s;
        s << "A:" << A->describe() << ", B:" << B.describe();
        return s.str();
    }
    int64_t detect_boundary(Seq& seq, sense_t sns, int64_t pos = 0,
            int tol1 = 0, double tol2 = 0.0) {
        const KmerDB * local_a;
        const KmerDB * local_b;
        // set up kmer db pointers
        if (sns == sense_AB) { local_a = &A; local_b = &B; }
        if (sns == sense_BA) { local_a = &B; local_b = &A; }
        // iterate stacked A & B queries from the start
        // A: -----------
        // B: -----------

        // for dissimilar, the chances that both are in the database depends on
        // the prior information we are given about the sense, but it is always
        // possible that both are.  if we are near the boundary, then both will
        // not be.  we are sliding from a range of one being found to the other
        // being found.

        // perhaps start with a rough query of A and B kmers along the entire
        // fragment then i would learn if there are A-and-not-B,
        // not-A-and-not-B, and B-and-not-A zones before starting

        // this still needs to cook a bit

    }
};   // class KmerBoundary
*/
}  // namespace jmers


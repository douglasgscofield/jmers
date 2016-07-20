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
#include <algorithm>
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
    std::vector<int> kmer_content;
    std::vector<int> base_support;
    
    Seq error_correct(Seq& seq, int max_corrections_per_kmer = 4, bool exhaustive = false, float w_b = .5, float w_q = .5)
    {
        Seq out;
        out.comment     = seq.comment + " error corrected";
        out.has_quality = seq.has_quality;
        out.name        = seq.name;
        out.sequence    = seq.sequence;
        out.quality     = seq.quality;
        
        int last_hit = A->kmer;
        std::vector<int> positions;
        positions.resize( max_corrections_per_kmer, 0 );
        for (int i = 0; i < out.sequence.length(); i++)
        {
            if (i < this->kmer_content.size()-A->kmer)
            {
                int kmer = A->query( out.sequence.substr(i, A->kmer) );
                if ( kmer == 0 )
                {
                    std::string current_kmer = out.sequence.substr(i, A->kmer);
                    last_hit = last_hit < A->kmer ? last_hit+1 : A->kmer;
                    
                    /* ok - so... should I go through all possible positions?
                        or should I make a list of positions and only check those
                        with the lowest base_support? should I factor in quality?
                        some sort of a*base_support+b*quality weighted value in 
                        a list probably. then loop over the list and check 
                        substitutions for all positions until we find the one
                        with the highest score... or just one with "a" score?
                        exhaustive search of all posibilities, or just check 
                        all variables independent?
                        
                        I'm gonna need to recurse this a bit...
                    
                    */
                    
                    // make a list of all possible positions that could be "wrong"
                    std::map<int, double> position_values;
                    for (int p = 0; p < last_hit; p++)
                    {
                        int pos = i+A->kmer-last_hit+p;
                        if (out.has_quality)
                            position_values[pos] = w_b*base_support[pos] + w_q*((int)out.quality[pos]-33);
                        else
                            position_values[pos] = base_support[pos];
                    }
                    
                    // remove the highest support possible position until no more than
                    // "max_corrections_per_kmer" remain.
                    
                    while (position_values.size() > max_corrections_per_kmer)
                    {
                        auto x = std::max_element(position_values.begin(), position_values.end(),
                            [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                return p1.second < p2.second; });
                        position_values.erase(x);
                    }
                    // for (auto p: position_values)
                    //     std::cout << p.first << " " << p.second << "\t";
                    // std::cout << std::endl;
                    
                    if (exhaustive)
                    {
                        // TODO: implement this as a recursion.
                        1+1;
                    }
                    else
                    {
                        for ( auto p : position_values )
                        {
                            int pos = p.first-i;
                            for ( char c : {'A', 'T', 'C', 'G'} )
                            {
                                if ( c == out.sequence[p.first] )
                                    continue;
                                std::string new_kmer = current_kmer;
                                new_kmer[pos] = c;
                                
                                if ( A->query( new_kmer ) )
                                {
                                    out.sequence[p.first] = c;
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    last_hit = 0;
                }
            }
        }
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
    jmers::FosmidEndFragment detect_boundary(Seq& seq, bool error_correct = true)
    {
        /* Check kmer content of current sequence */
        kmer_content.clear();
        for (int i = 0; i < seq.sequence.length()-A->kmer; i++)
        {
            kmer_content.push_back( A->query( seq.sequence.substr(i, A->kmer) ) );
        }
        base_support.clear();
        base_support.resize( seq.sequence.length(), 0 );
        int o = 0;
        for ( auto k : kmer_content)
        {
            for (int i = 0; i < A->kmer; i++)
            {
                base_support[o+i] += k;
            }
            o++;
        }
        
        /* use the kmer content to error correct the sequence? */
        Seq corrected_seq = error_correct ? this->error_correct(seq) : seq;
        if (error_correct)
        {
            kmer_content.clear();
            /* Update the kmer content */
            for (int i = 0; i < seq.sequence.length()-A->kmer; i++)
            {
                kmer_content.push_back(  A->query( corrected_seq.sequence.substr(i, A->kmer) ) );
            }
            /* Update the base support */
            base_support.clear();
            base_support.resize( seq.sequence.length(), 0 );
            o = 0;
            for ( auto k : kmer_content)
            {
                for (int i = 0; i < A->kmer; i++)
                {
                    base_support[o+i] += k;
                }
                o++;
            }
        }
        
        /* print kmer content and base support */
        std::cout << ">" << seq.name << std::endl;
        for ( auto c : kmer_content )
            std::cout << c << " ";
        std::cout << std::endl;
        for ( auto b : base_support )
            std::cout << b << " ";
        std::cout << std::endl;
        
        /* create fosmid end data structure */
        jmers::FosmidEndFragment fosmid_end(corrected_seq);
        //fosmid_end.infer_fragment_structure(kmer_content, base_support);
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


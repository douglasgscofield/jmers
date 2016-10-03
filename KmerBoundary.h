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
#include <map>
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

struct extension {
    char base = '-';
    int kmer_support = 0;
    bool complete = false;
    std::vector<extension> next;
};

// should this be a class - so that you can do these things from within the data structure?
// this would probably allow for better data access...

struct range {
    int start = -1;
    int end = -1;
    int max_k = -1;
    int max_i = -1;
    int longest_forward_score = 0;
    int longest_reverse_score = 0;
    std::string longest_forward = "";
    std::string longest_reverse = "";
    std::vector<extension> forward;
    std::vector<extension> reverse;
};

bool compare_range(const range &a, const range &b)
{
    return a.max_k > b.max_k;
}


void print_extension( std::vector<extension> extension, int indent = 0 )
{
    if ( extension.size() == 0 )
    {
        printf("\n");
        return;
    }
    bool first = true;
    for ( auto e : extension )
    {
        if ( !first )
        {
            for (int i = 0; i < indent; i++)
                printf(" ");
        }
        printf(" -> %c", e.base);
        if (e.complete)
            printf("*");
        print_extension(e.next, indent + 5);
    }
}

class KmerBoundary {
    
  public:
    KmerBoundary() {}
    ~KmerBoundary() {}

};   // class KmerBoundary

class KmerBoundarySimple : KmerBoundary {
    JellyfishDatabase *db;
    enum sense_t { sense_none = 0, sense_AB, sense_BA };
    std::vector<int> kmer_content;
    std::vector<int> base_support;
    std::vector<range> start_positions;
    
  public:
    KmerBoundarySimple(JellyfishDatabase *a) { db = a; }
    ~KmerBoundarySimple() {}
    std::string describe() const {
        std::stringstream s;
        s << "A:" << db->describe();
        return s.str();
    }
    jmers::FosmidEndFragment detect_boundary(Seq& seq, bool error_correct = true, bool print_kmer_content = false)
    {
        
        /* Check kmer content of current sequence */
        kmer_content.clear();
        for (uint64_t i = 0; i < seq.sequence.length()-db->kmer; i++)
        {
            kmer_content.push_back( db->query( seq.sequence.substr(i, db->kmer) ) );
        }
        base_support.clear();
        base_support.resize( seq.sequence.length(), 0 );
        int o = 0;
        for ( auto k : kmer_content)
        {
            for (int i = 0; i < db->kmer; i++)
            {
                base_support[o+i] += k;
            }
            o++;
        }
        if ( print_kmer_content )
        {
            std::cout << ">" << seq.name << std::endl;
            for ( auto c : kmer_content )
                std::cout << c << " ";
            std::cout << std::endl;
        }
        
        if (error_correct)
        {
            this->error_correct(seq);
            kmer_content.clear();
            
            /* Update the kmer content */
            for (uint64_t i = 0; i < seq.sequence.length()-db->kmer; i++)
            {
                kmer_content.push_back(  db->query( seq.sequence.substr(i, db->kmer) ) );
            }
            // /* Update the base support */
            base_support.clear();
            base_support.resize( seq.sequence.length(), 0 );
            o = 0;
            for ( auto k : kmer_content)
            {
                for (int i = 0; i < db->kmer; i++)
                {
                    base_support[o+i] += k;
                }
                o++;
            }
        }

        if ( print_kmer_content )
        {
            for ( auto c : kmer_content )
                std::cout << c << " ";
            std::cout << std::endl;
            // for ( auto b : base_support )
            //     std::cout << b << " ";
            // std::cout << std::endl;
        }
        
        /* create fosmid end data structure */
        jmers::FosmidEndFragment fosmid_end(seq);
        //fosmid_end.infer_fragment_structure(kmer_content, base_support);
        return fosmid_end;
    }
private:
    
    int find_start_positions(int number_of_positions = -1)
    {
        /* This algorithm is basically:
         * 
         * 1) find all "islands" of non-zero kmer content stretches
         *
         * 2) find the "peak" (highest kmer content value) of each island
         *
         * 3) return the max(#peaks, number_of_positions) highest peaks
         */
        
        // 1) Find all non-zero kmer content stretches
        std::vector <range> non_zero_stretches;
        range r;
        int i = 0;
        for ( auto k : kmer_content )
        {
            if (k)
            {
                if (r.start < 0)
                {
                    r.start = i;
                }
            }
            else
                if (r.start >= 0)
                {
                    r.end = i-1;
                    non_zero_stretches.push_back(r);
                    r.start = -1;
                    r.end = -1;
                }
            i++;
        }
        if (r.start > 0 && r.end <= 0)
        {
            r.end = i-1;
            non_zero_stretches.push_back(r);
        }
        
        // 2) find the highest kmer content value of each range
        
        for (uint64_t c = 0; c < non_zero_stretches.size(); c++ )
        {
            range r = non_zero_stretches[c];
            for (int i = r.start; i <= r.end; i++)
                if ( kmer_content[i] > r.max_k )
                {
                    r.max_k = kmer_content[i];
                    r.max_i = i;
                }
            non_zero_stretches[c].max_k = r.max_k;
            non_zero_stretches[c].max_i = r.max_i;
        }
        
        // return the max(#peaks, number_of_positions) highest peaks
        std::sort(non_zero_stretches.begin(), non_zero_stretches.end(), compare_range);
        int p = 0;
        for ( auto &r : non_zero_stretches)
        {
            start_positions.push_back(r);
            if (number_of_positions > 0 && ++p >= number_of_positions)
                break;
        }
        return 0;
    }
    
    int extend_ranges( Seq& seq )
    {
        for ( auto &r : this->start_positions )
        {
            this->extend( r, r.forward, seq, 1 );
            this->extend( r, r.reverse, seq, -1 );
        }
        return 0;
    }
    
    int extend( range r, std::vector<extension> &extensions, Seq& seq, int step = 1, std::string path = "" )
    {
        // first, check if we're running out nof sequence
        if ( (step > 0 && path.length() + r.end >= seq.sequence.length()) || (step <= 0 && (long)r.start-(long)path.length() <= 0) )
            return 0;
        
        // second, check if we have any extensions
        for ( auto base : {'A','C','T','G'} )
        {
            // make some convenience variables
            int start = step > 0 ? r.end+1+path.length() : r.start-1-path.length();
            std::string kmer_base = seq.sequence.substr(start, db->kmer-1);
            // then check the current extensions kmer content
            int base_kmer_content = db->query( step > 0 ? kmer_base + base : base + kmer_base );
            
            // if we have no extension, just continue
            if ( base_kmer_content <= 0 )
                continue;
            
            // first - check if we can extend back to the original sequence, if we manage to do this, then we consider the 
            // position in the path "complete" and continue
            if ( base == seq.sequence[r.end+path.length()] )
            {
                extension e;
                e.base = base;
                e.kmer_support = base_kmer_content;
                e.complete = true;
                extensions.push_back(e);
                return 0;
            }
            else // ... otherwise - continue the path somewhere else
            {
                extension e;
                e.base = base;
                e.kmer_support = base_kmer_content;
                extensions.push_back(e);
            }
        }
        
        // if the path continues, but isn't 'complete', extend all exensions
        for ( auto &e : extensions)
        {
            this->extend( r, e.next, seq, step, path + e.base);
        }
        
        
        return 0;
    }
    
    int apply_extensions(Seq &seq)
    {
        for ( auto &r : this->start_positions )
        {
            // Find best forward extension
            std::string max_path = "";
            int max_support = 0;
            this->find_max_extension(r.forward, &max_path, &max_support);
            for (uint64_t i = 0; i < max_path.length(); i++)
                seq.sequence[r.end + i] = max_path[i];
            
            max_path = "";
            max_support = 0;
            // Find best reverse extension
            this->find_max_extension(r.reverse, &max_path, &max_support);
            for (uint64_t i = 0; i < max_path.length(); i++)
                seq.sequence[r.start - i] = max_path[i];
            
        }
        fflush(stdout);
        return 0;
    }
    
    void find_max_extension( std::vector<extension> extensions, std::string *max_path, int *max_support, std::string path = "", int support = 0, int level = 0)
    {
        if ( extensions.size() == 0 )
            return;
        
        for ( auto e : extensions )
        {
            std::string new_path = path.c_str();
            new_path += e.base;
            if ( new_path.length() > max_path->length() || (new_path.length() == max_path->length() && (support+e.kmer_support) > *max_support) )
            {
                *max_path = new_path.c_str();
                *max_support = support+e.kmer_support;
            }
            find_max_extension(e.next, max_path, max_support, new_path, support+e.kmer_support, level +1);
        }
    }
    
    int error_correct(Seq& seq)
    {
        // 1) - find some good 'reliable' positions to start error correcting
        this->find_start_positions();
        
        // 2)  - Attempt to fill in the regions without kmer support
        // 2.1 - first find alternative paths
        this->extend_ranges( seq );
        
        // 2.2 - then choose the 'best' path and apply it to the sequence
        this->apply_extensions( seq );
        
        return 0;
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
        s << "A:" << db->describe() << ", B:" << B.describe();
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


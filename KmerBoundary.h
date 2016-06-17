// GPLv3 header


#include <string>
#include <cctype>
#include <iostream>
#include <sstream>

namespace jmers {

class KmerDB {
    DB;
    bool        is_empty;
    std::string description;
    KmerDB(...)  // initialise with jellyfish db?
    { }
    KmerDB()
        : is_empty(true)
    { }
    std::string describe() const {
        if (is_empty)
            return "empty";
        return description;
    }
    bool is_present(std::string& k, int tol1, double tol2) const {
        // query DB for presence of k with tolerances tol1 and tol2
        return false;
    }
};  // class KmerDB

KmerDB empty_db = KmerDB();

// perhaps a better design would be KmerBoundaryChimeric/KmerBoundarySimilar,
// where either side of the boundary is from the same kmer database, and
// KmerBoundaryDissimilar/KmerBoundaryDistinct, where the two sides are from
// different databases, one of them possibly being no database

// KmerBoundarySimilar would not have a sense
// KmerboundaryDissimilar would

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
        stringstream s;
        s << "A:" << A.describe() << ", B:" << B.describe();
        return s.str();
    }
    int64_t detect_boundary(Seq& seq, sense_t sns, int64_t pos = 0,
            int tol1 = 0, double tol2 = 0.0d0) {
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

}  // namespace jmers


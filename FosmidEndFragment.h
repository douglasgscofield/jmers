// TODO: GPLv3 copyright notice header

// TODO: infer_fragment_structure()
// TODO: split_fragment()
// TODO: perhaps adjust base qualities in vicinity of inferred pos-es? 
// TODO: split_fragment(): perhaps create sequential read name? 
// TODO: think through placement of end_pad and other namespace globals
// TODO: abstract write_pair_fast{a,q} so output format needn't be known here

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>
#include <iostream>
#include <sstream>

#include "Seq.h"

namespace jmers {

class FosmidEndFragment {
    Seq fragment;
    int64_t end1_pos;      // leftmost inferred start of genomic sequence in fragment
    // if inference of ligation was 'perfect', then ligation1_pos and ligation2_pos are
    // identical
    int64_t ligation1_pos; // leftmost site of ligation between fosmid ends in fragment
    int64_t ligation2_pos; // rightmost site of ligation between fosmid ends in fragment
    int64_t end2_pos;      // rightmost inferred end of genomic sequence in fragment
    Seq read1;  // filled by split_fragment
    Seq read2;  // filled by split_fragment

  public:

    FosmidEndFragment(const Seq& s)
        : fragment(s), end1_pos(-1), ligation1_pos(-1), ligation2_pos(-1), end2_pos(-1)
    {
        // is default copy constructor OK?
    }
    ~FosmidEndFragment()  // default destructor OK
    { }

    void infer_fragment_structure() {
        // fill end1_pos
        // fill ligation1_pos
        // fill ligation2_pos
        // fill end2_pos
    }

    void split_fragment() {
        // use inferred *_pos to fill read1 and read2 from fragment
        if (end1_pos < 0 || ligation1_pos < 0 || ligation2_pos < 0 || end2_pos < 0
                || end1_pos >= ligation1_pos || ligation2_pos >= end2_pos
                || ligation1_pos > ligation2_pos) {
            std::cerr << "jmers::FosmidEndFragment::split_fragment: pos values incorrect" << std::endl;
            dump(); exit(1);
        }
        // read name
        std::stringstream nm_ss;
        nm_ss << fragment.name << ":FFE";
        // read boundaries
        int64_t read1_beg = end1_pos      + end_pad;
        int64_t read1_len = ligation1_pos - end_pad - read1_beg;
        int64_t read2_beg = ligation2_pos + end_pad;
        int64_t read2_len = end2_pos      - end_pad - read2_beg;
        if (fragment.has_quality) {
            read1.fill(nm_ss.str(), "1", fragment.sequence.substr(read1_beg, read1_len),
                       fragment.quality.substr(read1_beg, read1_len));
            read2.fill(nm_ss.str(), "2", fragment.sequence.substr(read2_beg, read2_len),
                       fragment.quality.substr(read2_beg, read2_len));
        } else {
            read1.fill(nm_ss.str(), "1", fragment.sequence.substr(read1_beg, read1_len));
            read2.fill(nm_ss.str(), "2", fragment.sequence.substr(read2_beg, read2_len));
        }
    }
    void write_pair_fastq(ostream& os = std::cout) const {
        // where is output going...
        read1.write_fastq(read1_ostream);
        read2.write_fastq(read2_ostream);
    }
    void write_pair_fasta(ostream& os = std::cout) const {
        // where is output going...
        read1.write_fasta(read1_ostream);
        read2.write_fasta(read2_ostream);
    void dump(std::ostream& os = std::cerr) const {
        os << "jmers::FosmidEndFragment::dump:" << std::endl;
        os << "    fragment="; fragment.dump(os);
        os << "    end1_pos=" << end1_pos << " ligation1_pos=" << ligation1_pos
            << " ligation2_pos=" << ligation2_pos << " end2_pos=" << end2_pos << std::endl;
        os << "    read1="; read1.dump(os);
        os << "    read2="; read2.dump(os);
    }
}; // class FosmidEndFragment

} // namespace jmers


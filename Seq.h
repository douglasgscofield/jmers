// TODO: GPLv3 header
// TODO: add pairing?
// TODO: linewidth for write_fasta()

#include <string>
#include <cctype>
#include <iostream>
extern "C" {
// Heng Li's kseq.h for reading Fasta and FastQ
// http://lh3lh3.users.sourceforge.net/kseq.shtml
// carries MIT License
// here we just need kseq_t
#include "kseq/kseq.h"
}
namespace jmers {

struct Seq {
    std::string name, comment, sequence, quality;
    size_t l;
    bool has_quality;
    Seq()  // default constructor
        : l(0), has_quality(false)
    { }
    Seq(kseq_t *kseq_seq)  // constructor from kseq.h seq
    {
        fill(kseq_seq);
    }
    ~Seq() { }
    void fill(kseq_t *kseq_seq) {
        if (kseq_seq->name.l) name.assign(kseq_seq->name.s);
        if (kseq_seq->comment.l) comment.assign(kseq_seq->comment.s);
        if (kseq_seq->seq.l) {
            sequence.assign(kseq_seq->seq.s);
            l = sequence.length();
        } else l = 0;
        if (kseq_seq->qual.l) {
            quality.assign(kseq_seq->qual.s);
            has_quality = true;
        } else has_quality = (l == 0);
    }
    void dump(std::ostream& os) const {
        os << "jmers::Seq :";
        os << " l=" << l << " has_quality=" << has_quality;
        os << " name=" << name << " comment=" << comment;
        os << " sequence=" << sequence << " quality=" << quality;
        os << std::endl;
    }
    void write_fastq(std::ostream& os) const {
        if (! has_quality) { std::cerr << "jmers::seq.write_fastq: no quality string" << std::endl; exit(1); }
        os << "@" << name;
        if (comment.length()) os << " " << comment;
        os << std::endl;
        os << sequence << std::endl;
        os << "+" << std::endl;
        os << quality << std::endl;
    }
    void write_fasta(std::ostream& os, int linewidth = 0) const {
        os << ">" << name;
        if (comment.length()) os << " " << comment;
        os << std::endl;
        os << sequence << std::endl;
    }
}; // struct seq

} // namespace jmers


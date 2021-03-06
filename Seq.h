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

// TODO: add pairing?
// TODO: linewidth for write_fasta()

#ifndef __JMERS_SEQ_H__
#define __JMERS_SEQ_H__

#include "jmers.h"
#include "zlib.h"
#include <string>
#include <cctype>
#include <iostream>
extern "C" {
// Heng Li's kseq.h for reading Fasta and FastQ
// http://lh3lh3.users.sourceforge.net/kseq.shtml
// carries MIT License
// here we just need kseq_t
#include "kseq/kseq.h"
KSEQ_INIT(gzFile, gzread);
}
namespace jmers {


class Seq {
  public:
    std::string name, comment, sequence, quality;
    size_t l;          // sequence length shortcut, so l == 0 if no sequence
    bool has_quality;  // false when sequence with no quality, so true if no sequence
    Seq()  // default constructor
        : l(0), has_quality(true)
    { }
    Seq(kseq_t *kseq_seq)  // constructor using kseq.h kseq_t
    {
        fill(kseq_seq);
    }
    ~Seq() { }
    void extend(char base, int quality = 20, int quality_offset = 33)
    {
        sequence += base;
        quality += (char)(quality + quality_offset);
    }
    void extend(std::string bases, int quality = 20, int quality_offset = 33)
    {
        sequence += bases;
        for (auto c : bases)  // this assumes each base has default quality
            quality += (char)(quality + quality_offset);
        l += bases.length();
    }
    
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
    void fill(const std::string& n = empty, const std::string& c = empty,
              const std::string& s = empty, const std::string& q = empty) {
        name = n;
        comment = c;
        sequence = s;
        quality = s;
        l = sequence.length();
        if (quality.length())
            has_quality = true;
        else has_quality = (l == 0);
    }
    void dump(std::ostream& os = std::cerr) const {
        os << "jmers::Seq::dump:";
        os << " l=" << l << " has_quality=" << has_quality;
        os << " name=" << name << " comment=" << comment;
        os << " sequence=" << sequence << " quality=" << quality;
        os << std::endl;
    }
    void write_fastq(std::ostream& os = std::cout) const {
        if (! has_quality) {
            std::cerr << "jmers::Seq::write_fastq: no quality string" << std::endl;
            exit(1);
        }
        os << "@" << name;
        if (comment.length()) os << " " << comment;
        os << std::endl << sequence << std::endl << "+" << std::endl << quality << std::endl;
    }
    void write_fasta(std::ostream& os = std::cout, int linewidth = 0) const {
        os << ">" << name;
        if (comment.length()) os << " " << comment;
        os << std::endl << sequence << std::endl;
    }
}; // class Seq

} // namespace jmers

#endif

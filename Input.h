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

// TODO: read bzip2 
// TODO: read xz
// TODO: read fqzcomp ?
// TODO: check read from stdin
// TODO: check read from fifo
// TODO: learn exceptions :-)
// TODO: lift kseq.h from subdirectory

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>
#include <iostream>
#include <fstream>

#include "Seq.h"

namespace jmers {

extern "C" {
#include <zlib.h>
// Heng Li's kseq.h for reading Fasta and FastQ
// http://lh3lh3.users.sourceforge.net/kseq.shtml
// carries MIT License
#include "kseq/kseq.h"
}

class Input {
    std::string filename;  // filename
    gzFile fp;             // filepointer
    kseq_t *kseq_seq;      // handle for kseq.h
  public:
    Input(std::string f)  { open(f.c_str()); }
    Input(const char* fn) { open(fn); }
    ~Input()              { close(); }
    void open(const char* fn) {
        if (! fn) { std::cerr << "no filename" << std::endl; exit(1); }
        filename.assign(fn);
        fp = gzopen(fn, "r");
        if (! fp) { std::cerr << "could not open file " << filename << std::endl; exit(1); }
        kseq_seq = kseq_init(fp);
    }
    bool read(Seq& s) {  // jmers::Seq from Seq.h
        int64_t l = kseq_read(kseq_seq);
        if (l < 0) return(false);
        s.fill(kseq_seq);
        return(true);
    }
    void close() {
        kseq_destroy(kseq_seq);
        gzclose(fp);
    }
}; // class Input

} // namespace jmers


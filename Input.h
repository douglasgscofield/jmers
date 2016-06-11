// TODO: GPLv3 copyright notice header

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


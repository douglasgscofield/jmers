// GPLv3 header

// global config info within namespace jmers
// should this also include the other jmers includes to be the primary include file?

#include <string>
#include <cctype>
#include <iostream>

namespace jmers {

int64_t end_pad = 3; // when splitting fragment, trim this many bp in from each end pos

const std::string empty = std::string();  // empty string utility constant

std::ostream read1_ostream = std::cout;  // for read1 output
std::ostream read2_ostream = std::cout;  // for read2 output

}  // namespace jmers

/*!
 * \file   BubbleDescription.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include "MGIS/Raise.hxx"
#include "OperaHPC/BubbleDescription.hxx"

namespace opera_hpc {

  std::vector<BubbleDescription> BubbleDescription::read(const std::string& f) {
    auto tokenize = [](const std::string& line) {
      std::istringstream tokenizer(line);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(tokenizer),
                std::istream_iterator<std::string>(),
                std::back_inserter(tokens));
      return tokens;
    };
    auto bubbles = std::vector<BubbleDescription>{};
    std::ifstream in(f);
    if (!in) {
      mgis::raise("can't read file '" + std::string{f} + "'");
    }
    auto ln = mfem_mgis::size_type{};
    std::string line;
    while (getline(in, line)) {
      ++ln;
      auto tokens = tokenize(line);
      auto throw_if = [&f, &line, ln](const bool b,
                                      const std::string_view msg) {
        if (b) {
          mgis::raise("error at line '" + std::to_string(ln) +
                      "' while reading file '" + f + "' (" + std::string{msg} +
                      ")");
        }
      };
      if (tokens.empty()) {
        continue;
      }
      if (tokens[0][0] == '#') {
        continue;
      }
      // id center_x center_y center_z r
      throw_if(tokens.size() != 5, "ill-formed line '" + line + "'");
      //
      auto convert = [throw_if](auto& v, const std::string& w) {
        std::istringstream converter(w);
        converter >> v;
        throw_if(!converter || (!converter.eof()), "conversion failed");
      };
      auto bubble = BubbleDescription{};
      convert(bubble.boundary_identifier, tokens[0]);
      convert(bubble.center[0], tokens[1]);
      convert(bubble.center[1], tokens[2]);
      convert(bubble.center[2], tokens[3]);
      convert(bubble.radius, tokens[4]);
      bubbles.push_back(bubble);
    }
    return bubbles;
  }  // end of read

  mfem_mgis::real distance(const BubbleDescription& b,
                           const std::array<mfem_mgis::real, 3u>& p) noexcept {
    auto square = [](const mfem_mgis::real x) { return x * x; };
    return std::sqrt(square(b.center[0] - p[0]) +  //
                     square(b.center[1] - p[1]) +  //
                     square(b.center[2] - p[2]));
  }

}  // end of namespace opera_hpc

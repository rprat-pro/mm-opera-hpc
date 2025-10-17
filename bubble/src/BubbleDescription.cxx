/*!
 * \file   BubbleDescription.cxx
 * \brief  Implementation of bubble description utilities for reading and manipulating bubble data
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

  /**
   * \brief Read bubble descriptions from a text file
   * 
   * File format: Each line contains 5 space-separated values:
   * boundary_id center_x center_y center_z radius
   * 
   * Lines starting with '#' are treated as comments and ignored.
   * Empty lines are also ignored.
   * 
   * \param f Path to the file containing bubble descriptions
   * \return Vector of BubbleDescription objects parsed from the file
   * \throw std::runtime_error if file cannot be read or line format is invalid
   * 
   * Example file content:
   * # bubble_id x y z radius
   * 1 0.0 0.0 0.0 0.5
   * 2 1.0 1.0 1.0 0.3
   */
  std::vector<BubbleDescription> BubbleDescription::read(const std::string& f) {
    /**
     * \brief Function to tokenize a string into words
     * 
     * Splits a line by whitespace into individual tokens
     * 
     * \param line Input string to tokenize
     * \return Vector of string tokens
     */
    auto tokenize = [](const std::string& line) {
      std::istringstream tokenizer(line);
      std::vector<std::string> tokens;
      // Copy all whitespace-separated tokens into the vector
      std::copy(std::istream_iterator<std::string>(tokenizer),
                std::istream_iterator<std::string>(),
                std::back_inserter(tokens));
      return tokens;
    };
    
    // Container for all bubble descriptions
    auto bubbles = std::vector<BubbleDescription>{};
    
    // Open the input file
    std::ifstream in(f);
    if (!in) {
      mgis::raise("can't read file '" + std::string{f} + "'");
    }
    
    // Line counter for error reporting
    auto ln = mfem_mgis::size_type{};
    std::string line;
    
    // Process file line by line
    while (getline(in, line)) {
      ++ln;
      auto tokens = tokenize(line);
      
      /**
       * \brief Lambda to throw error with context information
       * 
       * \param b Boolean condition - throws if true
       * \param msg Error message to display
       */
      auto throw_if = [&f, &line, ln](const bool b,
                                      const std::string_view msg) {
        if (b) {
          mgis::raise("error at line '" + std::to_string(ln) +
                      "' while reading file '" + f + "' (" + std::string{msg} +
                      ")");
        }
      };
      
      // Skip empty lines
      if (tokens.empty()) {
        continue;
      }
      
      // Skip comment lines (starting with '#')
      if (tokens[0][0] == '#') {
        continue;
      }
      
      // Validate line format: must have exactly 5 tokens
      // Format: id center_x center_y center_z radius
      throw_if(tokens.size() != 5, "ill-formed line '" + line + "'");
      
      /**
       * \brief Lambda to convert string token to numeric value
       * 
       * Uses string stream for safe conversion with error checking
       * 
       * \param v Reference to variable to store converted value
       * \param w String token to convert
       */
      auto convert = [throw_if](auto& v, const std::string& w) {
        std::istringstream converter(w);
        converter >> v;
        // Ensure conversion succeeded and entire token was consumed
        throw_if(!converter || (!converter.eof()), "conversion failed");
      };
      
      // Parse bubble data from tokens
      auto bubble = BubbleDescription{};
      convert(bubble.boundary_identifier, tokens[0]);  // Boundary ID (integer)
      convert(bubble.center[0], tokens[1]);            // X coordinate of center
      convert(bubble.center[1], tokens[2]);            // Y coordinate of center
      convert(bubble.center[2], tokens[3]);            // Z coordinate of center
      convert(bubble.radius, tokens[4]);               // Bubble radius
      
      bubbles.push_back(bubble);
    }
    
    return bubbles;
  }  // end of read

  /**
   * \brief Calculate Euclidean distance from bubble center to a point
   * 
   * Computes the 3D Euclidean distance between the center of a bubble
   * and an arbitrary point in space. This is used to determine if a
   * point is associated with a particular bubble.
   * 
   * \param b Bubble description containing center coordinates
   * \param p 3D coordinates of the point [x, y, z]
   * \return Euclidean distance from bubble center to point
   * 
   * Formula: sqrt((x_b - x_p)^2 + (y_b - y_p)^2 + (z_b - z_p)^2)
   */
  mfem_mgis::real distance(const BubbleDescription& b,
                           const std::array<mfem_mgis::real, 3u>& p) noexcept {
    
    auto square = [](const mfem_mgis::real x) { return x * x; };
    
    // Compute 3D Euclidean distance
    return std::sqrt(square(b.center[0] - p[0]) +  // (x_b - x_p)^2
                     square(b.center[1] - p[1]) +  // (y_b - y_p)^2
                     square(b.center[2] - p[2]));  // (z_b - z_p)^2
  }

}  // end of namespace opera_hpc
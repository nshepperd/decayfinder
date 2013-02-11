#include <framework/logging/Logger.h>
#include <generators/dataobjects/MCParticle.h>
#include <framework/datastore/StoreArray.h>
#include <TDatabasePDG.h>
#include <analysis/modules/DecayFinder/Descriptor.h>
using namespace Belle2;

#include <utility>
#include <string>
#include <list>
#include <set>

bool matching(std::set<int> left, std::set<int> right, std::list< std::pair<int, int> > edges);

DecayDescriptor::DecayDescriptor(Descriptor* origin, std::vector<Descriptor*> decays, std::string arrow, bool inclusive) {
  this->origin = origin;
  this->decays = decays;
  this->arrow = arrow;
  this->inclusive = inclusive;
}

std::string DecayDescriptor::repr() const {
  std::string result;
  result += origin->repr() + " " + arrow;
  for(unsigned int i = 0; i < decays.size(); i++) {
    result += " " + decays[i]->repr();
  }
  return "(" + result + ")";
}


bool DecayDescriptor::match(int index) const {
  StoreArray<MCParticle> mcparticles;
  MCParticle* p = mcparticles[index];

  if(!this->origin->match(index)) {
    return false;
  }

  int desc_num = this->decays.size();

  // getFirstDaughter/getLastDaughter return the particle index + 1
  // getFirstDaughter() is 0 if there are no children
  int mcp_first = p->getFirstDaughter() - 1;
  int mcp_num = 0;
  if(mcp_first >= 0) {
    mcp_num = p->getLastDaughter() - mcp_first;
  }

  // If we have different numbers of particles, it clearly can't match
  // (until we do inclusive matchings)
  if(mcp_num < desc_num) {
    return false;
  }

  // Edges connects mcp children to matching descriptor children.
  std::list< std::pair<int, int> > edges;
  for(int i = 0; i < mcp_num; i++) {
    for(int j = 0; j < desc_num; j++) {
      if(this->decays[j]->match(mcp_first + i)) {
	edges.push_front(std::pair<int, int>(i, j));
      }
    }
  }

  // There might be a more efficient way to represent left and right than with set objects...
  std::set<int> left, right;
  if(this->inclusive) {
    // For an inclusive decay, leave 'left' empty. None of the monte carlo
    // children are 'required' to be matched.
  } else if(this->arrow == "=>") {
    // Gamma-inclusive. Only non-gammas must be matched.
    for(int i = 0; i < mcp_num; i++) {
      if(mcparticles[mcp_first + i]->getPDG() != 22) {
	left.insert(i);
      }
    }
  } else {
    // Exclusive decay. Everything has to be matched.
    for(int i = 0; i < mcp_num; i++) {
      left.insert(i);
    }
  }

  // So far, no descriptor children are optional.
  for(int j = 0; j < desc_num; j++) {
    right.insert(j);
  }

  return matching(left, right, edges);
}

// Return whether the edge `value` is compatible to `base` (doesn't share any edges).
// Constructor takes `base` as a parameter.
class CompatibleFilter {
private:
  const std::pair<int, int> base;
public:
  CompatibleFilter(const std::pair<int, int> base) : base(base) {}
  bool operator()(const std::pair<int, int> value) const {
    return (value.first != base.first) && (value.second != base.second);
  }
};

bool matching(std::set<int> left, std::set<int> right, std::list< std::pair<int, int> > edges) {
  // The problem of matching up mcparticle children to the decay descriptor tree
  // is actually the 'Bipartite Graph Matching' problem. 'left' and 'right'
  // are a partitioning of a graph into two disjoint subsets such that all edges
  // go between an element of 'left' and an element of 'right'.

  // In this case, 'left' is the monte carlo children, 'right' is the descriptor
  // children, and 'edges' connects those that match.

  // I could use a real algorithm with decent asymptotic properties, but since my
  // input is usually of size < 5, something like the Hopcroft-Karp algorithm would
  // probably be slower than naive brute force, due to complexity.

  // So essentially I do a depth first search, for each edge in some set of mutually
  // incompatible edges, selecting that edge, and recursing on the same problem with
  // that edge, its endpoints, and any other edges connected to those endpoints removed.

  // The procedure yields each subset of 'edges' that forms a valid bijection between
  // 'left' and 'right'. That is, all valid ways of matching the monte carlo children
  // to the descriptor. In this case, I only care about whether such a subset exists.

  if(left.empty() && right.empty()) {
    // All vertices paired up.
    return true;
  } else if(edges.empty()) {
    // Out of compatible edges
    return false;
  }

  std::list< std::pair<int, int> > incompatible(edges);
  while(!incompatible.empty()) {
    // pop an edge at random
    std::pair<int, int> e = incompatible.front();
    incompatible.pop_front();

    // remove incompatible edges and vertices
    std::set<int> new_left(left);
    new_left.erase(e.first);
    std::set<int> new_right(right);
    new_right.erase(e.second);

    std::list< std::pair<int, int> > new_edges;
    std::list< std::pair<int, int> >::iterator it;
    for(it = edges.begin(); it != edges.end(); it++) {
      if((it->first != e.first) && (it->second != e.second)) {
	new_edges.push_front(*it);
      }
    }

    // recurse on the reduced problem
    if(matching(new_left, new_right, new_edges)) {
      return true;
    }

    // then we try the next edge that is incompatible to the ones we tried
    incompatible.remove_if(CompatibleFilter(e));
  }

  return false;
}

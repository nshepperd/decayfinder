#include <framework/logging/Logger.h>
#include <generators/dataobjects/MCParticle.h>
#include <framework/datastore/DataStore.h>
#include <framework/datastore/StoreArray.h>
#include <TDatabasePDG.h>
#include <analysis/modules/DecayFinder/Descriptor.h>
using namespace Belle2;

#include <string>

AtomicDescriptor::AtomicDescriptor(std::string name) {
  this->name = name;

  // special names
  if(name == "X" || name == "X+" || name == "X-") {
    this->pdgcode = 0;
    return;
  }

  const char* str = name.c_str();

  TDatabasePDG* pdb = TDatabasePDG::Instance();
  TParticlePDG* entry = pdb->GetParticle(str);
  if(entry) {
    this->pdgcode = entry->PdgCode();
  } else {
    B2FATAL("Unrecognised particle name: " << name);
  }
}

bool AtomicDescriptor::match(int index) const {
  StoreArray<MCParticle> mcparticles;
  if(this->pdgcode == 0) {
    // Special descriptors.
    if(this->name == "X") {
      return true;
    } else if(this->name == "X+") {
      return mcparticles[index]->getCharge() > 0;
    } else if(this->name == "X-") {
      return mcparticles[index]->getCharge() < 0;
    }
  }

  // Just match particle id.
  return (this->pdgcode == mcparticles[index]->getPDG());
}

std::string AtomicDescriptor::repr() const {
  return this->name;
}

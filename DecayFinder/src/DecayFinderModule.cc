#include <framework/logging/Logger.h>
#include <generators/dataobjects/MCParticle.h>
#include <framework/datastore/DataStore.h>
#include <framework/datastore/StoreArray.h>
#include <TDatabasePDG.h>

#include <jansson.h>
#include <set>
#include <list>
#include <utility>
#include <vector>

#include <analysis/modules/DecayFinder/DecayFinderModule.h>


using namespace Belle2;

static void print_tree(const int i, const int level);
static Descriptor* convert_json(json_t* desc);

// Register with framework
REG_MODULE(DecayFinder)

// DecayFinderModule
DecayFinderModule::DecayFinderModule() : Module() {
  this->setDescription("New decay finder.");
  this->addParam("pattern", this->m_pattern, "Decay tree pattern.", std::string(""));
}

DecayFinderModule::~DecayFinderModule() { }

void DecayFinderModule::initialize() {
  json_error_t error;
  json_t* json = json_loads(this->m_pattern.c_str(), 0, &error);
  if(json == NULL) {
    B2FATAL("JSON decoding error.");
  }

  this->m_descriptor = convert_json(json);
  printf("descriptor: %s\n", this->m_descriptor->repr().c_str());
}

void DecayFinderModule::terminate() { }

void DecayFinderModule::event() {
  StoreArray<MCParticle> mcparticles;
  for(int i = 0; i < mcparticles.getEntries(); i++) {
    if(this->m_descriptor->match(i)) {
      print_tree(i, 0);
    }
  }
}

// print_tree:
// Recursively print out in pretty colours the decay tree starting at particle i.
// 'level' gives a level of indentation for the output, only really used in internal recursion 
static void print_tree(const int i, const int level) {
  StoreArray<MCParticle> mcparticles;
  TDatabasePDG* pdb = TDatabasePDG::Instance();

  MCParticle* p = mcparticles[i];
  int j;

  int pid = p->getPDG();
  double mass = p->getMass();
  double charge = p->getCharge();
  double energy = p->getEnergy();

  TParticlePDG* entry = pdb->GetParticle(pid);
  std::string name = "unknown";
  if(entry) {
    name = std::string(entry->GetName());
  }

  std::string prefix, suffix;
  if(0 < level && level < 7) {
    char buffer[6];
    snprintf(buffer, 6, "\033[3%im", level);
    prefix = std::string(buffer);
    suffix = "\033[m";
  }

  for(j = 0; j < level; j++) {
    printf("    ");
  }

  printf("%s[%i] %s mass=%f energy=%f charge=%f%s\n", prefix.c_str(), pid, name.c_str(), mass, energy, charge, suffix.c_str());

  if(p->getFirstDaughter() > 0) {
    int m = p->getFirstDaughter();
    int n = p->getLastDaughter();
    for(j = m-1; j <= n-1; j++) {
      print_tree(j, level + 1);
    }
  }
}

static Descriptor* convert_json(json_t* desc) {
  // Turn a json object into a C++ Descriptor object

  if(!json_is_object(desc)) {
    B2FATAL("Descriptor is not object.");
  }

  std::string type = json_string_value(json_object_get(desc, "type"));
  if(type == "Atomic") {
    // Atomic descriptor: particle name / property
    return new AtomicDescriptor(json_string_value(json_object_get(desc, "name")));
  } else if(type == "Decay") {
    // Decay descriptor: (A -> B C D [...])
    Descriptor* origin = convert_json(json_object_get(desc, "origin"));

    // convert all the decay children
    json_t* jsondecay = json_object_get(desc, "decays");
    int num_decays = json_array_size(jsondecay);
    std::vector<Descriptor*> decays(num_decays);
    for(int i = 0; i < num_decays; i++) {
      decays[i] = convert_json(json_array_get(jsondecay, i));
    }

    return new DecayDescriptor(origin, decays,
			       json_string_value(json_object_get(desc, "arrow")),
			       json_is_true(json_object_get(desc, "inclusive")));
  } else if(type == "Logical") {
    std::string op = json_string_value(json_object_get(desc, "op"));
    Descriptor* left = convert_json(json_object_get(desc, "left"));
    Descriptor* right = convert_json(json_object_get(desc, "right"));
    if(op == "||") {
      return new OrDescriptor(left, right);
    } else if(op == "&&") {
      return new AndDescriptor(left, right);
    } else {
      B2FATAL("Unknown logical operator: " << op);
    }
  }
  // else if(type == "List") {
  //   return match_list(index, descriptor);
  // }

  return NULL;
}


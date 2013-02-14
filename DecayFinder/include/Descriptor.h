/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2010 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Neil Shepperd <nshepperd@gmail.com>                      *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#ifndef DECAYFINDER_DESCRIPTOR_H
#define DECAYFINDER_DESCRIPTOR_H

#include <vector>
#include <string>

namespace Belle2 {

  // Matches a certain kind of particle/decay.
  class Descriptor {
  public:
    virtual bool match(int index) const = 0;
    virtual std::string repr() const = 0;
    virtual ~Descriptor() {}
  };

  // Matches a named particle, or class of particle (eg. B0, X+, baryon, spin=1).
  class AtomicDescriptor : public Descriptor {
  public:
    AtomicDescriptor(std::string name);
    virtual std::string repr() const;
    virtual bool match(int index) const;
    virtual ~AtomicDescriptor() {}
  private:
    std::string name;
    int pdgcode;
  };

  // Matches a particular decay, A -> B C ..., where A, B, C, ... are Descriptors
  class DecayDescriptor : public Descriptor {
  public:
    DecayDescriptor(Descriptor* origin, std::vector<Descriptor*> decays, std::string arrow, bool inclusive);
    virtual std::string repr() const;
    virtual bool match(int index) const;
    virtual ~DecayDescriptor() {}
  private:
    Descriptor* origin;
    std::vector<Descriptor*> decays;
    std::string arrow;
    bool inclusive;
  };

  // A || B, A && B
  class LogicalDescriptor : public Descriptor {
  public:
    LogicalDescriptor(Descriptor* left, Descriptor* right) : left(left), right(right) {}
    virtual ~LogicalDescriptor() {}
  protected:
    Descriptor* left;
    Descriptor* right;
  };

  class OrDescriptor : public LogicalDescriptor {
  public:
    OrDescriptor(Descriptor* left, Descriptor* right) : LogicalDescriptor(left, right) {}
    virtual bool match(int index) const;
    virtual std::string repr() const;
  };

  class AndDescriptor : public LogicalDescriptor {
  public:
    AndDescriptor(Descriptor* left, Descriptor* right) : LogicalDescriptor(left, right) {}
    virtual bool match(int index) const;
    virtual std::string repr() const;
  };

}

#endif // DECAYFINDER_DESCRIPTOR_H

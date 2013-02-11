/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2010 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Neil Shepperd <nshepperd@gmail.com>                      *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#ifndef DECAYFINDER_MODULE_H
#define DECAYFINDER_MODULE_H

#include <framework/core/Module.h>
#include <analysis/modules/DecayFinder/Descriptor.h>

#include <string>

namespace Belle2 {

  class DecayFinderModule : public Module {

  public:

    DecayFinderModule();

    /** Destructor. */
    virtual ~DecayFinderModule();

    virtual void initialize();
    virtual void terminate();

    /** Method is called for each event. */
    virtual void event();

  private:
    std::string m_pattern;

    Descriptor* m_descriptor;
  };

}

#endif // DECAYFINDER_MODULE_H

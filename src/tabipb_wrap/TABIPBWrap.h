#ifndef H_TABIPB_WRAP_H
#define H_TABIPB_WRAP_H

#include "TABIPBStruct.h"

#ifdef TABIPB_APBS
#include "generic/valist.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef TABIPB_APBS
TABIPBOutput runTABIPBWrapAPBS(TABIPBInput tabipbIn, Valist* APBSMolecule);
#endif

#ifdef __cplusplus
}
#endif

#endif

#ifndef H_TABIPB_WRAP_H
#define H_TABIPB_WRAP_H

#include "TABIPBStruct.h"
#include "generic/valist.h"

#ifdef __cplusplus
extern "C" {
#endif

TABIPBOutput runTABIPBWrapAPBS(TABIPBInput tabipbIn, Valist* APBSMolecule);

#ifdef __cplusplus
}
#endif

#endif

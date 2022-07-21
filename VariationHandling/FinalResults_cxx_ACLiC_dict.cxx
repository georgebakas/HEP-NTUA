// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIgeorgebakasdIDocumentsdIHEPmINTUA_uldIVariationHandlingdIFinalResults_cxx_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling/./FinalResults.cxx"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_FinalResults_cxx_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./FinalResults.cxx",
0
    };
    static const char* includePaths[] = {
"/Applications/root_v6.22.02/include",
"/Applications/root_v6.22.02/etc/",
"/Applications/root_v6.22.02/etc//cling",
"/Applications/root_v6.22.02/include/",
"/Applications/root_v6.22.02/include",
"/Applications/root_v6.22.02/include/",
"/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FinalResults_cxx_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FinalResults_cxx_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./FinalResults.cxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"AddSystematicToErrorBar", payloadCode, "@",
"DrawWithRatio", payloadCode, "@",
"FinalResults", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FinalResults_cxx_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FinalResults_cxx_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FinalResults_cxx_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FinalResults_cxx_ACLiC_dict() {
  TriggerDictionaryInitialization_FinalResults_cxx_ACLiC_dict_Impl();
}

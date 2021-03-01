#ifndef _BASE_H_
#define _BASE_H_
#define _UL_

#ifdef _UL_
#include "AnalysisConstants_UL.h"
namespace AnalysisConstants = AnalysisConstants_UL;
#else
#include "AnalysisConstants.h"
#endif

#include "AnalysisUtilities.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"

#endif
#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <vector>
#include <map>

#include <TString.h>

namespace Combination
{
  std::map<TString, std::vector<float>> legendCoordinates = {{"leadingJetPt", {0., 0., 0., 0.}}};
}

#endif
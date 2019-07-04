#include "BoostedAnalysisInitializer.h"
R__LOAD_LIBRARY(BoostedAnalysisInitializer_cxx)

void test(float deepAK8Cut = 0.6, float tTaggerCut = 0.2, float mvaCut = 0.8)
{
  BoostedAnalysisInitializer utils("2017");
  
  for(int i=0; i<utils.listOfFiles.size(); i++)
  {
    std::cout<<utils.listOfFiles[i]<<std::endl;
  }
}
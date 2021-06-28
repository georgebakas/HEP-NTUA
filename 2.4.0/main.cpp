#include "CombineMeasurements.cxx"

int main()
{
  try
  {
    std::cout << "Start!!" << std::endl;

    //TFile *outFile = new TFile("outFileParton.root", "RECREATE");
    TFile *outFile = new TFile("outFileFiducial.root", "RECREATE");
    CombineMeasurements(outFile);
    //CombineFiducialMeasurements(outFile, true);
    outFile->Close();
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << '\n';
  }
  return 0;
}
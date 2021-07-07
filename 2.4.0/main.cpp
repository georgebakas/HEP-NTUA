#include "CombineMeasurements.cxx"
#include "CombineUnfoldedMeasurements.cxx"
int main()
{
  try
  {
    std::cout << "Start!!" << std::endl;

    TFile *outFile = new TFile("outFileParton.root", "RECREATE");
    //TFile *outFile = new TFile("outFileFiducial.root", "RECREATE");
    //CombineFiducialMeasurements(outFile, true);
    //CombineFiducialMeasurements(outFile);
    CombineUnfoldedMeasurements(outFile, "Parton");

    outFile->Close();
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << '\n';
  }
  return 0;
}

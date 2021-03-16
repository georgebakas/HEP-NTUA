#include "CombineMeasurements.cxx"

int main()
{
  try
  {
    std::cout << "Start!!" << std::endl;

    TFile *outFile = new TFile("outFile.root", "RECREATE");
    CombineMeasurements(outFile);
    CombineMeasurements(outFile, true);
    outFile->Close();
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << '\n';
  }
  return 0;
}
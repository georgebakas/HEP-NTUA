std::vector<TString> listFiles(const char *dirname="2016/Unfolding_PSWeights", const char *ext=".root")
{
  std::vector<TString> list_of_files;
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
     TSystemFile *file;
     TString fname;
     TIter next(files);
     while ((file=(TSystemFile*)next())) {
       fname = file->GetName();
       if (!file->IsDirectory() && fname.EndsWith(ext)) {
         //cout << fname.Data() << endl;
         list_of_files.push_back(fname.Data());
       }
     }
   }
   return list_of_files;
}

void test()
{
  std::vector<TString> files = listFiles();

  for (int i=0; i<files.size(); i++)
  {
    cout<<files[i].Data()<<endl;
  }
}

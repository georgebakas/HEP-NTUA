#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

////////////////////
//  constructors  //
////////////////////
summary_order::summary_order(){}

summary_order::summary_order(string _resultdirectory, summary_generic & _generic){
  Logger logger("summary_order::summary_order");
  logger << LOG_DEBUG << "called" << endl;

  // only one 'summary_generic'
  ygeneric = &_generic;
  // only one 'observable_set'
  osi = &ygeneric->oset;

  resultdirectory = _resultdirectory;
  contribution_file = vector<string> ();
  combination.resize(1);
  combination_type.resize(1, 0);
  accuracy_relative = 0.;
  accuracy_normalization = "";

  logger << LOG_DEBUG << "finished" << endl;
}




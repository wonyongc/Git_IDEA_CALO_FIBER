#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>
#include <algorithm>
#include "TMath.h"

double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity);
double FindSmallestInterval(std::vector<double>* vals, const double fraction, const bool verbosity);


double poissonf(double*x,double*par);

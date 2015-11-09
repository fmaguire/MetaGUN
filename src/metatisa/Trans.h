
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Methods for format transformation
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifndef TRNAS_H
#define TRNAS_H

#include <string>
#include <vector>
#include <assert.h>
#include <math.h>
#include "TypeDefBase.h"
using namespace std;

// convert letter format of DNA sequence to digital format
char Letter2Digital(const char res);
string Letter2Digital(const string& seq);

// convert digital format of DNA sequence to letter format
char Digital2Letter (const char res);
string Digital2Letter(const string& seq);

// transfer to opposite strand
char revs(const char res);
string revs(const string& s);

// transfermation between motif and digital
int motif2Digital(const string& seq);
string digital2Motif(const int index, const int len);

#endif

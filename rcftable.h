#ifndef _RCFTABLE_H
#define _RCFTABLE_H
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include "finfield.h"
using namespace std;

void getRCFs2(vector<int> &rcfs, vector<int> &rcftypes, finfield &F);
void getRCFs3(vector<int> &rcfs, vector<int> &rcftypes, finfield &F);

// we have that invtable[matindex(M)] = 0 if M is not invertible
// only okay if q^(d^2) < 2^31
void make_invtable(int d, finfield &F, int *invtable);

void make_rcftable(int d, finfield &F, int *invtable, int *rcftable, int *simtable);

void init_tables(int d, finfield &F, int **invtable, int **rcftable, int **simtable, 
				 vector<int> &rcfs, vector<int> &rcftypes);

void free_tables(int *invtable, int *rcftable, int *simtable);
#endif

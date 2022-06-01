/* File parsing and basic geometry operations */
#ifndef TMalign_basic_fun_h
#define TMalign_basic_fun_h 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <malloc.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <map>

#include "pstream.h" // For reading gzip and bz2 compressed files

using namespace std;


void PrintErrorAndQuit(const string sErrorString);

template <typename T> T getmin(const T &a, const T &b)
{
    return b<a?b:a;
}

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

//template <typename T> T getmin(const T &a, const T &b);
/*
template <typename T> inline T getmin(const T &a, const T &b)
{
    return b<a?b:a;
}
*/
//template <class A> void NewArray(A *** array, int Narray1, int Narray2);
//template <class A> void DeleteArray(A *** array, int Narray);

string AAmap(char A);

char AAmap(const string &AA);

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void split(const string &line, vector<string> &line_vec,
    const char delimiter=' ');

/* strip white space at the begining or end of string */
string Trim(const string &inputString);

size_t get_PDB_lines(const string filename,
    vector<vector<string> >&PDB_lines, vector<string> &chainID_list,
    vector<int> &mol_vec, const int ter_opt, const int infmt_opt,
    const string atom_opt, const int split_opt, const int het_opt);

/* read fasta file from filename. sequence is stored into FASTA_lines
 * while sequence name is stored into chainID_list.
 * if ter_opt >=1, only read the first sequence.
 * if ter_opt ==0, read all sequences.
 * if split_opt >=1 and ter_opt ==0, each sequence is a separate entry.
 * if split_opt ==0 and ter_opt ==0, all sequences are combined into one */
size_t get_FASTA_lines(const string filename,
    vector<vector<string> >&FASTA_lines, vector<string> &chainID_list,
    vector<int> &mol_vec, const int ter_opt=3, const int split_opt=0);

int read_PDB(const vector<string> &PDB_lines, double **a, char *seq,
    vector<string> &resi_vec, const int read_resi);

double dist(double x[3], double y[3]);
double dot(double *a, double *b);
void transform(double t[3], double u[3][3], double *x, double *x1);
void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3]);

/* read user specified pairwise alignment from 'fname_lign' to 'sequence'.
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void read_user_alignment(vector<string>&sequence, const string &fname_lign,
    const int i_opt);

/* read list of entries from 'name' to 'chain_list'.
 * dir_opt is the folder name (prefix).
 * suffix_opt is the file name extension (suffix_opt).
 * This function should only be called by main function, as it will
 * terminate a program if wrong alignment is given */
void file2chainlist(vector<string>&chain_list, const string &name,
    const string &dir_opt, const string &suffix_opt);

#endif

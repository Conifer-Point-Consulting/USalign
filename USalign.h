#ifndef USalign_h
#define USalign_h 1

#include <vector>
#include <string>

using namespace std;

struct TMatrix {
    double transform[3];
    double rotation_matrix[3][3];
};

int TMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int cp_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const string &dir1_opt, const string &dir2_opt, const int byresi_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const bool se_opt,
    TMatrix *tmatrix=nullptr
);

int MMalign(const string &xname, const string &yname,
    const string &fname_super, const string &fname_lign,
    const string &fname_matrix, vector<string> &sequence,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool d_opt, const bool full_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list);

int MMdock(const string &xname, const string &yname, const string &fname_super, 
    const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool u_opt, const bool d_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list);

int mTMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_matrix,
    vector<string> &sequence, double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    bool u_opt, const bool d_opt, const double TMcut,
    const int infmt_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, bool fast_opt,
    const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const int byresi_opt,
    const vector<string> &chain_list);

int main_original(int argc, char *argv[]);

#endif

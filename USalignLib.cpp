#include "USalignLib.h"

using namespace std;

int doAlignment(string pdb1, string pdb2, TMatrix &tmatrix) {
    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = pdb1;
    string yname       = pdb2;
    string fname_super = ""; // file name for superposed structure
    string fname_lign  = ""; // file name for user alignment
    string fname_matrix= ""; // file name for output matrix
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    int  i_opt = 0;     // 1 for -i, 3 for -I
    int  o_opt = 0;     // 1 for -o, 2 for -rasmol
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    bool   full_opt  = false;// do not show chain level alignment
    double TMcut     =-1;
    bool   se_opt    =false;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =2;     // END, or different chainID
    int    split_opt =2;     // split each chains
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    cp_opt    =0;     // do not check circular permutation
    int    mirror_opt=0;     // do not align mirror
    int    het_opt=0;        // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    int    byresi_opt=0;     // set -byresi to 0
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set


    if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    /* read initial alignment file from 'align.txt' */
    if (i_opt) read_user_alignment(sequence, fname_lign, i_opt);

    if (byresi_opt) i_opt=3;

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    if (dir1_opt.size()+dir_opt.size()==0) chain1_list.push_back(xname);
    else file2chainlist(chain1_list, xname, dir_opt+dir1_opt, suffix_opt);

    int i; 
    if (dir_opt.size())
        for (i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (dir2_opt.size()==0) chain2_list.push_back(yname);
    else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);

    /* real alignment. entry functions are MMalign_main and 
     * TMalign_main */
    TMalign(xname, yname, fname_super, fname_lign, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir_opt, dir1_opt, dir2_opt, byresi_opt,
        chain1_list, chain2_list, se_opt, &tmatrix);

    /* clean up */
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(sequence);

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    if (outfmt_opt<2) printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}

TMatrix getMatrix(string pdb_str1, string pdb_str2)
{
    TMatrix tmatrix;
    doAlignment(pdb_str1, pdb_str2, tmatrix); 
    return tmatrix;
}

#include "se.h"
#include "TMalign.h"

struct TMatrix {
    double transform[3];
    double rotation_matrix[3][3];

    void print_transform() {
        printf("Transform: %f\t%f\t%f\n", transform[0], transform[1], transform[2]);
    }
    void print_rotation() {
        printf("Rotation Matrix:\n");
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                printf("%f\t", rotation_matrix[i][j]);
            }
            printf("\n");
        }
    }
};

/* TMalign, RNAalign, CPalign, TMscore */
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
    const bool se_opt, TMatrix &tmatrix)
{
    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    char   *secx, *secy;       // for the secondary structure 
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=byresi_opt;  // whether to read residue index
    if (byresi_opt==0 && o_opt) read_resi=2;

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt1_opt, atom_opt, split_opt, het_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, 
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx,xa, xlen, secx,atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, split_opt,
                        het_opt);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (byresi_opt) extract_aln_from_resi(sequence,
                        seqx,seqy,resi_vec1,resi_vec2,byresi_opt);

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;

                    /* entry function for structure alignment */
                    if (cp_opt) CPalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);
                    else if (se_opt)
                    {
                        int *invmap = new int[ylen+1];
                        u0[0][0]=u0[1][1]=u0[2][2]=1;
                        u0[0][1]=         u0[0][2]=
                        u0[1][0]=         u0[1][2]=
                        u0[2][0]=         u0[2][1]=
                        t0[0]   =t0[1]   =t0[2]   =0;
                        se_main(
                            xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j], 
                            outfmt_opt, invmap);
                        if (outfmt_opt>=2) 
                        {
                            Liden=L_ali=0;
                            int r1,r2;
                            for (r2=0;r2<ylen;r2++)
                            {
                                r1=invmap[r2];
                                if (r1<0) continue;
                                L_ali+=1;
                                Liden+=(seqx[r1]==seqy[r2]);
                            }
                        }
                        delete [] invmap;
                    }
                    else TMalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);

                    /* print result */
                    // if (outfmt_opt==0) print_version();
                    if (cp_opt) output_cp(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        seqxA,seqyA,outfmt_opt);
                    output_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix+chainID_list1[chain_i]:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    
                    tmatrix.transform[0] = t0[0];
                    tmatrix.transform[1] = t0[1];
                    tmatrix.transform[2] = t0[2];
                    tmatrix.rotation_matrix[0][0] = u0[0][0];
                    tmatrix.rotation_matrix[0][1] = u0[0][1];
                    tmatrix.rotation_matrix[0][2] = u0[0][2];
                    tmatrix.rotation_matrix[1][0] = u0[1][0];
                    tmatrix.rotation_matrix[1][1] = u0[1][1];
                    tmatrix.rotation_matrix[1][2] = u0[1][2];
                    tmatrix.rotation_matrix[2][0] = u0[2][0];
                    tmatrix.rotation_matrix[2][1] = u0[2][1];
                    tmatrix.rotation_matrix[2][2] = u0[2][2];


                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    delete [] secy;
                    resi_vec2.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete [] seqx;
            delete [] secx;
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

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
        chain1_list, chain2_list, se_opt, tmatrix);

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

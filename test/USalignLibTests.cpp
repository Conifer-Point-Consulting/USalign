#include "../USalignLib.h"
#include "test_data.cpp"

/*
    Runs tests on:
    getMatrix()
*/


int main(int argc, char *argv[])
{
    TMatrix tmatrix1;
    TMatrix tmatrix2;
    TMatrix tmatrix3;
    string p1, p2;
    if (argc > 2) {
        p1 = argv[1];
        p2 = argv[2];
        tmatrix1 = getMatrix(p1, p2);

        printf("\n%s aligned to %s\n", p1.c_str(), p2.c_str());
        tmatrix1.print_transform();
        tmatrix1.print_rotation();
    } else {
        // run other tests:
        tmatrix1 = getMatrix(FOUR_RMG_BUILTIN, J3_SIRT1);     // pdb strings
        tmatrix2 = getMatrix("PDB1.pdb", "PDB2.pdb");     // pdb files
        tmatrix3 = getMatrix(FOUR_RMG_BUILTIN, "PDB1.pdb");     // pdb string aligned to pdb file

        printf("\n4rmg_builtin aligned to J3_SIRT1\n");
        tmatrix1.print_transform();
        tmatrix1.print_rotation();
        printf("\nPDB1.pdb aligned to PDB2.pdb\n");
        tmatrix2.print_transform();
        tmatrix2.print_rotation();
        printf("\n4rmg_builtin aligned to PDB1.pdb\n");
        tmatrix3.print_transform();
        tmatrix3.print_rotation();
    }

    return 0;
}

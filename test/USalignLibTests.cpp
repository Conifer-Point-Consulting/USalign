#include <string>
#include <stdio.h>
#include "../USalignLib.h"

/*
    Runs tests on:
    getMatrix()
*/


void print_transform(TMatrix &tmatrix) {
    printf("Transform: %f\t%f\t%f\n", tmatrix.transform[0], tmatrix.transform[1], tmatrix.transform[2]);
}

void print_rotation(TMatrix &tmatrix) {
    printf("Rotation Matrix:\n");
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            printf("%f\t", tmatrix.rotation_matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    TMatrix tmatrix1;
    TMatrix tmatrix2;
    TMatrix tmatrix3;
    std::string p1, p2;
    if (argc > 2) {
        p1 = argv[1];
        p2 = argv[2];
        tmatrix1 = getMatrix(p1, p2);

        printf("\n%s aligned to %s\n", p1.c_str(), p2.c_str());
        print_transform(tmatrix1);
        print_rotation(tmatrix1);
    } else {
        // run built-in test
        tmatrix1 = getMatrix("PDB1.pdb", "PDB2.pdb");
        printf("\nPDB1.pdb aligned to PDB2.pdb\n");
        print_transform(tmatrix1);
        print_rotation(tmatrix1);
    }

    return 0;
}

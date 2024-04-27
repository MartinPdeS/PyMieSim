extern "C" {
    void jdzo(int& nt, int* n, int* m, char* p, double* zo);
}

int main() {
    int nt = 1400; // Number of zeros
    int n[1400], m[1400];
    char p[1400 * 4];
    double zo[1400];

    // Call the Fortran subroutine
    jdzo(nt, n, m, p, zo);

}
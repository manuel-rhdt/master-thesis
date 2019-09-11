#include "Gillespie.h"
#include "numtools/numtools.h"

/*------------------------Globally defined functions-------------------------*/

void analyse(int, int, int, int);


void analyse(int run, int block, int n_steps_per_block, int step) {
    int i, ns;
    char filename[40];
    FILE *fp;

    ns = block * n_steps_per_block + step;

    for (i = 0; i < sys.Ncomp; i++) {

        sprintf(filename, "%s.%d.%s", sys.name, run, Xname[i]);
        fp = fopen(filename, "a");

        if (ns == 0) {
            fprintf(fp, "#time\t\t");
            fprintf(fp, "%s\t", Xname[i]);
            fprintf(fp, "step\n");
        }
        fprintf(fp, "%e\t", sys.tau_run + sys.tau_blk);
        fprintf(fp, "%d\t", X[i]);
        fprintf(fp, "%d\n", ns);

        fclose(fp);
    }
    return;
}

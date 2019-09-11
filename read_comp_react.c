#include "Gillespie.h"
#include "stats.h"

/*---------------------Globally defined functions----------------------------*/

void read_components();

void read_reactions();

void read_components() {
    int i, Ncomp;
    char filename[40];
    FILE *fp;

    sprintf(filename, "%s.components", sys.name);
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("could not open %s\n", filename);
        abort();
    } else {
        fscanf(fp, "%d\n", &Ncomp);
        if (Ncomp != sys.Ncomp) {
            printf("The number of components is %d\n", Ncomp);
            printf("The number of components should be %d\n", sys.Ncomp);
            abort();
        } else
            for (i = 0; i < sys.Ncomp; i++) fscanf(fp, "%d\t\t%s\n", &X[i], Xname[i]);
    }
    return;
}

void read_reactions() {
    int i, j, Nreact;
    char filename[40], dummy[40];
    FILE *fp;

    sprintf(filename, "%s.reactions", sys.name);
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("could not open %s\n", filename);
        abort();
    } else {
        fscanf(fp, "%d%*s\n", &Nreact);
        printf("Nreact is %d\n", Nreact);
        if (Nreact != sys.Nreact) {
            printf("The number of reactions is %d\n", Nreact);
            printf("The number of reactions should be %d\n", sys.Nreact);
            abort();
        } else {

            for (i = 0; i < sys.Nreact; i++) {

                fscanf(fp, "%lg %d %d %s\n", &R[i].k, &R[i].Nreact, &R[i].Nprod, &dummy);

                if (R[i].Nreact == 0)
                    fscanf(fp, "%s", &dummy);
                else
                    fscanf(fp, "%s %d", &dummy, &R[i].react[0].index);
                for (j = 1; j < R[i].Nreact; j++)
                    fscanf(fp, "%s %s %d", &dummy, &dummy, &R[i].react[j].index);

                fscanf(fp, "%s", &dummy);

                if (R[i].Nprod == 0)
                    fscanf(fp, "%s", &dummy);
                else
                    fscanf(fp, "%d %s %d\n", &R[i].prod[0].change,
                           &dummy, &R[i].prod[0].index);

                for (j = 1; j < R[i].Nprod; j++)
                    fscanf(fp, "%s %d %s %d\n", &dummy, &R[i].prod[j].change,
                           &dummy, &R[i].prod[j].index);
            }
        }
    }
    return;
}
	
		 

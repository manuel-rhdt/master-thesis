/* numtools.c
   Common tools for numerical code
   Sander Pronk, 1996,97,98
   Joel Wijngaarde, 1998
   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#ifdef MEM_CHECK
#include <dmalloc.h>
#endif

#define NR_END 1
#define FREE_ARG char*

#include "numtools.h"

char buffer[1024];

long int seed = -58947482;

static FILE *olog;
static boolean auto_flush_log;

/* local functions */

double h_norm(hist *h);

#ifdef __GNUC__

__inline__ double sqr(double x)
#else
double sqr(double x)
#endif
{
    return x * x;
}

/***************************** 
** memory allocation functions 
*/

void *malloc1D(int sizex, size_t sizedata) {
    char *data;

    data = (char *) malloc((size_t) (sizex * sizedata));
    if (!data) {
        log_error("malloc problem in malloc1D\n");
        exit(1);
    }

    return (void *) data;
}

void *realloc1D(void *ptr, int sizex, size_t sizedata) {
    char *data;

    data = (char *) realloc(ptr, (size_t) (sizex * sizedata));
    if (!data) {
        log_error("realloc problem is realloc1D\n");
        exit(1);
    }

    return (void *) data;

}

void **malloc2D(int sizex, int sizey, size_t sizedata) {
    char **data;
    int i;                     /* memory location indices */

    data = (char **) malloc1D(sizex, sizeof(void *));
    data[0] = (char *) malloc1D(sizex * sizey, sizedata);

    /* put all pointers in the right spot */
    for (i = 1; i < sizex; i++)
        data[i] = data[i - 1] + sizey * sizedata;

    /* return the pointer */
    return (void **) data;
}


void ***malloc3D(int sizex, int sizey, int sizez, size_t sizedata) {
    char ***data;
    int i, j;

    data = (char ***) malloc2D(sizex, sizey, sizeof(char *));
    data[0][0] = (char *) malloc1D(sizex * sizey * sizez, sizedata);

    /* put all pointers in the right spot */

    for (j = 1; j < sizey; j++) {
        data[0][j] = data[0][j - 1] + sizez * sizedata;
    }
    for (i = 1; i < sizex; i++) {
        data[i][0] = data[i - 1][0] + sizez * sizey * sizedata;
        for (j = 1; j < sizey; j++) {
            data[i][j] = data[i][j - 1] + sizez * sizedata;
        }
    }

    /* return the pointer */
    return (void ***) data;
}


double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double ***t;

    /* allocate pointers to pointers to rows */
    t = (double ***) malloc((size_t) ((nrow + NR_END) * sizeof(double **)));
    if (!t) printf("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double **) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double *)));
    if (!t[nrl]) printf("allocation failure 2 in df3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double *) malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(double)));
    if (!t[nrl][ncl]) printf("allocation failure 3 in df3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

int ***i3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    int ***t;

    /* allocate pointers to pointers to rows */
    t = (int ***) malloc((size_t) ((nrow + NR_END) * sizeof(int **)));
    if (!t) printf("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (int **) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int *)));
    if (!t[nrl]) printf("allocation failure 2 in df3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (int *) malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(int)));
    if (!t[nrl][ncl]) printf("allocation failure 3 in df3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

void free1D(void *array) {
    /* NULL check not neccesary, 
       because ANSI does it for us */
    free(array);
}

void free2D(void **array) {
    if (array) {
        free1D(array[0]);
        free1D(array);
    }
}

void free3D(void ***array) {
    if (array) {
        free1D(array[0][0]);
        free1D(array[0]);
        free1D(array);
    }
}

void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
                   int ndl, int ndh)
/* free a double df3tensor allocated by df3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
}

void free_i3tensor(int ***t, int nrl, int nrh, int ncl, int nch,
                   int ndl, int ndh)
/* free a double df3tensor allocated by df3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
}


/***************************
** data read/write functions 
*/

FILE *open_data(const char *filename) {
    FILE *theFile;
    /*printf("Opening: %s\n",filename);
    fflush(stdout);*/
    theFile = fopen(filename, "r");

    return theFile;
}

int read_data_int(FILE *theFile) {
    int read;

    do {
        if (feof(theFile)) {
            fprintf(stderr, "ERROR: read_data_int: End of file reached\n");
            return 0;
        }
        if (ferror(theFile)) {
            fprintf(stderr, "ERROR: read_data_int: File error\n");
            perror("read_data_int");
            return 0;
        }
        fscanf(theFile, "%s", buffer);
    } while (!(is_numeric(buffer)));

    sscanf(buffer, "%d", &read);
    return read;
}

int is_numeric(char *buf) {
    double dum;

    if (sscanf(buf, "%lg", &dum) == 0)
        return FALSE;
    else
        return TRUE;
}

double read_data_dbl(FILE *theFile) {
    double read;

    do {
        if (feof(theFile)) {
            fprintf(stderr, "ERROR: read_data_dbl: End of file reached\n");
            return 0;
        }
        if (ferror(theFile)) {
            fprintf(stderr, "ERROR: read_data_dbl: File error\n");
            perror("read_data_dbl");
            return 0;
        }
        fscanf(theFile, "%s", buffer);
    } while (!(is_numeric(buffer)));

    sscanf(buffer, "%lg", &read);
    return read;
}

void read_data_string(FILE *theFile, char *string) {
    int dum;

    do {
        dum = getc(theFile);
        if (dum <= 0) {
            if (feof(theFile)) {
                fprintf(stderr, "ERROR: read_string: end of file reached\n");
                return;
            } else {
                fprintf(stderr,
                        "ERROR: read_string: File error before string\n");
                perror("read_data_string");
                return;
            }
        }
    } while (dum != ':');
    if (!fscanf(theFile, "%s", string)) {
        fprintf(stderr, "File error\n");
        perror("read_data_string");
        return;
    }
}

void close_data(FILE *theFile) {
    fclose(theFile);
}

/******************** 
** log-file functions 
*/

static struct tms tmsstart, tmsend, tmsprev;/* stores process time */
static clock_t start, end, prev;

int log_init(const char *name, boolean auto_flush) {
    char filename[FILENAME_SIZE];
    int pid;
    time_t tp;

    if (!strstr(name, ".log")) {
        strcpy(filename, name);
        strcat(filename, ".log");
    }
    olog = fopen(filename, "w");

    if (!olog)
        return 1;

    if (gethostname(filename, FILENAME_SIZE))
        perror("log_init");

    pid = (int) getpid();

    time(&tp);

    start = times(&tmsstart);
    prev = times(&tmsprev);

    auto_flush_log = auto_flush;

    printf("Log file of run: %s\n", name);
    printf("Creation date: %s", ctime(&tp));
    printf("Hostname: %s\n", filename);
    printf("Process ID: %d\n", pid);

    fprintf(olog, "# Log file of run: %s\n", name);
    fprintf(olog, "# Creation date: %s", ctime(&tp));
    fprintf(olog, "# Hostname: %s\n", filename);
    fprintf(olog, "# Process ID: %d\n\n", pid);

    if (auto_flush_log)
        fflush(olog);

    return 0;
}

void log_exit(void) {
    long clktck;
    time_t tp;

    time(&tp);

    clktck = sysconf(_SC_CLK_TCK);
    end = times(&tmsend);

    printf("\nTime values: %.2fuser %.2fsys %.2freal",
           (tmsend.tms_utime - tmsstart.tms_utime) / (double) clktck,
           (tmsend.tms_stime - tmsstart.tms_stime) / (double) clktck,
           (end - start) / (double) clktck);
    printf("\nLog closed at: %s", ctime(&tp));

    fprintf(olog, "\n# Time values: %.2fuser %.2fsys %.2freal",
            (tmsend.tms_utime - tmsstart.tms_utime) / (double) clktck,
            (tmsend.tms_stime - tmsstart.tms_stime) / (double) clktck,
            (end - start) / (double) clktck);
    fprintf(olog, "\n# Log closed at: %s", ctime(&tp));
    fclose(olog);
}

void log_write(const char *output, ...) {
    va_list argp;

    va_start(argp, output);
    vfprintf(olog, output, argp);
    va_end(argp);
    fprintf(olog, "\n");
    if (auto_flush_log)
        fflush(olog);
}

void log_flush(void) {
    fflush(olog);
}

void log_error(const char *output, ...) {
    va_list argp;

    fprintf(olog, "** ERROR **: ");
    va_start(argp, output);
    vfprintf(olog, output, argp);
    va_end(argp);
    fprintf(olog, "\n");
    if (auto_flush_log)
        fflush(olog);
}

void log_warning(const char *output, ...) {
    va_list argp;

    fprintf(olog, "** WARNING **: ");
    va_start(argp, output);
    vfprintf(olog, output, argp);
    va_end(argp);
    fprintf(olog, "\n");
    if (auto_flush_log)
        fflush(olog);
}

void log_time(int mask, const char *output, ...) {
    va_list argp;
    long clktck;

    clktck = sysconf(_SC_CLK_TCK);
    end = times(&tmsend);

    fprintf(olog, "#");

    if (LOG_USER_TIME & mask)
        fprintf(olog, " %.2fuser",
                (tmsend.tms_utime - tmsstart.tms_utime) / (double) clktck);
    if (LOG_SYSTEM_TIME & mask)
        fprintf(olog, " %.2fsys",
                (tmsend.tms_stime - tmsstart.tms_stime) / (double) clktck);
    if (LOG_REAL_TIME & mask)
        fprintf(olog, " %.2freal", (end - start) / (double) clktck);

    fprintf(olog, ": ");

    va_start(argp, output);
    vfprintf(olog, output, argp);
    va_end(argp);

    fprintf(olog, "\n");

    if (auto_flush_log)
        fflush(olog);
}

void log_delta_time(int mask, const char *output, ...) {
    va_list argp;
    long clktck;

    clktck = sysconf(_SC_CLK_TCK);
    end = times(&tmsend);

    fprintf(olog, "# ");

    va_start(argp, output);
    vfprintf(olog, output, argp);
    printf(output, argp);
    va_end(argp);

    if (LOG_USER_TIME & mask) {
        fprintf(olog, " %.2fuser",
                (tmsend.tms_utime - tmsprev.tms_utime) / (double) clktck);
        printf(" %.2fuser",
               (tmsend.tms_utime - tmsprev.tms_utime) / (double) clktck);
    }
    if (LOG_SYSTEM_TIME & mask) {
        fprintf(olog, " %.2fsys",
                (tmsend.tms_stime - tmsprev.tms_stime) / (double) clktck);
        printf(" %.2fsys",
               (tmsend.tms_stime - tmsprev.tms_stime) / (double) clktck);
    }
    if (LOG_REAL_TIME & mask) {
        fprintf(olog, " %.2freal", (end - prev) / (double) clktck);
        printf(" %.2freal", (end - prev) / (double) clktck);
    }

    fprintf(olog, "\n");
    printf("\n");

    if (auto_flush_log)
        fflush(olog);

    prev = times(&tmsprev);

}


/************************************* 
** one-dimensional histogram functions 
*/

hist *h_create(int nbin, double min, double max) {
    hist *newh;
    int i;

    if (min > max)
        return NULL;
    if (nbin < 2)
        return NULL;

    newh = (hist *) malloc(sizeof(hist));
    if (newh == NULL)
        return NULL;
    newh->val = (int *) malloc(sizeof(int) * nbin);
    if (newh->val == NULL) {
        free(newh);
        return NULL;
    }
    newh->fact = (double *) malloc(sizeof(double) * nbin);
    if (newh->fact == NULL) {
        free(newh->val);
        free(newh);
        return NULL;
    }
    newh->nbin = nbin;
    newh->min = min;
    newh->width = max - min;
    newh->total = 0;

    for (i = 0; i < nbin; i++) {
        newh->val[i] = 0;
        newh->fact[i] = 1.;
    }
    return newh;
}

void h_destroy(hist *h) {
    free(h->val);
    free(h->fact);
    free(h);
}

int h_bin(hist *h, double x) {
    return (h->nbin * ((x - h->min) / h->width));
}

double h_x(hist *h, int bin) {
    return (h->width * (double) bin / (double) h->nbin) + h->min;
}


void h_add(hist *h, double x) {
    int bin = h_bin(h, x);

    if ((bin >= 0) && (bin < h->nbin)) {
        h->val[bin]++;
        h->total++;
    }
}

void h_set_count(hist *h, double x, int n) {
    int bin = h_bin(h, x);

    if ((bin > 0) && (bin < h->nbin)) {
        h->total += (n - (h->val[bin]));
        h->val[bin] = n;
    }
}

int h_get_count(hist *h, double x) {
    int bin = h_bin(h, x);

    if ((bin > 0) && (bin < h->nbin))
        return h->val[bin];
    else
        return 0;
}

double h_norm(hist *h) {
    return 1. / ((h->width / h->nbin) * (double) h->total);
}

double h_average(hist *h) {
    int i;
    double tmp = 0;

    for (i = 0; i < h->nbin; i++) {
        tmp += h->val[i] * h_x(h, i);
    }
    return tmp / (double) h->total;
}

void h_reset(hist *h) {
    int i;

    h->total = 0;

    for (i = 0; i < h->nbin; i++) {
        h->val[i] = 0;
        h->fact[i] = 0;
    }
}

double h_bin_width(hist *h) {
    return h->width / h->nbin;
}

int h_n_binned(hist *h) {
    return h->total;
}

int h_n_bins(hist *h) {
    return h->nbin;
}

void h_set_fact(hist *h, int bin, double fact) {
    if ((bin >= 0) && (bin < h->nbin))
        h->fact[bin] = fact;
}

double h_get_fact(hist *h, int bin) {
    if ((bin >= 0) && (bin < h->nbin))
        return h->fact[bin];
    return 0;
}

double h_get_norm_count(hist *h, double x) {
    int bin = h_bin(h, x);
    double norm = h_norm(h);

    if ((bin > 0) && (bin < h->nbin))
        return norm * h->val[bin];
    else
        return 0;
}

int h_write(hist *h, const char *filename) {
    FILE *out;
    int i;

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nbin; i++) {
        if (!fprintf(out, "%g\t%g\n", h_x(h, i), h->fact[i] * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
    }
    fclose(out);
    return TRUE;
}

int h_norm_write(hist *h, const char *filename) {
    FILE *out;
    int i;
    double norm = h_norm(h);

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nbin; i++) {
        if (!fprintf(out, "%g\t%g\n", h_x(h, i), h->fact[i]
                                                 * norm * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
    }
    fclose(out);
    return TRUE;
}

int h_plot_write(hist *h, const char *filename) {
    FILE *out;
    int i;

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nbin; i++) {
        if (!fprintf(out, "%g\t%g\n", h_x(h, i), h->fact[i] * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
        if (!fprintf(out, "%g\t%g\n", h_x(h, i + 1), h->fact[i] * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
    }
    fclose(out);
    return TRUE;
}

int h_plot_norm_write(hist *h, const char *filename) {
    FILE *out;
    int i;
    double norm = h_norm(h);

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nbin; i++) {
        if (!fprintf(out, "%g\t%g\n", h_x(h, i), h->fact[i] * norm
                                                 * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
        if (!fprintf(out, "%g\t%g\n", h_x(h, i + 1), h->fact[i] * norm
                                                     * (double) h->val[i])) {
            fclose(out);
            return FALSE;
        }
    }
    fclose(out);
    return TRUE;
}

/************************************* 
** two-dimensional histogram functions 
*/

hist2D *h2D_create(int nxbin, int nybin,
                   double xmin, double xmax, double ymin, double ymax) {
    hist2D *newh;
    int i, j;

    if (xmin > xmax || ymin > ymax)
        return NULL;
    if (nxbin < 2 || nybin < 2)
        return NULL;

    newh = (hist2D *) malloc(sizeof(hist2D));
    if (newh == NULL)
        return NULL;
    newh->val = (int **) malloc2D(nxbin, nybin, sizeof(int));
    if (newh->val == NULL) {
        free(newh);
        return NULL;
    }
    newh->fact = (double **) malloc2D(nxbin, nybin, sizeof(double));
    if (newh->fact == NULL) {
        free2D((void **) (newh->val));
        free(newh);
        return NULL;
    }
    newh->nxbin = nxbin;
    newh->nybin = nybin;
    newh->xmin = xmin;
    newh->ymin = ymin;
    newh->xwidth = xmax - xmin;
    newh->ywidth = ymax - ymin;
    newh->total = 0;

    for (i = 0; i < nxbin; i++) {
        for (j = 0; j < nybin; j++) {
            newh->val[i][j] = 0;
            newh->fact[i][j] = 1.0;
        }
    }
    return newh;
}

void h2D_destroy(hist2D *h) {
    free2D((void **) (h->val));
    free2D((void **) (h->fact));
    free(h);
}

int h2D_xbin(hist2D *h, double x) {
    return (h->nxbin * ((x - h->xmin) / h->xwidth));
}

int h2D_ybin(hist2D *h, double y) {
    return (h->nybin * ((y - h->ymin) / h->ywidth));
}

double h2D_x(hist2D *h, int xbin) {
    return (h->xwidth * (double) xbin / (double) h->nxbin) + h->xmin;
}

double h2D_y(hist2D *h, int ybin) {
    return (h->ywidth * (double) ybin / (double) h->nybin) + h->ymin;
}

void h2D_add(hist2D *h, double x, double y) {
    int xbin = h2D_xbin(h, x);
    int ybin = h2D_ybin(h, y);

    if ((xbin >= 0) && (xbin < h->nxbin)) {
        if ((ybin >= 0) && (ybin < h->nybin)) {
            h->val[xbin][ybin]++;
            h->total++;
        }
    }
}

void h2D_set_count(hist2D *h, double x, double y, int n) {
    int xbin = h2D_xbin(h, x);
    int ybin = h2D_ybin(h, y);

    if ((xbin >= 0) && (xbin < h->nxbin)) {
        if ((ybin >= 0) && (ybin < h->nybin)) {
            h->total += (n - (h->val[xbin][ybin]));
            h->val[xbin][ybin] = n;
        }
    }
}

int h2D_get_count(hist2D *h, double x, double y) {
    int xbin = h2D_xbin(h, x);
    int ybin = h2D_ybin(h, y);

    if ((xbin >= 0) && (xbin < h->nxbin)) {
        if ((ybin >= 0) && (ybin < h->nybin)) {
            return h->val[xbin][ybin];
        }
    }
    return 0;
}

double h2D_norm(hist2D *h) {
    return 1.0 / ((h->xwidth / h->nxbin) * (h->ywidth / h->nybin) * h->total);
}

double h2D_average() {
    /* not implemented yet */
    /* maybe two funtions h2D_average_x & h2D_average_y */
    /* width an x cq. y value as argument */
    return 0;
}

void h2D_reset(hist2D *h) {
    int i, j;

    h->total = 0;

    for (i = 0; i < h->nxbin; i++) {
        for (j = 0; j < h->nybin; j++) {
            h->val[i][j] = 0;
            h->fact[i][j] = 0;
        }
    }
}

double h2D_xbin_width(hist2D *h) {
    return h->xwidth / h->nxbin;
}

double h2D_ybin_width(hist2D *h) {
    return h->ywidth / h->nybin;
}

int h2D_n_binned(hist2D *h) {
    return h->total;
}

int h2D_nx_bins(hist2D *h) {
    return h->nxbin;
}

int h2D_ny_bins(hist2D *h) {
    return h->nybin;
}

void h2D_set_fact() {
    /* not implemented yet */
}

double h2D_get_fact(hist2D *h, int xbin, int ybin) {
    if ((xbin >= 0) && (xbin < h->nxbin)) {
        if ((ybin >= 0) && (ybin < h->nybin)) {
            return h->fact[xbin][ybin];
        }
    }
    return 0;
}

double h2D_get_norm_count(hist2D *h, double x, double y) {
    int xbin = h2D_xbin(h, x);
    int ybin = h2D_ybin(h, y);
    double norm = h2D_norm(h);

    if ((xbin >= 0) && (xbin < h->nxbin)) {
        if ((ybin >= 0) && (ybin < h->nybin)) {
            return norm * h->val[xbin][ybin];
        }
    }

    return 0;
}

int h2D_write(hist2D *h, char seperator, const char *filename) {
    FILE *out;
    int i, j;

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nxbin; i++) {
        if (seperator == TRUE)
            if (!fprintf(out, "\n")) {
                fclose(out);
                return FALSE;
            }
        for (j = 0; j < h->nybin; j++) {
            if (!fprintf(out, "%g\t%g\t%g\n", h2D_x(h, i), h2D_y(h, j),
                         h->fact[i][j] * (double) h->val[i][j])) {
                fclose(out);
                return FALSE;
            }
        }
    }
    fclose(out);
    return TRUE;
}

int h2D_norm_write(hist2D *h, char seperator, const char *filename) {
    FILE *out;
    int i, j;
    double norm = h2D_norm(h);

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nxbin; i++) {
        if (seperator == TRUE)
            if (!fprintf(out, "\n")) {
                fclose(out);
                return FALSE;
            }
        for (j = 0; j < h->nybin; j++) {
            if (!fprintf(out, "%g\t%g\t%g\n", h2D_x(h, i), h2D_y(h, j),
                         h->fact[i][j] * norm * (double) h->val[i][j])) {
                fclose(out);
                return FALSE;
            }
        }
    }
    fclose(out);
    return TRUE;
}

int h2D_plot_write(hist2D *h, const char *filename)
/* not implemented yet */
{
    FILE *out;
    int i, j;

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nxbin; i++) {
        for (j = 0; j < h->nybin; j++) {
            /* not implemented yet */
        }
    }
    fclose(out);
    return TRUE;
}

int h2D_plot_norm_write(hist2D *h, const char *filename)
/* not implemented yet */
{
    FILE *out;
    int i, j;

    out = fopen(filename, "w");
    if (!out)
        return FALSE;

    for (i = 0; i < h->nxbin; i++) {
        for (j = 0; j < h->nybin; j++) {
            /* not implemented yet */
        }
    }
    fclose(out);
    return TRUE;
}

/*****************************
** create seed
*/

#define MSEED 161803398

long create_seed(void) {
    time_t tp;

    seed = -(long int) time(&tp);

#ifdef __alpha__        /* correct for 64bit time values */
    seed %= MSEED;
#endif

    return seed;
}

/***************************** 
** code from numerical recipes 
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum) {
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


double factrl(int n) {
    static int ntop = 6;
    static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0};
    int j;

    if (n < 0) perror("Negative factorial in routine factrl");
    if (n > 32) return exp(gammln(n + 1.0));
    while (ntop < n) {
        j = ntop++;
        a[ntop] = a[j] * ntop;
    }
    return a[n];
}

double gammln(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}


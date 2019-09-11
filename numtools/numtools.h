/* numtools.h
   Common tools for numerical code
   Sander Pronk, 1996,97,98
   Joel Wijngaarde, 1998
   */

#ifndef _NUMTOOLS_H_
#define _NUMTOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define DEG2RAD (PI/180.)

#define FILENAME_SIZE 256

#define LOG_USER_TIME    (1<<1)
#define LOG_SYSTEM_TIME  (1<<2)
#define LOG_REAL_TIME    (1<<3)
#define NO_TIME          (1<<4)

typedef enum {
    false = 0, true = 1
} boolean;

/*static double sqtmp;*/

extern long int seed;

/*#define sqr(x) (sqtmp=(x), sqtmp*sqtmp)*/

#ifdef __GNUC__

__inline__ double sqr(double x);

#else
double sqr(double x);
#endif

void *malloc1D(int sizex, size_t sizedata);

/* allocates a dynamic array of size sizex with data
   size sizedata. If the log is open any errors will be logged.

   sizex = number of elements
   sizedata = size of element
   returns: pointer to allocated data  */

void *realloc1D(void *ptr, int sizex, size_t sizedata);

/* reallocates a dynamic array with a new size sizex and data
   size sizedata. If the log is open any errors will be logged

   ptr = dynamically alloced array
   sizex = number of elements
   sizedata = size of element
   returns: pointer to realloced data */

void **malloc2D(int sizex, int sizey, size_t sizedata);

/* allocates a 2D dynamic array of size sizex and sizey, with
   data size sizedata. If the log is open any errors will be logged.

   sizex = number of elements in first subscript
   sizey = number of elements in second subscript
   sizedata = size of element
   returns: pointer to allocated data  */

void ***malloc3D(int sizex, int sizey, int sizez, size_t sizedata);

/* allocates a 3D dynamic array cq matrix of size (sizex X sizey X sizez),
   with size datatype.  If the log is open any errors will be logged.

   sizex = number of elements in first subscript
   sizey = number of elements in second subscript
   sizez = number of elements in third subscript
   sizedata = size of element
   returns: pointer to allocated data  */

double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

/*From numerical recipes. A 3D tensor, which ranges from nrl to nrh,
  ncl, nch, ndl, ndh.*/

int ***i3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

/*From numerical recipes. A 3D tensor, which ranges from nrl to nrh,
  ncl, nch, ndl, ndh.*/

void free1D(void *array);

/* frees memory held by 1d array.

   array = 1D array to free */

void free2D(void **array);

/* frees memory held by 2d array
 
   array = 2D array to free */

void free3D(void ***array);

/* frees memory held by 3D array/matrix

   array = 3D array to free */

void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
                   int ndl, int ndh);

/* frees memory from a 3D tensor.
   This routine is taken from numerical recipes.*/

void free_i3tensor(int ***t, int nrl, int nrh, int ncl, int nch,
                   int ndl, int ndh);

/* frees memory from a 3D tensor.
   This routine is taken from numerical recipes.*/

int is_numeric(char *buf);

/* checks whether string buf is numeric, returns TRUE when it is*/

FILE *open_data(const char *filename);

/* opens a data file named filename and returns the file pointer */

void close_data(FILE *theFile);

/* closes the data file */

int read_data_int(FILE *theFile);

/* reads the next integer from a line in the datafile. It
   skips any non-numbers in the file and reads the first number
   it can read */

double read_data_dbl(FILE *theFile);

/* reads a double from a line in the datafile. See readDataInt for further
   explanation */

void read_data_string(FILE *theFile, char *string);
/* reads a string */


/* log tools */

int log_init(const char *name, boolean auto_flush);

/* opens a log file, with the name "name.log", returns 0 if
   everything went OK, nonzero if an error occurred 
   if auto_flush is true, then the log file buffer will be flushed
   after every log_write*/

void log_write(const char *output, ...);

/* writes to log file (as in printf) */

void log_error(const char *output, ...);

/* writes to log file with prepended message ERROR: */

void log_warning(const char *output, ...);

/* writes to log file with prepended message WARNING: */

void log_time(int mask, const char *output, ...);

/* writes intermediate timing values, 'or' the defines declared above
   to determine the output form and the type of time that
   should be logged */

void log_delta_time(int mask, const char *output, ...);

/* writes "delta" timing values, 'or' the defines declared above
   to determine the output form and the type of time that
   should be logged */

void log_flush(void);

/* flushes the log buffer (Note that with log_error and log_warning
 the log is automatically flushed ) */

void log_exit(void);
/* closes the log */

/* histogram functions */

typedef struct {
    int *val;
    double *fact;
    int nbin;
    double min;
    double width;
    int total;
} hist;

typedef struct {
    int **val;
    double **fact;
    int nxbin;
    int nybin;
    double xmin;
    double ymin;
    double xwidth;
    double ywidth;
    int total;
} hist2D;

hist *h_create(int nbin, double min, double max);

/* create & initialize histogram with nbin bins, between min and max
   and return pointer to it */

void h_destroy(hist *h);

/* destroy hist h */

void h_add(hist *h, double x);

/* add one count to the histogram */

void h_set_count(hist *h, double x, int n);

/* set count value of bin at value x at n */

int h_get_count(hist *h, double x);

/* get count of bin at value x */

double h_get_norm_count(hist *h, double x);

/* get normalized count (does not work with multipilcation factors)
   of bin at value x */

double h_x(hist *h, int bin);

/* get x-value corresponding to bin nr. */

int h_bin(hist *h, double x);

/* get bin nr. corresponding to x-value */

double h_norm(hist *h);

/* get the normalization factor for the histogram, based on
   the number of binned values and bin width/number */

double h_average(hist *h);

/* calculate the plain average of the given histogram */

void h_reset(hist *h);

/* reset/clear the histogram */

double h_bin_width(hist *h);

/* return the bin width */

int h_n_binned(hist *h);

/* return the number of binned values */

int h_n_bins(hist *h);

/* return the number of bins */

void h_set_fact(hist *h, int bin, double fact);

/* set multiplication factor for bin bin */

double h_get_fact(hist *h, int bin);

/* get multiplication factor for bin bin */

int h_write(hist *h, const char *filename);

/* write the histogram values to file filename, returning TRUE if
   everything went OK, false otherwise */

int h_norm_write(hist *h, const char *filename);

/* write the normalized histogram values to filename */

int h_plot_write(hist *h, const char *filename);

/* write the histogram to file, suitable for plotting */

int h_plot_norm_write(hist *h, const char *filename);

/* write the normalized histogram to file, suitable for plotting */

hist2D *h2D_create(int nxbin, int nybin,
                   double xmin, double xmax, double ymin, double ymax);

/* create & initialize histogram with nxbin bins in the x-direction and nybin
   bins in the y-direction, between the min and max and return pointer to it */

void h2D_destroy(hist2D *h);

/* destroy hist2D h */

int h2D_xbin(hist2D *h, double x);

/* get x-bin nr. corresponding to x-value */

int h2D_ybin(hist2D *h, double y);

/* get y-bin nr. corresponding to y-value */

double h2D_x(hist2D *h, int xbin);

/* get x-value corresponding to xbin nr. */

double h2D_y(hist2D *h, int ybin);

/* get y-value corresponding to ybin nr. */

void h2D_add(hist2D *h, double x, double y);

/* add one count to the histogram */

void h2D_set_count(hist2D *h, double x, double y, int n);

/* set count value of bin at value x,y at n */

int h2D_get_count(hist2D *h, double x, double y);

/* get count of bin at value x,y */

double h2D_norm(hist2D *h);

/* get the normalization factor for the histogram, based on
   the number of binned values and bin width/number */

double h2D_average();

/* calculate the plain average of the given histogram */

void h2D_reset(hist2D *h);

/* reset/clear the histogram */

double h2D_xbin_width(hist2D *h);

/* return the xbin width */

double h2D_ybin_width(hist2D *h);

/* return the ybin width */

int h2D_n_binned(hist2D *h);

/* return the number of binned values */

int h2D_nx_bins(hist2D *h);

/* return the number of bins in x-direction */

int h2D_ny_bins(hist2D *h);

/* return the number of bins in x-direction */

void h2D_set_fact();

/* set multiplication factor for bin bin */

double h2D_get_fact(hist2D *h, int xbin, int ybin);

/* get multiplication factor for bin bin */

double h2D_get_norm_count(hist2D *h, double x, double y);

int h2D_write(hist2D *h, char seperator, const char *filename);

/* write the histogram values to file filename, returning TRUE if
   everything went OK, false otherwise */

int h2D_norm_write(hist2D *h, char seperator, const char *filename);

/* write the normalized histogram values to filename */

int h2D_plot_write(hist2D *h, const char *filename);

/* write the histogram to file, suitable for plotting */

int h2D_plot_norm_write(hist2D *h, const char *filename);

/* write the normalized histogram to file, suitable for plotting */

long create_seed(void);

/* create seed with the help of the clock */

double ran1(long *idum);

extern int iff, ma[56], inext, inextp, mj, mz, mbig;

extern double ran3(void);

double factrl(int n);

double gammln(double xx);

/* from numerical recipes */

/*double sqr(double x);*/


#endif

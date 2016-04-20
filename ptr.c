#include <math.h>
#include <netpbm/pam.h>

// https://engineering.purdue.edu/kak/computervision/ECE661/HOUGH.C

#define TAN_TOO_BIG 263.0       /* tan so big that angle is PI/2 */
#define PI 3.1415926536

static void houghtransform (int * edgedata, int * hough, int w, int h)
  {
    int rhoWid = w;
    int thetaHt = h;
    int rhoWidM1 = rhoWid - 1;
    int border = 1;
    /* initialize accumulator bins */
    for (int i = 0; i < thetaHt; i++)
      for (int j = 0; j < rhoWid; j++)
        hough [i * w + j] = 0;

    double rhoNorm = rhoWidM1 / sqrt ((double) (w * w) + (double) h * h);

    double tantable [h];
    for (int i = 0; i < h; i ++)
      {
        double theta = (double) i * PI / (double) thetaHt;
        double tanTheta = tan (theta);
        tantable [i] = tanTheta;
      }

    int yEnd = h - border;
    int xEnd = w - border;
    for (int y = border + 1; y < yEnd; y ++)
      {
        for (int x = border + 1; x < xEnd; x ++)
          {
            //if (imgIn [y] [x] == ON)
            if (edgedata [y * w + x])
              {
                for (int i = 0; i < thetaHt; i ++)
                  {
                    //double theta = (double) i * PI / (double) thetaHt;
                    //double tanTheta = tan (theta);
                    double tanTheta = tantable [i];
                    double rho;
                    if (tanTheta > TAN_TOO_BIG)
                      {
                        /*printf ("too big: x = %d, i = %d, j = %d\n", x, i, j);*/
                        rho = (double) x;
                      }
                    else
                      {
                        double denom = tanTheta * tanTheta + 1.0;
                        double y1 = ((double) y - (double) x * tanTheta) / denom;
                        double x1 = ((double) x * tanTheta * tanTheta - (double) y * tanTheta) / denom;
                        rho = sqrt (x1 * x1 + y1 * y1);
                      }
                    int j = (long) (rho * rhoNorm + 0.5);
                    //if (hough [i * w + j] < 255)
                      hough [i * w + j] ++;
                  }
              }
          }
      }
  }

static void dumpImage (int * img, int h, int w, char * name, int n)
  {
#if 0
    for (int r = 0; r < h; r ++)
      {
        for (int c = 0; c < w; c ++)
          printf (" %3d", img [r] [c]);
        printf ("\n");
      }
#endif
    char fname [132];
    sprintf (fname, "process/%s%05d.pam", name, n);
    FILE * fout = fopen (fname, "w");
    struct pam outpam;
    memset (& outpam, 0, sizeof (outpam));
    outpam.size = w * h;
    outpam.len = w * h;
    outpam.file = fout;
    outpam.format = PGM_FORMAT;
    //outpam.plainformat = -1;
    outpam.plainformat = 0;
    outpam.height = h;
    outpam.width = w;
    outpam.depth = 1;
    outpam.maxval = 255;
    outpam.bytes_per_sample = 1;
    outpam.allocation_depth = 0; 
    strcpy (outpam.tuple_type, "GRAYSCALE");

    pnm_writepaminit (& outpam);

    tuple * tuplerow;
    tuplerow = pnm_allocpamrow (& outpam);

    for (int row = 0; row < outpam.height; row ++)
      {
        unsigned int column;
        //pnm_readpamrow (& outpam, tuplerow);
        for (int column = 0; column < outpam.width; column ++)
          {
            //tuplerow [column] [0] = img [row] [column];
            tuplerow [column] [0] = img [row * w + column];
          }
        pnm_writepamrow (& outpam, tuplerow);
      }
    pnm_freepamrow (tuplerow);
    fclose (fout);
  }

// http://www.cglabprograms.com/2008/10/dda-line-drawing-algorithm.html
static void drawline (int * img, int h, int w, int x1, int y1, int x2, int y2, int color)
  {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int s;
    if (abs (dx) > abs (dy))
      s = abs (dx);
    else
      s = abs (dy);
 
    float xi = dx / (float) s;
    float yi = dy / (float) s;
 
    float x = x1;
    float y = y1;
 
    //img [y1 * w + x1] = color;
 
    for (int m = 0; m < s; m ++)
      {
        x += xi;
        y += yi;
        //img [(int) y * w + (int) x] = color;
      }
  }

static int nwrite = 1;

static int image (int n)
  {

    struct pam inpam;
    tuple * tuplerow;
    unsigned int row;

    char fname [132];
    sprintf (fname, "frames/output%05d.pam", n);
    //FILE * fin = fopen ("output00316.pam", "r");
    FILE * fin = fopen (fname, "r");
    if (! fin)
      {
        perror ("fin");
        return 1;
      }

    //pm_init (argv [0], 0);
    pm_init ("ptr", 0);

    //pnm_readpaminit (stdin, & inpam, PAM_STRUCT_SIZE (tuple_type));
    pnm_readpaminit (fin, & inpam, PAM_STRUCT_SIZE (tuple_type));

    // Read and convert to grayscale

    int h = inpam.height;
    int w = inpam.width;
    //int img [h] [w];
    //int img [h * w];
    fprintf (stderr, "%05d h %d w %d\n", n, h, w);
    static int * img = NULL;
    if (! img)
      img = malloc (h * w * sizeof (int));

    tuplerow = pnm_allocpamrow (& inpam);
    for (row = 0; row < inpam.height; row ++)
      {
        unsigned int column;
        pnm_readpamrow (& inpam, tuplerow);
        for (column = 0; column < inpam.width; column ++)
          {
{ int ii = row * w + column; if (ii < 0 || ii >= h * w) printf ("wha ??\n");}
// Assuming RGB
// 0.21 R + 0.72 G + 0.07 B
            img [row * w + column] =  tuplerow [column] [0] * 0.21 + tuplerow [column] [1] * 0.72 + tuplerow [column] [2] * 0.07;
            if (img [row * w + column] < 0)
              img [row * w + column] = 0;
            if (img [row * w + column] > 255)
              img [row * w + column] = 255;
          }
      }
    pnm_freepamrow (tuplerow);
    fclose (fin);
    //dumpImage (img, h, w, "1_mono", n);

    //int edgedata [h * w];
    //int markup [h * w];
    //int hough [h * w];
    static int * edgedata = NULL;
    if (! edgedata)
      edgedata = malloc (h * w * sizeof (int));
    static int * markup = NULL;
    if (! markup)
      markup = malloc (h * w * sizeof (int));
    static int * hough = NULL;
    if (! hough)
      hough = malloc (h * w * sizeof (int));
    static int * straight = NULL;
    if (! straight)
      straight = malloc (h * w * sizeof (int));

// Edge detection

// http://dasl.mem.drexel.edu/alumni/bGreen/www.pages.drexel.edu/_weg22/edge.html
    int X, Y, I;
    //int * imgp = & img [0];
    int * imgp = img;
    //int edgedata [h * w];

   int GX[3][3];
   int GY[3][3];

   /* 3x3 GX Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
   GX[0][0] = -1; GX[0][1] = 0; GX[0][2] = 1;
   GX[1][0] = -2; GX[1][1] = 0; GX[1][2] = 2;
   GX[2][0] = -1; GX[2][1] = 0; GX[2][2] = 1;

   /* 3x3 GY Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
   GY[0][0] =  1; GY[0][1] =  2; GY[0][2] =  1;
   GY[1][0] =  0; GY[1][1] =  0; GY[1][2] =  0;
   GY[2][0] = -1; GY[2][1] = -2; GY[2][2] = -1;
    for(Y = 0; Y <= (h - 1); Y ++)
      {
        for(X = 0; X <= (w - 1); X ++)
          {
            int sumX = 0;
            int sumY = 0;
            int SUM;

            /* image boundaries */
            if (Y == 0 || Y == h - 1)
              SUM = 0;
            else if (X == 0 || X == w - 1)
              SUM = 0;

            /* Convolution starts here */
            else
              {

                /*-------X GRADIENT APPROXIMATION------*/
                for (int I = -1; I <= 1; I ++)
                  {
                    for(int J = -1; J <= 1; J ++)
                      {
{ int ii = X + I + (Y + J) * w; if (ii < 0 || ii >= h * w) printf ("wha ??\n");}
                        sumX = sumX + (int) ((* (imgp + X + I + 
                             (Y + J) * w)) * GX [I + 1] [J + 1]);
                      }
                  }

                /*-------Y GRADIENT APPROXIMATION-------*/
                for (int I = -1; I <= 1; I ++)
                  {
                    for(int J = -1; J <= 1; J ++)
                      {
                        sumY = sumY + (int) ((* (imgp + X + I + 
                                (Y + J) * w)) * GY [I + 1] [J + 1]);
                      }
                  }

                /*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
                SUM = abs (sumX) + abs (sumY);
              }

            if (SUM > 255)
              SUM = 255;
            if (SUM < 0)
              SUM = 0;

            //edgedata [X + Y * w] = 255 - (unsigned char) (SUM);
            //edgedata [X + Y * w] = (unsigned char) (SUM);
            edgedata [X + Y * w] = SUM == 255 ? 255 : 0;
          }
      }
    dumpImage (edgedata, h, w, "2_edge", n);

// Hough transform

    //int hough [h * w];
    houghtransform (edgedata, hough, w, h);

    //dumpImage (hough, h, w, "3_hough", n);
#if 0
// Scale hough down for debug image
    int shough [h * w];
    int max = hough [0];
    for (int i = 0; i < h * w; i ++)
      {
        shough [i] = hough [i];
        if (hough [i] > max)
          max = hough [i];
      }
    if (max > 255)
      {
        int d = max - 255;
        for (int i = 0; i < h * w; i ++)
          {
            shough [i] = hough [i] * 255 / max;
            if (shough [i] < 0)
              shough [i] = 0;
          }
      }
    //dumpImage (shough, h, w, "4_shough", n);
#endif

    int rhoWidM1 = w - 1;
    double rhoNorm = rhoWidM1 / sqrt ((double) (w * w) + (double) h * h);

#if 0
#define threshold 10

// Find candidates above threshold
    for (int i = 0; i < h * w; i ++)
      {
        int v = hough [i];
        if (v > threshold)
          {
            int x = i % w;
            int y = i / w;
            double rho = x / rhoNorm;
            double theta = y * PI / (double) h;
            fprintf (stderr, "%05d x %d y %d v %d rho %6.1f theta %0.4f\n", n, x, y, v, rho, theta);
          }
      }
#endif
#if 0
    // Find the two strongest lines.
    int s1 = hough [0];
    int si1 = 0;
    int s2 = s1;
    int si2 = 0;
    for (int i = 0; i < h * w; i ++)
      {
        int v = hough [i];
        if (v > s1)
          {
            s1 = v;
            si1 = i;
          }
      }

    for (int i = 0; i < h * w; i ++)
      {
        int v = hough [i];
        if (v > s2 && i != si1)
          {
            s2 = v;
            si2 = i;
          }
      }

    int rhoWidM1 = w - 1;
    double rhoNorm = rhoWidM1 / sqrt ((double) (w * w) + (double) h * h);
    int s1x = si1 % w;
    int s1y = si1 / w;
    double rho1 = s1x / rhoNorm;
    double theta1 = s1y * PI / (double) h;
    int s2x = si2 % w;
    int s2y = si2 / w;
    double rho2 = s2x / rhoNorm;
    double theta2 = s2y * PI / (double) h;
    fprintf (stderr, "%05d s1 %d si1 %d rho %f theta %f\n", n, s1, si1, rho1, theta1);
    fprintf (stderr, "%05d s2 %d si2 %d rho %f theta %f\n", n, s2, si2, rho2, theta2);
#endif

// image is in top-left corner;\n");
//  rho increases along horizontal axis to
//  maximum rho equal to image diagonal length;
//  theta increases downward from 0 radians to PI radians.

//
//   +-----------------------------
//   |\
//   | \
//   |  \    .
//   |   .
//   .
//

// http://stackoverflow.com/questions/7613955/hough-transform-equation
//
// And to transform between the two, use the equation 
//   y = -(cos(theta)/sin(theta))x + r/sin(theta).
// Thus
//   m = -(cos(theta)/sin(theta))
// and
//   b = r/sin(theta). 
// These obviously break down when sin(theta)=0 or theta=0, which is why the
// rotational coordinate system is preferred (there aren't any problems with
// lines with infinite slopes).

#define solve(x, theta, r) ( -(cos(theta)/sin(theta)) * x + r/sin(theta))

#if 0
    // Solve for x = 0
    double s1y_x0 = solve (0, PI/2 - theta1, rho1);
    // Solve for x = w - 1 (right edge)
    double s1y_x1 = solve (w - 1, PI/2 - theta1, rho1);

    // Solve for x = 0
    double s2y_x0 = solve (0, PI/2 - theta2, rho2);
    // Solve for x = w - 1 (right edge)
    double s2y_x1 = solve (w - 1, PI/2 - theta2, rho2);

    fprintf (stderr, "s1y_x0 %f s1y_x1 %f\n", s1y_x0, s1y_x1);
    fprintf (stderr, "s2y_x0 %f s2y_x1 %f\n", s2y_x0, s2y_x1);
    // Debug; dump an copy of the original with the edges drawn in black
#endif

#define idx2x(idx) (idx % w)
#define idx2y(idx) (idx / w)
#define idx2rt(idx, rho, theta) \
  { \
    int rhoWidM1 = w - 1; \
    double rhoNorm = rhoWidM1 / sqrt ((double) (w * w) + (double) h * h); \
 \
    int x = idx2x (idx); \
    int y = idx2y (idx); \
    rho = x / rhoNorm; \
    theta = y * PI / (double) h; \
    if (theta > PI / 2) \
      theta -= PI; \
  }

// width 512
//#define min_ht_v 300
// width 1080
#define min_ht_v 400

// Find the strongest line

    int l1v = hough [0];
    int l1i = 0;
    for (int i = 0; i < h * w; i ++)
      {
        int v = hough [i];
        if (v > l1v)
          {
            l1v = v;
            l1i = i;
          }
      }

    if (l1v < min_ht_v)
      {
        fprintf (stderr, "%05d can't find good edge\n", n);
        return 0;
      }

    double rho1, theta1;
    idx2rt (l1i, rho1, theta1);

    fprintf (stderr, "%05d strongest v %d rho %6.1f theta %0.4f\n", n, l1v, rho1, theta1);

// Find the strongest parallel line at least 100 pixels away

// the x index is the de-normalized rho; ~= to pixels
#define dx_threshold 100
// the y index is the de-normalized theta; defining parallel as dy < 5
#define dy_threshold 5

    int l1x = idx2x (l1i);
    int l1y = idx2y (l1i);

    int l2v = hough [0];
    int l2i = -1;
    for (int i = 0; i < h * w; i ++)
      {
        int x = idx2x (i);
        int y = idx2y (i);
        // not parallle?
        if (abs (y - l1y) >= dy_threshold)
          continue;
        // too close?
        if (abs (x - l1x) < dx_threshold)
          continue;
        int v = hough [i];
        if (v > l2v)
          {
            l2v = v;
            l2i = i;
          }
      }

    fprintf (stderr, "%05d l2i %d l2v %d\n", n, l2i, l2v);
    if (l2i == -1)
      {
        fprintf (stderr, "%05d can't find other edge\n", n);
        return 0;
      }

    if (l2v < min_ht_v)
      {
        fprintf (stderr, "%05d can't find good parallel edge; l2v %d\n", n, l2v);
        return 0;
      }
    double rho2, theta2;
    idx2rt (l2i, rho2, theta2);

    fprintf (stderr, "%05d parallel  v %d rho %6.1f theta %0.4f\n", n, l2v, rho2, theta2);


    // Solve for x = 0
    double l1y_x0 = solve (0, PI/2 + theta1, rho1);
    // Solve for x = w - 1 (right edge)
    double l1y_x1 = solve (w - 1, PI/2 + theta1, rho1);

    // Solve for x = 0
    double l2y_x0 = solve (0, PI/2 + theta2, rho2);
    // Solve for x = w - 1 (right edge)
    double l2y_x1 = solve (w - 1, PI/2 + theta2, rho2);

    fprintf (stderr, "%05d l1y_x0 %d l1y_x1 %d l2y_x0 %d l2y_x1 %d\n", n, (int) l1y_x0, (int) l1y_x1, (int) l2y_x0, (int) l2y_x1);


#if 1
// Debug: write marked-up originals

    //int markup [h * w];

    for (int i = 0; i < h * w; i ++)
      markup [i] = img [i];

    drawline (markup, h, w, 0, l1y_x0, w - 1, l1y_x1, 0);
    drawline (markup, h, w, 0, l2y_x0, w - 1, l2y_x1, 255);
    dumpImage (markup, h, w, "5_markup", n);
#endif

// http://www.cprogramto.com/c-program-for-dda-algorithm/


// Straighten the tape by moving columns up or down and center.

    // // midpoint of the left end of the tape
    // int my = (l1y_x0 + l2y_x0) / 2;
    // // number of pixels beyond the center line
    // int cy = my - (h / 2);
    // difference in height between the left and right ends
    // Which edge is lower
    int low = l1y_x0;
    if (l2y_x0 > low)
      low = l2y_x0;
    int dlow = low - (h / 2) - 150;
    // And down 150 pixels
    dlow += 150;
    int dy = l1y_x1 - l1y_x0;
//printf ("l1y_x0 %d l2y_x0 %d low %d dlow %d %d dy %d\n", (int) l1y_x0, (int) l2y_x0, low, low - (h / 2), dlow, dy);
    //int straight [h * w];
    for (int x = 0; x < w; x ++)
      {
        int xdy = x * dy / w;
        //fprintf (stderr, "%d %d\n", x, xdy);
        for (int y = 0; y < h; y ++)
          {
            int from = y + xdy + dlow;
            if (from < 0)
              from = 0;
            if (from > h - 1)
              from = h - 1;
            straight [y * w + x] = img [from * w + x];
          }
      }
    fprintf (stderr, "%05d nwrite %05d\n", n, nwrite);
    dumpImage (straight, h, w, "6_straight", nwrite);
    nwrite ++;
    return 0;
  }

int main (int argc, char * argv [])
  {
#if 0
    int n;
    if (argc != 2)
      {
        fprintf (stderr, "Usage 'ptr n' where n is the image number\n");
        exit (1);
      }
    n = atoi (argv [1]);
#endif
    for (int i = 1; ; i ++)
    //for (int i = 180; i <= 190; i ++)
      if (image (i))
        break;
   }

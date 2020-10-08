/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/**                                                                   **/
/** This suite is a modification of the UCR Suite and, consequently,  **/
/** follows the same copyright (Transcribed below, without any change)**/
/**                                                                   **/
/** This modified version is the responsability of Matthieu Herrmann  **/
/** and Geoff I. Webb.                                                **/
/**                                                                   **/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>

#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code
#define POSITIVE_INFINITY 1e20       //Pseudo Infitinte number for this code

using namespace std;

/// Data structure for sorting the query
typedef struct Index
{   double value;
  int    index;
} Index;


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
  Index* y = (Index*)b;
  return abs(y->value) - abs(x->value);   // high to low
}

/** Square distances between two doubles */
inline double square_dist(double v1, double v2) {
  double d = v1 - v2;
  return d * d;
}

/** Unsigned arithmetic:
 * Given an 'index' and a 'window', get the start index corresponding to std::max(0, index-window) */
inline std::size_t cap_start_index_to_window(std::size_t index, std::size_t window){
  if(index>window){ return index-window; } else { return 0; }
}

/** Unsigned arithmetic:
 * Given an 'index', a 'window' and an 'end', get the stop index corresponding to std::min(end, index+window+1).
 * The expression index+window+1 is illegal for any index>0 as window could be MAX-1
 * */
inline std::size_t cap_stop_index_to_window_or_end(std::size_t index, std::size_t window, std::size_t end){
  // end-window is valid when window<end
  if(window<end && index+1<end-window){ return index + window + 1; } else { return end; }
}


double dtw(
    const double *lines,
    const double *cols,
    int l,
    int w,
    double bsf
    )
{
  // 1) --- Create the upper bound from bsf and a margin
  const double UB = bsf;

  // 2) --- Alias in line/column concept: we only allocate for the columns, using the smallest possible dimension.
  const size_t nbcols = l;
  const size_t nblines = l;

  // 3) --- Cap the windows.
  if (w > nblines) { w = nblines; }

  // 4) --- Buffers allocations
  // Add an extra column for the "matrix border" condition, init to +INF.
  // Using an unique contiguous array. Base indices are:
  // 'c' for the current row,
  // 'p' for the previous one
  std::vector<double> buffers_v((1+nbcols) * 2, POSITIVE_INFINITY);
  double *buffers = buffers_v.data();
  size_t c{1}, p{nbcols+2};                 // Account for the extra column (in front)

  // 5) --- Computation of DTW
  buffers[c-1] = 0;
  size_t next_start{0};
  size_t pruning_point{0};

  for(size_t i=0; i<nblines; ++i) {
    // --- --- --- --- Swap and variables init
    std::swap(c, p);
    const double li = lines[i];
    const std::size_t jStop = cap_stop_index_to_window_or_end(i, w, nbcols);
    const std::size_t jStart = std::max(cap_start_index_to_window(i, w), next_start);
    std::size_t next_pruning_point = jStart; // Next pruning point init at the start of the line
    std::size_t j = jStart;
    next_start = jStart;
    // --- --- --- --- Init the first column
    buffers[c+j-1] = POSITIVE_INFINITY;
    double cost = POSITIVE_INFINITY;
    // --- --- --- --- Compute DTW up to the pruning point while advancing next_start: diag and top
    for(; j==next_start && j < pruning_point; ++j) {
      const auto d = square_dist(li, cols[j]);
      cost = std::min(buffers[p + j - 1], buffers[p + j]) + d;
      buffers[c + j] = cost;
      if(cost<=UB){ next_pruning_point = j + 1;} else { ++next_start; }
    }
    // --- --- --- --- Compute DTW up to the pruning point without advancing next_start: prev, diag, top
    for(; j < pruning_point; ++j) {
      const auto d = square_dist(li, cols[j]);
      cost = std::min(cost, std::min(buffers[p + j - 1], buffers[p + j])) + d;
      buffers[c + j] = cost;
      if(cost<=UB){ next_pruning_point = j + 1;}
    }
    // --- --- --- --- Compute DTW at "pruning_point": 2 cases
    if(j<jStop){
      const auto d = square_dist(li, cols[j]);
      if(j==next_start){ // Advancing next start: only diag. Done if v>UB.
        cost = buffers[p + j - 1] + d;
        buffers[c + j] = cost;
        if(cost<=UB){ next_pruning_point = j + 1;} else {return POSITIVE_INFINITY; }
      } else { // Not advancing next start: at least a path possible in previous cells.
        cost = std::min(cost, buffers[p + j - 1]) + d;
        buffers[c + j] = cost;
        if(cost<=UB){ next_pruning_point = j + 1;}
      }
      ++j;
    } else if(j==next_start) { return POSITIVE_INFINITY; }
    // --- --- --- --- Compute DTW after "pruning_point": prev. Go on while we advance the next pruning point.
    for(;j==next_pruning_point && j<jStop;++j){
      const auto d = square_dist(li, cols[j]);
      cost = cost + d;
      buffers[c + j] = cost;
      if(cost<=UB){ ++next_pruning_point; }
    }

    pruning_point=next_pruning_point;
  }// End for i loop


  // 6) --- If the pruning_point did not reach the number of columns, we pruned something
  if(pruning_point != nbcols){ return POSITIVE_INFINITY; } else {
    return buffers[c+nbcols-1];
  }
}




/// Print function for debugging
void printArray(double *x, int len)
{   for(int i=0; i<len; i++)
  printf(" %6.2lf",x[i]);
  printf("\n");
}

/// If expected error happens, teminated the program.
void error(int id)
{
  if(id==1)
    printf("ERROR : Memory can't be allocated!!!\n\n");
  else if ( id == 2 )
    printf("ERROR : File not Found!!!\n\n");
  else if ( id == 3 )
    printf("ERROR : Can't create Output File!!!\n\n");
  else if ( id == 4 )
  {
    printf("ERROR : Invalid Number of Arguments!!!\n");
    printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
    printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
  }
  exit(1);
}

/// Main Function
int main(  int argc , char *argv[] )
{
  FILE *fp;            /// data file pointer
  FILE *qp;            /// query file pointer
  double bsf;          /// best-so-far
  double *t, *q;       /// data array and query array
  int *order;          ///new order of the query
  double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;


  double d;
  long long i , j;
  double ex , ex2 , mean, std;
  int m=-1, r=-1;
  long long loc = 0;
  double t1,t2;
  int kim = 0,keogh = 0, keogh2 = 0;
  double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
  double *buffer, *u_buff, *l_buff;
  Index *Q_tmp;

  /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
  int EPOCH = 100000;

  /// If not enough input, display an error.
  if (argc<=3)
    error(4);

  /// read size of the query
  if (argc>3)
    m = atol(argv[3]);

  /// read warping windows
  if (argc>4)
  {   double R = atof(argv[4]);
    if (R<=1)
      r = floor(R*m);
    else
      r = floor(R);
  }

  fp = fopen(argv[1],"r");
  if( fp == NULL )
    error(2);

  qp = fopen(argv[2],"r");
  if( qp == NULL )
    error(2);

  /// start the clock
  t1 = clock();


  /// malloc everything here
  q = (double *)malloc(sizeof(double)*m);
  if( q == NULL )
    error(1);
  qo = (double *)malloc(sizeof(double)*m);
  if( qo == NULL )
    error(1);
  uo = (double *)malloc(sizeof(double)*m);
  if( uo == NULL )
    error(1);
  lo = (double *)malloc(sizeof(double)*m);
  if( lo == NULL )
    error(1);

  order = (int *)malloc(sizeof(int)*m);
  if( order == NULL )
    error(1);

  Q_tmp = (Index *)malloc(sizeof(Index)*m);
  if( Q_tmp == NULL )
    error(1);

  u = (double *)malloc(sizeof(double)*m);
  if( u == NULL )
    error(1);

  l = (double *)malloc(sizeof(double)*m);
  if( l == NULL )
    error(1);

  cb = (double *)malloc(sizeof(double)*m);
  if( cb == NULL )
    error(1);

  cb1 = (double *)malloc(sizeof(double)*m);
  if( cb1 == NULL )
    error(1);

  cb2 = (double *)malloc(sizeof(double)*m);
  if( cb2 == NULL )
    error(1);

  u_d = (double *)malloc(sizeof(double)*m);
  if( u == NULL )
    error(1);

  l_d = (double *)malloc(sizeof(double)*m);
  if( l == NULL )
    error(1);

  t = (double *)malloc(sizeof(double)*m*2);
  if( t == NULL )
    error(1);

  tz = (double *)malloc(sizeof(double)*m);
  if( tz == NULL )
    error(1);

  buffer = (double *)malloc(sizeof(double)*EPOCH);
  if( buffer == NULL )
    error(1);

  u_buff = (double *)malloc(sizeof(double)*EPOCH);
  if( u_buff == NULL )
    error(1);

  l_buff = (double *)malloc(sizeof(double)*EPOCH);
  if( l_buff == NULL )
    error(1);


  /// Read query file
  bsf = INF;
  i = 0;
  j = 0;
  ex = ex2 = 0;

  while(fscanf(qp,"%lf",&d) != EOF && i < m)
  {
    ex += d;
    ex2 += d*d;
    q[i] = d;
    i++;
  }
  fclose(qp);

  /// Do z-normalize the query, keep in same array, q
  mean = ex/m;
  std = ex2/m;
  std = sqrt(std-mean*mean);
  for( i = 0 ; i < m ; i++ )
    q[i] = (q[i] - mean)/std;

  /// Sort the query one time by abs(z-norm(q[i]))
  for( i = 0; i<m; i++)
  {
    Q_tmp[i].value = q[i];
    Q_tmp[i].index = i;
  }
  qsort(Q_tmp, m, sizeof(Index),comp);

  /// also create another arrays for keeping sorted envelop
  for( i=0; i<m; i++)
  {   int o = Q_tmp[i].index;
    order[i] = o;
    qo[i] = q[o];
    uo[i] = u[o];
    lo[i] = l[o];
  }
  free(Q_tmp);

  /// Initial the cummulative lower bound
  for( i=0; i<m; i++)
  {   cb[i]=0;
    cb1[i]=0;
    cb2[i]=0;
  }

  i = 0;          /// current index of the data in current chunk of size EPOCH
  j = 0;          /// the starting index of the data in the circular array, t
  ex = ex2 = 0;
  bool done = false;
  int it=0, ep=0, k=0;
  long long I;    /// the starting index of the data in current chunk of size EPOCH

  while(!done)
  {
    /// Read first m-1 points
    ep=0;
    if (it==0)
    {   for(k=0; k<m-1; k++)
      if (fscanf(fp,"%lf",&d) != EOF)
        buffer[k] = d;
    }
    else
    {   for(k=0; k<m-1; k++)
      buffer[k] = buffer[EPOCH-m+1+k];
    }

    /// Read buffer of size EPOCH or when all data has been read.
    ep=m-1;
    while(ep<EPOCH)
    {   if (fscanf(fp,"%lf",&d) == EOF)
      break;
      buffer[ep] = d;
      ep++;
    }

    /// Data are read in chunk of size EPOCH.
    /// When there is nothing to read, the loop is end.
    if (ep<=m-1) {   done = true; }
    else {

      /// Just for printing a dot for approximate a million point. Not much accurate.
      if (it%(1000000/(EPOCH-m+1))==0)
        fprintf(stderr,".");

      /// Do main task here..
      ex=0; ex2=0;
      for(i=0; i<ep; i++)
      {
        /// A bunch of data has been read and pick one of them at a time to use
        d = buffer[i];

        /// Calcualte sum and sum square
        ex += d;
        ex2 += d*d;

        /// t is a circular array for keeping current data
        t[i%m] = d;

        /// Double the size for avoiding using modulo "%" operator
        t[(i%m)+m] = d;

        /// Start the task when there are more than m-1 points in the current chunk
        if( i >= m-1 )
        {
          mean = ex/m;
          std = ex2/m;
          std = sqrt(std-mean*mean);

          /// compute the start location of the data in the current circular array, t
          j = (i+1)%m;
          /// the start location of the data in the current chunk
          I = i-(m-1);


          /// Take another linear time to compute z_normalization of t.
          /// Note that for better optimization, this can merge to the previous function.
          for(k=0;k<m;k++) {   tz[k] = (t[(k+j)] - mean)/std; }

          /// Compute DTW and early abandoning if possible
          dist = dtw(tz, q, m, r, bsf);

          if( dist < bsf )
          {   /// Update bsf
            /// loc is the real starting location of the nearest neighbor in the file
            bsf = dist;
            loc = (it)*(EPOCH-m+1) + i-m+1;
          }

          /// Reduce obsolute points from sum and sum square
          ex -= t[j];
          ex2 -= t[j]*t[j];
        }
      }

      /// If the size of last chunk is less then EPOCH, then no more data and terminate.
      if (ep<EPOCH)
        done=true;
      else
        it++;
    }
  }

  i = (it)*(EPOCH-m+1) + ep;
  fclose(fp);

  free(q);
  free(u);
  free(l);
  free(uo);
  free(lo);
  free(qo);
  free(cb);
  free(cb1);
  free(cb2);
  free(tz);
  free(t);
  free(l_d);
  free(u_d);
  free(l_buff);
  free(u_buff);

  t2 = clock();
  printf("\n");

  /// Note that loc and i are long long.
  cout << "Location : " << loc << endl;
  cout << "Distance : " << sqrt(bsf) << endl;
  cout << "Data Scanned : " << i << endl;
  cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

  /// printf is just easier for formating ;)
  printf("\n");
  printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
  printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
  printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
  printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
  return 0;
}

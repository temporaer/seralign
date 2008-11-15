// Original readdata routines are in orig_readdata.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR)

#define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

int quicksort5(int*, int*, int*, int*, double*, int, int);
int partition5(int*, int*, int*, int*, double*, int, int);
void skip_to_end_of_line(FILE*);
int get_line(FILE*, char*, int);
int max_line_length(FILE*);

extern void dsyev_();

int readdata_sdpa(char* datafilename, int* passed_m, int* passed_numblk, int** passed_blksz,
                  char** passed_blktype, double** passed_b, double** passed_CAent,
                  int** passed_CArow, int** passed_CAcol, int** passed_CAinfo_entptr,
                  int** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                  char** passed_CAinfo_storage)
{
  // declare variables 
  int i, ct, buflen, ret, num, pos1, pos2, blk, *nnz, *temp_num, *temp_blk, specindex;
  char *buf, c, *ptr1, *ptr2;
  double entry;
  FILE *datafile;

  int    m, numblk, *blksz;
  char   *blktype;
  double *b, *CAent;
  int    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;


  ////////////////////////////////
  // Step 0
  ////////////////////////////////

  // Open file, and determine the length of the longest line. 
  // Note: this can fail if last line doesn't have a new line character; how to fix? 

  datafile = fopen(datafilename, "r");
  if(datafile == NULL) {
    printf("Can't get file %s\n", datafilename);
    exit(0);
  }
  buflen = max_line_length(datafile) + 10;
  fclose(datafile);

  // Allocate buffer space, which will hold each line that is read in

  MYCALLOC(buf, char, buflen);


  ////////////////////////////////
  // Step 1
  ////////////////////////////////


  // Get m, numblk, *blksz, *blktype, *b;
  // size of partitions of (CAent, CArow, CAcol) into (data matrix)-(block) pairs

  // Open file

  datafile = fopen(datafilename, "r");

  // Read through the comment lines. 
 
  c = getc(datafile);
  while(c == '"' || c == '*') {
    skip_to_end_of_line(datafile);
    c = getc(datafile);
  }
  ungetc(c, datafile);

  // Get m

  ret = get_line(datafile, buf, buflen);
  sscanf(buf, "%d", &m);

  // Get numblk

  ret = get_line(datafile, buf, buflen);
  sscanf(buf, "%d", &numblk);

  // Prepare to get blksz and blktype
  MYCALLOC(blksz, int, numblk);
  MYCALLOC(blktype, char, numblk);

  // Get blksz
  ret = get_line(datafile, buf, buflen);
  ptr1 = buf;
	for(i = 0; i < numblk; i++) {
		blksz[i] = strtol(ptr1, &ptr2, 10);
	  ptr1 = ptr2;
  }

  // Get blktype (easy from blksz)
  // Note: In SDPA format, X nonneg variables constitutes
  // one (diagonal) block of size X, which is indicated by
  // a negative block size.

  for(i = 0; i < numblk; i++) {
    if(blksz[i] < 0) {
      blktype[i] = 'd';
      blksz[i] = -blksz[i];
    }
    else blktype[i] = 's';
  }

  // Prepare to get b

  MYCALLOC(b, double, m);

  // Get b 

  ret = get_line(datafile, buf, buflen);
  ptr1 = buf;
  for(i = 0; i < m; i++) {
	  b[i] = strtod(ptr1, &ptr2);
	  ptr1 = ptr2;
	}

  // Prepare for next step
  
  MYCALLOC(nnz, int, (m+1)*numblk);

  // Determine how many nnz entries are in each (data matrix)-(block) pair

  ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  do {
    if(fabs(entry) > 1.0e-15) nnz[ DATABLOCKIND(num,blk,numblk) ]++;
    ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  } while(ret == 5);

  // Finished with this pass through file

  fclose(datafile);


  ////////////////////////////////
  // Step 2
  ////////////////////////////////

  // Use information gathered from Step 1 to allocate space
  // for remaining structures and to set CAinfo...

  // Allocate all CAinfo...

  MYCALLOC(CAinfo_entptr,    int,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, int,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type,      char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage,   char, (m+1)*numblk);

  // Setup all CAinfo...

  ct = 0;

  for(i = 0; i <= m; i++)
    for(blk = 1; blk <= numblk; blk++) {

      specindex = DATABLOCKIND(i,blk,numblk);
      CAinfo_entptr[specindex] = ct;
      CAinfo_rowcolptr[specindex] = ct;
      ct += nnz[specindex];

      CAinfo_storage[specindex] = 's'; // SDPA format implies same for all
      if(blktype[blk-1] == 'd') CAinfo_type[specindex] = 'd';
      else CAinfo_type[specindex] = 's';

    }

  CAinfo_entptr[(m+1)*numblk] = ct;
  CAinfo_rowcolptr[(m+1)*numblk] = ct;

  // Allocate remaining space plus some temp space for use in Step 3

  MYCALLOC(CAent, double, ct);
  MYCALLOC(CArow, int, ct);
  MYCALLOC(CAcol, int, ct);
  MYCALLOC(temp_num, int, ct);
  MYCALLOC(temp_blk, int, ct);


  ////////////////////////////////
  // Step 3
  ////////////////////////////////

  // Repeat pass through file, this time only getting
  // specific entry information

  // Open file 
  datafile = fopen(datafilename, "r");

  // Lines: comments
  c = getc(datafile);
  while(c == '"' || c == '*') {
    skip_to_end_of_line(datafile);
    c = getc(datafile);
  }
  ungetc(c, datafile);

  // Line: number of constraints
  ret = get_line(datafile, buf, buflen);

  // Line: number of blocks.
  ret = get_line(datafile, buf, buflen);

  // Line: block sizes 
  ret = get_line(datafile, buf, buflen);

  // Line: b 
  ret = get_line(datafile, buf, buflen);

  // Now, loop through the entries, extracting info. 

  ct = 0;
  ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  do {
    if(fabs(entry) > 1.0e-15) {
      if(pos1 == 0 || pos2 == 0) {
        printf("Error (readdata_sdpa): Encountered '0' row or column index.\n");
        exit(0);
      }
      temp_num[ct] = num;
      temp_blk[ct] = blk;
      if(pos1 > pos2) { CArow[ct] = pos1; CAcol[ct] = pos2; }
      else { CArow[ct] = pos2; CAcol[ct] = pos1; }
      if(num == 0) CAent[ct] = -entry;
      else CAent[ct] = entry;
      ct++;
    }
    ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  } while(ret == 5);

  // Close data file
  fclose(datafile);

  ////////////////////////////////
  // Step 4
  ////////////////////////////////

  // Sort CAent, CArow, CAcol according to (data matrix) -> (block) order

  // Not totally sure this is correct but will go with it for now
  quicksort5(temp_num, temp_blk, CArow, CAcol, CAent, 0, CAinfo_entptr[(m+1)*numblk]-1);

  /////////////////////////////////
  // End Step 3
  /////////////////////////////////

  /////////////////////////////////
  // End Step 2
  /////////////////////////////////

  MYFREE(temp_num);
  MYFREE(temp_blk);

  /////////////////////////////////
  // End Step 1
  /////////////////////////////////

  MYFREE(nnz);

  /////////////////////////////////
  // End Step 0
  /////////////////////////////////

  MYFREE(buf);

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;

}

int quicksort5(int* A1, int* A2, int* A3, int* A4, double* A5, int p, int r)
{
   int q;

   if(p < r) 
   {
      q = partition5(A1, A2, A3, A4, A5, p, r);
      quicksort5(A1, A2, A3, A4, A5, p, q);
      quicksort5(A1, A2, A3, A4, A5, q+1, r);
   }

   return 1;
}

int partition5(int* A1, int* A2, int* A3, int* A4, double* A5, int p, int r)
{
   int i, j;
   int sv1, sv2;
   int t1, t2, t3, t4;
   double t5;

   sv1 = A1[p]; sv2 = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] > sv1 || (A1[j] == sv1 && A2[j] > sv2) );
      do i++;
      while(A1[i] < sv1 || (A1[i] == sv1 && A2[i] < sv2) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j]; t3 = A3[j]; t4 = A4[j]; t5 = A5[j];
         A1[j] = A1[i]; A2[j] = A2[i]; A3[j] = A3[i]; A4[j] = A4[i]; A5[j] = A5[i];
         A1[i] = t1; A2[i] = t2; A3[i] = t3; A4[i] = t4; A5[i] = t5;
      }
      else return j;
   }

   return 0;
}


// This routine skips to the end of the current line of input from the file datafile. 

void skip_to_end_of_line(FILE *datafile)
{
  char c;
  
  c = getc(datafile);
  while (c != '\n') c = getc(datafile);
}

// This routine reads a line of input into a buffer, and translates all
//   occurences of "," "{" "}" "(" ")" to spaces. 

int get_line(FILE *datafile, char *buffer, int bufsiz)
{
  int i, k=0;
  char c;
  
  c = getc(datafile);
  while (c != '\n') {
    buffer[k] = c;
    k++;
    c = getc(datafile);
    if(c == EOF) return(2);
    if(k >= bufsiz) {
      printf("Line too long in input file!  Adjust BUFFERSIZ in readprob.c\n");
      return(1);
    }
  }
  buffer[k] = '\n';
  buffer[k+1] = 0;

  for(i = 0; i <= k; i++) {
    if(buffer[i] == ',' || buffer[i] == '{' ||
       buffer[i] == '}' || buffer[i] == '(' ||
       buffer[i]==')') buffer[i]=' ';
  }

  return(0);
}

int max_line_length(FILE *datafile)
{
  int maxlen=0, k=0, c;
  
  c = getc(datafile);
  while(c != EOF) {
    k = 0;
    while(c != '\n') {
      c = getc(datafile);
      k++;
    }
    if(k > maxlen) maxlen=k;
    c = getc(datafile);
  }
  
  return(maxlen);
}



int readdata_sdplr(char* datafilename, int* passed_m, int* passed_numblk, int** passed_blksz,
                   char** passed_blktype, double** passed_b, double** passed_CAent,
                   int** passed_CArow, int** passed_CAcol, int** passed_CAinfo_entptr,
                   int** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                   char** passed_CAinfo_storage)
{
  // This code assumes a very specific structure of input file,
  // and no error checking is done at all.

  int i, j, k, *nnz, num, pos1, pos2, sz;
  int ct, blk, specindex;
  double entry;
  char type;
  FILE *datafile;

  int    m, numblk, *blksz;
  char   *blktype;
  double *b, *CAent;
  int    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;

  ////////////////////////////////
  // Step 1
  ////////////////////////////////

  // Open file
  datafile = fopen(datafilename, "r");
  if(datafile == NULL) {
    printf("Error (readdata_sdplr): Can't get file %s.\n", datafilename);
    exit(0);
  }

  // Read m
  fscanf(datafile, "%d", &m);

  // Read numblk
  fscanf(datafile, "%d", &numblk);

  //if(numblk != 1) { printf("SDPLR data file has more than one block!\n"); exit(0); }

  // Setup blksz and blktype
  MYCALLOC(blksz, int, numblk);
  MYCALLOC(blktype, char, numblk);

  // Read block size and set blktype
  for(blk = 0; blk < numblk; blk++) {
    fscanf(datafile, "%d", &(blksz[blk]));
    if(blksz[blk] > 0) blktype[blk] = 's';
    else if(blksz[blk] < 0) {
      blksz[blk] *= -1;
      blktype[blk] = 'd';
    }
    else { printf("Problem reading data. Block size 0!\n"); exit(0); }
  }    
  

  // Allocate space for b
  MYCALLOC(b, double, m);

  // Read b
  for(i = 0; i < m; i++)
    fscanf(datafile, "%lf", &(b[i]));

  // Read eval-adjust-val
  // Right now, we will just discard this
  fscanf(datafile, "%lf", &entry);

  // Prepare for next step  
  MYCALLOC(nnz, int, (m+1)*numblk);

  // Determine how many nnz entries are in each (data matrix)-(block) pair
  // This needs work!
  for(i = 0; i <= m; i++) {
    for(k = 0; k < numblk; k++) {
      fscanf(datafile, "%d %d %c %d", &num, &blk, &type, &sz);
      if(type == 's') {
        nnz[ DATABLOCKIND(num,blk,numblk) ] = 0;
        for(j = 1; j <= sz; j++) {
          fscanf(datafile, "%d %d %lf", &pos1, &pos2, &entry);
          if(fabs(entry) > 1.0e-15) nnz[ DATABLOCKIND(num,blk,numblk) ]++;
        }
      }
      else if(type == 'l') {
        nnz[ DATABLOCKIND(num,blk,numblk) ] = -sz*(blksz[blk-1] + 1); // neg num will allow to identify low-rank data matrices later
        for(j = 1; j <= sz*(blksz[blk-1] + 1); j++)
          fscanf(datafile, "%lf", &entry);
      }
    }
  }

  // Finished with this pass through file
  fclose(datafile);

  ////////////////////////////////
  // Step 2
  ////////////////////////////////

  // Allocate all CAinfo...
  MYCALLOC(CAinfo_entptr,    int,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, int,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type,      char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage,   char, (m+1)*numblk);

  // Setup all CAinfo...
  ct = 0;

  for(i = 0; i <= m; i++)
    for(blk = 1; blk <= numblk; blk++) {

      specindex = DATABLOCKIND(i,blk,numblk);
      CAinfo_entptr[specindex] = ct;
      CAinfo_rowcolptr[specindex] = ct;
      if(nnz[specindex] > 0) ct += nnz[specindex]; // sparse matrix
      else ct += -nnz[specindex]; // low-rank matrix

      if(nnz[specindex] > 0) CAinfo_storage[specindex] = 's'; // sparse matrix (sparse)
      else CAinfo_storage[specindex] = 'd'; // low-rank matrix (dense)
      if(blktype[blk-1] == 'd') CAinfo_type[specindex] = 'd'; // this should never occur in this input format
      else if(nnz[specindex] > 0) CAinfo_type[specindex] = 's';
      else CAinfo_type[specindex] = 'l';

    }

  CAinfo_entptr[(m+1)*numblk] = ct;
  CAinfo_rowcolptr[(m+1)*numblk] = ct;

  // Allocate remaining space plus some temp space for use in Step 3

  MYCALLOC(CAent, double, ct);
  MYCALLOC(CArow, int, ct);
  MYCALLOC(CAcol, int, ct);


  ////////////////////////////////
  // Step 3
  ////////////////////////////////

  // Open file
  datafile = fopen(datafilename, "r");

  // Redo first lines in file
  fscanf(datafile, "%d", &m);
  fscanf(datafile, "%d", &numblk);
  for(blk = 0; blk < numblk; blk++) {
    fscanf(datafile, "%d", &(blksz[blk])); 
    if(blksz[blk] < 0) blksz[blk] *= -1;
  }
  for(i = 0; i < m; i++) fscanf(datafile, "%lf", &(b[i]));
  fscanf(datafile, "%lf", &entry);

  // Read entries in file into data structures
  // This needs work.
  for(i = 0; i <= m; i++) {
    for(k = 0; k < numblk; k++) {
      fscanf(datafile, "%d %d %c %d", &num, &blk, &type, &sz);
      ct = CAinfo_entptr[DATABLOCKIND(i,blk,numblk)];
      if(type == 's') {
        for(j = 1; j <= sz; j++) {
          fscanf(datafile, "%d %d %lf", &pos1, &pos2, &entry);
          if(fabs(entry) > 1.0e-15) {
            if(pos1 == 0 || pos2 == 0) { printf("Error (readdata_sdplr): Encountered '0' row or column index.\n"); exit(0); }
            if(pos1 > pos2) { CArow[ct] = pos1; CAcol[ct] = pos2; }
            else { CArow[ct] = pos2; CAcol[ct] = pos1; }
            CAent[ct] = entry;
            ct++;
          }
        }
      }
      else if(type == 'l') {
        for(j = 1; j <= sz*(blksz[blk-1] + 1); j++) {
          fscanf(datafile, "%lf", &entry);
          if(j <= sz) { pos1 = j; pos2 = j; }
          else {
            pos1 = (j - sz)%blksz[blk-1];
            pos2 = (j - sz)/blksz[blk-1] + 1;
            if(pos1 == 0) { pos1 = blksz[blk-1]; pos2--; }
          }
          CArow[ct] = pos1; CAcol[ct] = pos2;
          CAent[ct] = entry;
          ct++;
        }
      }
    }
  }

  // Finished with this pass through file
  fclose(datafile);

  /////////////////////////////////
  // End Step 3
  /////////////////////////////////

  /////////////////////////////////
  // End Step 2
  /////////////////////////////////

  /////////////////////////////////
  // End Step 1
  /////////////////////////////////

  MYFREE(nnz);

  /////////////////////////////////
  // End Step 0
  /////////////////////////////////

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;
}


int readdata_raw(char* datafilename, int* passed_m, int* passed_numblk, int** passed_blksz,
                char** passed_blktype, double** passed_b, double** passed_CAent,
                int** passed_CArow, int** passed_CAcol, int** passed_CAinfo_entptr,
                int** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                char** passed_CAinfo_storage)
{
  int h, i, k;
  FILE *fid;

  int    m, numblk, *blksz;
  char   *blktype;
  double *b, *CAent;
  int    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;

  fid = fopen(datafilename, "r");
  if(fid == NULL) {
    printf("Can't get file %s\n", datafilename);
    exit(0);
  }

  fscanf(fid, "%d\n", &m);

  fscanf(fid, "%d\n", &numblk);

  MYCALLOC(blksz, int, numblk);
  MYCALLOC(blktype, char, numblk);
  MYCALLOC(b, double, m);

  for(k = 1; k <= numblk; k++)
    fscanf(fid, "%d %c\n", &(blksz[k-1]), &(blktype[k-1]));

  for(i = 1; i <= m; i++)
    fscanf(fid, "%lf\n", &(b[i-1]));

  MYCALLOC(CAinfo_entptr, int, (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, int, (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type, char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage, char, (m+1)*numblk);

  for(h = 0; h < (m+1)*numblk; h++)
    fscanf(fid, "%d %d %c %c\n", &(CAinfo_entptr[h]), &(CAinfo_rowcolptr[h]), &(CAinfo_type[h]), &(CAinfo_storage[h]));
  
  fscanf(fid, "%d %d\n", &(CAinfo_rowcolptr[(m+1)*numblk]), &(CAinfo_entptr[(m+1)*numblk]));

  MYCALLOC(CArow, int, CAinfo_rowcolptr[(m+1)*numblk]);
  MYCALLOC(CAcol, int, CAinfo_rowcolptr[(m+1)*numblk]);
  MYCALLOC(CAent, double, CAinfo_entptr[(m+1)*numblk]);

  for(h = 0; h < CAinfo_rowcolptr[(m+1)*numblk]; h++)
    fscanf(fid, "%d %d\n", &(CArow[h]), &(CAcol[h]));

  for(h = 0; h < CAinfo_entptr[(m+1)*numblk]; h++)
    fscanf(fid, "%lf\n", &(CAent[h]));

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;

}


int writedata_raw(char* datafilename, int m, int numblk, int* blksz,
                  char* blktype, double* b, double* CAent,
                  int* CArow, int* CAcol, int* CAinfo_entptr,
                  int* CAinfo_rowcolptr, char* CAinfo_type,
                  char* CAinfo_storage)
{
  int h, i, k;
  FILE *fid;

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_raw): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++)
    fprintf(fid, "%d %c\n", blksz[k-1], blktype[k-1]);

  for(i = 1; i <= m; i++)
    fprintf(fid, "%.16e\n", b[i-1]);

  for(h = 0; h < (m+1)*numblk; h++)
    fprintf(fid, "%d %d %c %c\n", CAinfo_entptr[h], CAinfo_rowcolptr[h], CAinfo_type[h], CAinfo_storage[h]);
  
  fprintf(fid, "%d %d\n", CAinfo_rowcolptr[(m+1)*numblk], CAinfo_entptr[(m+1)*numblk]);

  for(h = 0; h < CAinfo_rowcolptr[(m+1)*numblk]; h++)
    fprintf(fid, "%d %d\n", CArow[h], CAcol[h]);

  for(h = 0; h < CAinfo_entptr[(m+1)*numblk]; h++)
    fprintf(fid, "%.16e\n", CAent[h]);

  fclose(fid);

  return 0;
}

int writedata_sdplr(char* datafilename, int m, int numblk, int* blksz,
                    char* blktype, double* b, double* CAent,
                    int* CArow, int* CAcol, int* CAinfo_entptr,
                    int* CAinfo_rowcolptr, char* CAinfo_type,
                    char* CAinfo_storage)
{
  int h, i, j, k;
  int r, c, nnz, sz, info, lwork, rank=0, maxsz;
  char jobz, uplo;
  double *MM, *w, *work, maxeval=0.0, eps=1.0e-13;
  FILE *fid;

//   if(numblk != 1) {
//     printf("Error (writedata_sdplr): Format currently only supports one block.\n");
//     return 0;
//   }

  // Anything having to do with eigen assumes only one block
  
  maxsz = -1;
  for(k = 0; k < numblk; k++)
    if(blksz[k] > maxsz) maxsz = blksz[k];

  jobz = 'V';
  uplo = 'L';
  lwork = 3*maxsz-1;
  MYCALLOC(MM, double, maxsz*maxsz);
  MYCALLOC(w, double, maxsz);
  MYCALLOC(work, double, lwork);

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_sdplr): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++) {
    if(blktype[k-1] == 's') fprintf(fid, "%d\n", blksz[k-1]);
    else if(blktype[k-1] == 'd') fprintf(fid, "%d\n", -blksz[k-1]);

  }

  for(i = 1; i <= m; i++)
    fprintf(fid, "%.16e  ", b[i-1]);
  fprintf(fid, "\n");

  fprintf(fid, "-1.0\n");

  for(h = 0; h <= m; h++)
  {
    for(k = 1; k <= numblk; k++)
    {
      if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 's')
      {
        // Experimental: check for low rank if density is high
        sz  = blksz[k-1];
        nnz = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        if((double)2*nnz/(sz*(sz+1)) > 0.75) {
          for(j = 0; j < sz*sz; j++) MM[j] = 0.0;
          for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++) {
            r = CArow[j] - 1;
            c = CAcol[j] - 1;
            MM[c*sz + r] = MM[r*sz + c] = CAent[j];
          }
          dsyev_(&jobz, &uplo, &sz, MM, &sz, w, work, &lwork, &info);
          if(info == 0) {
            maxeval = -1.0e10;
            for(j = 0; j < sz; j++) if(fabs(w[j]) > maxeval) maxeval = fabs(w[j]);
            rank = 0;
            for(j = 0; j < sz; j++) if(fabs(w[j])/maxeval > eps) rank++;
            printf("(h,k) = (%d,%d) : rank %d\n", h, k, rank);
          }
//           else printf("h = %d : eval computation bad\n", h);

        }

        if(rank <= sz/10 && (double)2*nnz/(sz*(sz+1)) > 0.75) {
          fprintf(fid, "%d %d l %d\n", h, k, rank);
          for(j = 0; j < sz; j++) if(fabs(w[j]/maxeval) > eps)
            fprintf(fid, "%.15e\n", w[j]);
          for(j = 0; j < sz; j++) if(fabs(w[j]/maxeval) > eps)
            for(i = 0; i < sz; i++)
              fprintf(fid, "%.15e\n", MM[j*sz + i]);

        }
        else {
        
          // Right now, this matrix is assumed sparse
          j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
          fprintf(fid, "%d %d s %d\n", h, k, j);
          for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
            fprintf(fid, "%d %d %.16e\n", CArow[j], CAcol[j], CAent[j]);

        }

      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'l')
      {
        // Right now, this matrix is assumed dense
        j = (CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ])/(blksz[k-1]+1);
        fprintf(fid, "%d %d l %d\n", h, k, j);
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
          fprintf(fid, "%.16e\n", CAent[j]);
      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'd')
      {
        // Right now this matrix is assumed sparse
        j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        fprintf(fid, "%d %d s %d\n", h, k, j);
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
          fprintf(fid, "%d %d %.16e\n", CArow[j], CAcol[j], CAent[j]);
      }
      else {
        printf("Error (writedata_sdplr): Encountered data matrix not of type 's' or 'l' or 'd'.\n");
        fclose(fid);
        return 0;
      }
    }
  }

  fclose(fid);

  MYFREE(MM);
  MYFREE(w);
  MYFREE(work);

  return 0;

}

int writedata_sdpa(char* datafilename, int m, int numblk, int* blksz,
                   char* blktype, double* b, double* CAent,
                   int* CArow, int* CAcol, int* CAinfo_entptr,
                   int* CAinfo_rowcolptr, char* CAinfo_type,
                   char* CAinfo_storage)
{
  int h, i, j, k;
  FILE *fid;

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_sdpa): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++) {
    if(blktype[k-1] == 's') fprintf(fid, "%d ", blksz[k-1]);
    else if(blktype[k-1] == 'd') fprintf(fid, "%d ", -blksz[k-1]);
  }
  fprintf(fid, "\n");

  for(i = 1; i <= m; i++)
    fprintf(fid, "%.16e  ", b[i-1]);
  fprintf(fid, "\n");

  for(h = 0; h <= m; h++)
  {
    for(k = 1; k <= numblk; k++)
    {
      if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 's')
      {
        // Right now, this matrix is assumed sparse
        j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++) {
          if(h == 0) fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], -CAent[j]);
          else       fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], CAent[j]);
        }
      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'l')
      {
        printf("error: Low-rank matrices not supported in SDPA format.\n");
        exit(0);
      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'd')
      {
        // Right now this matrix is assumed sparse
        j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
          fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], CAent[j]);
      }
      else {
        printf("Error (writedata_sdplr): Encountered data matrix not of type 's' or 'l' or 'd'.\n");
        fclose(fid);
        return 0;
      }
    }
  }

  fclose(fid);

  return 0;

}

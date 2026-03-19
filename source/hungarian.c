//
//(Minimum) Assignment Problem by Hungarian Algorithm
//taken from Knuth's Stanford Graphbase
// https://ranger.uta.edu/~weems/NOTES5311/hungarian.c


#include <stdio.h>
#include <stdlib.h>
#define INF (0x7FFFFFFF)
#define verbose (0)

// 

void hungarian(long long int size1, long long int size2, long long int Array[][size2], char Result[][size2]);


void HUNGMAIN (long long int size1, long long int size2, long long int COST[][size2], long int **MATCH)
{

//long int Array[size1][size2];
//long int (*Array)[size2] = malloc(size1*sizeof(*Array)); 
long long int (*Array)[size2] = COST;
char (*Result)[size2] = malloc(size1*sizeof(*Result));  // used as boolean

long int i;
/*
for (i=0;i<size1;++i)
  for (j=0;j<size2;++j)
    Array[i][j]=999999999999;
*/
long int y,x;
/*
//printf ("I am here %li %li \n", size1, size2);
//
for (y=0;y<size1;y++) {
  for (x=0;x<size2;x++) {
      //printf ("%li ", COST[y][x]);
      Array[y][x]=COST[y][x];
    }
  //printf ("\n");
}
//printf ("I am here\n");
*/
hungarian(size1, size2, Array, Result);
i=0;
for (y=0;y<size1;++y)
  for (x=0;x<size2;++x)
    if (Result[y][x]) {
      //printf("%lld and %lld are connected in the assignment\n",y,x);
      (*MATCH)[i]=x;
      i++;
    }
free (Result);
}

void hungarian(long long int size1, long long int size2, long long int Array[][size2], char Result[][size2])
{
long long int i,j;
int False=0;
int True=1;

unsigned long int m=size1,n=size2;
long long int k;
long long int l;
long long int s;
//long int col_mate[size1];
//long int unchosen_row[size1];
//long int row_dec[size1];
long long int *col_mate = malloc(size1*sizeof (long int));
long long int *unchosen_row = malloc(size1*sizeof (long int));
long long int *row_dec = malloc(size1*sizeof (long int));
for (i=0;i<size1;i++) {
  col_mate[i]=0;
  unchosen_row[i]=0;
  row_dec[i]=0;
}
//long int row_mate[size2];
//long int parent_row[size2];
//long int col_inc[size2];
//long int slack[size2];
//long int slack_row[size2];
long long int *row_mate = malloc(size2*sizeof(long int));
long long int *parent_row= malloc(size2*sizeof(long int));
long long int *col_inc= malloc(size2*sizeof(long int));
long long int *slack= malloc(size2*sizeof(long int));
long long int *slack_row= malloc(size2*sizeof(long int));

for (i=0;i<size2;i++) {
  row_mate[i]=0;
  parent_row[i]=0;
  col_inc[i]=0;
  slack[i]=0;
  slack_row[i]=0;
}
long long int t;
long long int q;
long long int unmatched;
long long int cost=0;

for (i=0;i<size1;++i)
  for (j=0;j<size2;++j)
    Result[i][j]=False;

// Begin subtract column minima in order to start with lots of zeroes 12
//printf("Using heuristic\n");
for (l=0;l<n;l++)
{
  s=Array[0][l];
  for (k=1;k<n;k++)
    if (Array[k][l]<s)
      s=Array[k][l];
  cost+=s;
  if (s!=0)
    for (k=0;k<n;k++)
      Array[k][l]-=s;
}
// End subtract column minima in order to start with lots of zeroes 12

// Begin initial state 16
t=0;
for (l=0;l<n;l++)
{
  row_mate[l]= -1;
  parent_row[l]= -1;
  col_inc[l]=0;
  slack[l]=INF;
}
for (k=0;k<m;k++)
{
  s=Array[k][0];
  for (l=1;l<n;l++)
    if (Array[k][l]<s)
      s=Array[k][l];
  row_dec[k]=s;
  for (l=0;l<n;l++)
    if (s==Array[k][l] && row_mate[l]<0)
    {
      col_mate[k]=l;
      row_mate[l]=k;
      if (verbose)
        printf("matching col %lld==row %lld\n",l,k);
      goto row_done;
    }
  col_mate[k]= -1;
  if (verbose)
    printf("node %lld: unmatched row %lld\n",t,k);
  unchosen_row[t++]=k;
row_done:
  ;
}
// End initial state 16
 
// Begin Hungarian algorithm 18
if (t==0)
  goto done;
unmatched=t;
while (1)
{
  if (verbose)
    printf("Matched %lld rows.\n",m-t);
  q=0;
  while (1)
  {
    while (q<t)
    {
      // Begin explore node q of the forest 19
      {
        k=unchosen_row[q];
        s=row_dec[k];
        for (l=0;l<n;l++)
          if (slack[l])
          {
            long long int del;
            del=Array[k][l]-s+col_inc[l];
            if (del<slack[l])
            {
              if (del==0)
              {
                if (row_mate[l]<0)
                  goto breakthru;
                slack[l]=0;
                parent_row[l]=k;
                if (verbose)
                  printf("node %lld: row %lld==col %lld--row %lld\n",
                    t,row_mate[l],l,k);
                unchosen_row[t++]=row_mate[l];
              }
              else
              {
                slack[l]=del;
                slack_row[l]=k;
              }
          }
        }
      }
      // End explore node q of the forest 19
      q++;
    }
 
    // Begin introduce a new zero into the matrix 21
    s=INF;
    for (l=0;l<n;l++)
      if (slack[l] && slack[l]<s)
        s=slack[l];
    for (q=0;q<t;q++)
      row_dec[unchosen_row[q]]+=s;
    for (l=0;l<n;l++)
      if (slack[l])
      {
        slack[l]-=s;
        if (slack[l]==0)
        {
          // Begin look at a new zero 22
          k=slack_row[l];
          if (verbose)
            printf(
              "Decreasing uncovered elements by %lld produces zero at [%lld,%lld]\n",
              s,k,l);
          if (row_mate[l]<0)
          {
            for (j=l+1;j<n;j++)
              if (slack[j]==0)
                col_inc[j]+=s;
            goto breakthru;
          }
          else
          {
            parent_row[l]=k;
            if (verbose)
              printf("node %lld: row %lld==col %lld--row %lld\n",t,row_mate[l],l,k);
            unchosen_row[t++]=row_mate[l];
          }
          // End look at a new zero 22
        }
      }
      else
        col_inc[l]+=s;
    // End introduce a new zero into the matrix 21
  }
breakthru:
  // Begin update the matching 20
  if (verbose)
    printf("Breakthrough at node %lld of %lld!\n",q,t);
  while (1)
  {
    j=col_mate[k];
    col_mate[k]=l;
    row_mate[l]=k;
    if (verbose)
      printf("rematching col %lld==row %lld\n",l,k);
    if (j<0)
      break;
    k=parent_row[j];
    l=j;
  }
  // End update the matching 20
  if (--unmatched==0)
    goto done;
  // Begin get ready for another stage 17
  t=0;
  for (l=0;l<n;l++)
  {
    parent_row[l]= -1;
    slack[l]=INF;
  }
  for (k=0;k<m;k++)
    if (col_mate[k]<0)
    {
      if (verbose)
        printf("node %lld: unmatched row %lld\n",t,k);
      unchosen_row[t++]=k;
    }
  // End get ready for another stage 17
}
done:

// Begin doublecheck the solution 23
for (k=0;k<m;k++)
  for (l=0;l<n;l++)
    if (Array[k][l]<row_dec[k]-col_inc[l])
      exit(0);
for (k=0;k<m;k++)
{
  l=col_mate[k];
  if (l<0 || Array[k][l]!=row_dec[k]-col_inc[l])
    exit(0);
}
k=0;
for (l=0;l<n;l++)
  if (col_inc[l])
    k++;
if (k>m)
  exit(0);
// End doublecheck the solution 23
// End Hungarian algorithm 18

for (i=0;i<m;++i)
{
  Result[i][col_mate[i]]=True;
 /*TRACE("%lld - %lld\n", i, col_mate[i]);*/
}
for (k=0;k<m;++k)
{
  for (l=0;l<n;++l)
  {
    /*TRACE("%lld ",Array[k][l]-row_dec[k]+col_inc[l]);*/
    Array[k][l]=Array[k][l]-row_dec[k]+col_inc[l];
  }
  /*TRACE("\n");*/
}
for (i=0;i<m;i++)
  cost+=row_dec[i];
for (i=0;i<n;i++)
  cost-=col_inc[i];
//printf("Cost is %lld\n",cost);
free (col_mate);
free (unchosen_row); 
free (row_dec);
free (row_mate); 
free (parent_row); 
free (col_inc); 
free (slack); 
free (slack_row); 

}


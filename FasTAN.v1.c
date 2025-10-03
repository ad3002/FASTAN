/*******************************************************************************************
 *
 *  Convert a .fasta or .1seq file to a .gdb:
 *
 *  Author:   Gene Myers
 *            Modified by Richard Durbin to also read .1seq files
 *  Origin:   Complete rework/specialization from fasta2DAM of DAZZ_DB for FASTGA module
 *  Creation: Jan 2024
 *  Last Mod: Mar 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <dirent.h>

#include "GDB.h"

#define SAT_MAX 3000

static char *Usage = "[-v] <source:path>[<fa_extn>|<1_extn>]";

void analyze_block(uint8 *seq, int len)
{ int    i, p, x;
  int    d, e, f, a;
  uint16 kmer;
  uint32 u;
  uint16 count[0x10000];
  uint16 index[0x10000];
  uint16 antis[0x10000];
  uint16 diags[0x10000];
  uint32 prelm[0x10000];
  uint32 final[0x10000];
  uint8 *s7;
  int    l7;

  for (i = 0; i < 0x10000; i++)          //  Init counters
    antis[i] = count[i] = diags[i] = 0;

  s7 = seq+7;
  l7 = len-7;

  kmer = seq[0];                 //  count # of each 8-mer
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    { kmer = (kmer << 2) | s7[i];
      count[kmer] += 1;
    }

  p = 0;                         //  turn counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  kmer = seq[0];                 //  place positions in index in order of 8-mer
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    { kmer = (kmer << 2) | s7[i];
      index[count[kmer]++] = i;
    }

  p = 0;                          //  examine pos-pairs within each 8-mer bucket
  for (i = 0; i < 0x10000; i++)   //    and count anti-diag's of each pair within
    { x = count[i];               //    band of width MAX_SAT
      if (p == x)
        continue;
      f = index[p];
      for (p++; p < x; p++)
        { e = index[p];
          d = e-f;
          if (d < 3000)
            antis[(f+e)>>1] += 1; // push d,(f+e)>>1
          f = e;
        }
    }

  p = 0;                          //  turn anti-d counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = antis[i];
      antis[i] = p;
      p += x;
    }

  p = 0;                         // sweep pos-pairs again this time placing
  for (i = 0; i < 0x10000; i++)  //   anti-d + diag in prelim in order of anti-d.
    { x = count[i];              //   count diag's.
      if (p == x)
        continue;
      f = index[p];
      for (p++; p < x; p++)
        { e = index[p];
          d = e-f;
          if (d < SAT_MAX)
            { a = ((f+e)>>1);
              prelm[antis[a]++] = (a << 16) | d;
              diags[d>>3] += 1;
            }
          f = e;
        }
    }

  p = 0;                        //  turn diag counts into sort ptrs
  for (i = 0; i < SAT_MAX; i++)
    { x = diags[i];
      diags[i] = p;
      p += x;
    }

  for (i = 0; i < p; i++)       //  place anti-d/d pairs in final in order of
    { u = prelm[i];  //    diagonal
      final[diags[(u&0xffff)>>3]++] = u;
    }

  // for (i = 0; i < p; i++)
    // printf(" %4d : %5d\n",final[i]&0xffff,final[i]>>16);

  printf("HITS = %d\n",p);
}


int main(int argc, char *argv[])
{ char *spath, *tpath;
  int   ftype;
  GDB  _gdb, *gdb = &_gdb;

  int VERBOSE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("FasTAN")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  if (argc == 2)
    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,1);
  else
    ftype = Get_GDB_Paths(argv[1],argv[2],&spath,&tpath,1);

  Create_GDB(gdb,spath,ftype,1,NULL);

  free(tpath);
  free(spath);

  { uint8 *buffer;
    int    clen;
    int    i, p;

    buffer = (uint8 *) malloc(gdb->maxctg + 4);
    for (i = 0; i < gdb->ncontig; i++)
      { Get_Contig(gdb,i,NUMERIC,(char *) buffer);
        clen = gdb->contigs[i].clen;
        if (clen < 0x10000)
          analyze_block(buffer,clen);
        else
          for (p = 0; p+0x8000 <= clen; p += 0x8000)
            if (p+0x10000 > clen)
              analyze_block(buffer+p,clen-p);
            else
              analyze_block(buffer+p,0x10000);
      }
  }

  exit (0);
}

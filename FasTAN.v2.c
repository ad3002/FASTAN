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
#include "align.h"

#define SAT_MAX 3000

#undef  SHOW_SEEDS
#define SHOW_CHAINS

static char *Usage = "[-v] <source:path>[<fa_extn>|<1_extn>]";

void Print_Seq(uint8 *seq)
{ static char dna[4] = { 'a', 'c', 'g', 't' };
  int j;

  for (j = 0; j < 8; j++)
    printf("%c",dna[seq[j]]);
}

typedef struct
  { void       *block;
    Work_Data  *work;
    Align_Spec *spec;
    Alignment  *align;
  } Work_Packet;

typedef struct
  { uint16  diag;
    uint16  anti;
  } Seed;

void filter_block(uint8 *seq, int off, int len, Work_Packet *packet)
{ void       *block = packet->block;
  Work_Data  *work  = packet->work;
  Align_Spec *spec  = packet->spec;
  Alignment  *align = packet->align;
  Path       *path  = align->path;

  int    i, p, x;
  int    d, e, f, a;
  uint16 kmer;
  uint16 *count;  // 0x10000
  uint16 *index;  // 0x08000
  uint16 *diags;  // 0x08000  (less, actually "sm")
  Seed   *post;   // 0x08000
  Seed   *hits;   // 0x08000
  int     nhits;
  uint8 *s7;
  int    l7, sm;

  count = (uint16 *) block;
  index = count + 0x10000;
  diags = index +  0x8000; 
  post  = (Seed *) (diags + 0x8000);
  hits  = (Seed *) count;

  for (i = 0; i < 0x10000; i++)          //  Init counters
    count[i] = 0;

  s7 = seq+7;
  l7 = len-7;
  sm = ((SAT_MAX-1) >> 3) + 1;

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
  index[l7] = 0;

  index[0] |= 0x8000;
  for (i = 0; i < 0xffff; i++)   //  mark bucket ends and reset count
    { index[count[i]] |= 0x8000;
      count[i] = 0;
    }
  count[0xffff] = 0;

  e = (index[0] & 0x7fff);        //  count anti-diagonals of adjacent pairs in a bucket
  for (i = 1; i < l7; i++)
    { f = index[i];
      if (f < 0x8000)
        { if (f-e < SAT_MAX)
            count[f+e] += 1;
          e = f;
        }
      else
        e = (f & 0x7fff);
    }

  p = 0;                          //  turn counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  for (i = 0; i < sm; i++)         //  Init counters
    diags[i] = 0;

  e = (index[0] & 0x7fff);         //  place diag/anti-d hits in post in order of anti-d
  for (i = 1; i < l7; i++)
    { f = index[i];
      if (f < 0x8000)
        { d = f-e;
          if (d < SAT_MAX)
            { a = e+f;
              x = count[a]++;
              post[x].anti = a;
              post[x].diag = d;
              diags[d>>3] += 1;
            }
          e = f;
        }
      else
        e = (f & 0x7fff);
    }

  p = 0;                        //  turn diag counts into sort ptrs
  for (i = 0; i < sm; i++)
    { x = diags[i];
      diags[i] = p;
      p += x;
    }
  nhits = p;

  for (i = 0; i < nhits; i++)       //  place diag/anti-d pairs in hits in lex order
    { x = (post[i].diag >> 3);
      hits[diags[x]++] = post[i];
    }
  hits[nhits].diag = 0;
  hits[nhits].anti = 0;

  printf("\n\nPanel %d-%d  ::  Hiths = %d\n",off,off+len,nhits); fflush(stdout);

#ifdef SHOW_CHAINS

  { int p1, p2;
    int t1, t2;
    int a1, a2;
    int d, a;
    int cap;
    int cbeg, cend, cov, pure;
    int   dgmin, dgmax, low, hgh;
    int64 dave, dnum;
#ifdef CHAIN
    int p;
    int s1, s2;
#endif

    t1 = diags[1];
    t2 = diags[2];
    for (i = 3; i <= sm; i++)
      { p1 = t1;
        t1 = p2 = t2;
        if (i == sm)
          t2 = t1;
        else
          t2 = diags[i];
        cap  = t1+t2;
        a1 = hits[p1].anti;
        a2 = hits[p2].anti;
        cend = -129;
        cbeg = -129;
        cov  = 0;
        while (p1+p2 < cap)
          { if (p2 >= t2 || a1 < a2) 
              { d = hits[p1].diag;
                a = a1;
              }
            else
              { d = hits[p2].diag;
                a = a2;
              }

            if (a > cend+128)
              { if (dnum*(cend-cbeg) >= 2*dave && cov > .5*(cend-cbeg) && !pure)
                  { dave = dave/dnum;
                    printf("\nChain @%5lld len = %4d cov = %4d in bands %d - %d  ::  %lld\n",
                           off + (cbeg-dave)/2,(cend-cbeg)/2,cov/2,8*(i-1),8*(i+1),dave);
#ifdef CHAIN
                    for (p = s1; p < p1; p++)
                      printf("   1(%d)  %4d %8d  (%8d,%8d)\n",p,hits[p].diag,hits[p].anti,
                             off+(hits[p].diag+hits[p].anti)/2,off+(hits[p].anti-hits[p].diag)/2);
                    for (p = s2; p < p2; p++)
                      printf("   2(%d)  %4d %8d  (%8d,%8d)\n",p,hits[p].diag,hits[p].anti,
                             off+(hits[p].diag+hits[p].anti)/2,off+(hits[p].anti-hits[p].diag)/2);
#endif

                    dgmin = dave-8;
                    dgmax = dave+8;
                    low   = dave/2;
                    hgh   = dave/2;
                    if (dgmin <= 0)
                      dgmin = 1;
                    if (low >= dgmin)
                      low = dgmin-1;
                    Local_Alignment(align,work,spec,dgmin,dgmax,2*off+(cbeg+cend)/2,low,hgh);
                    if (path->aepos - path->abpos < .45*(cend-cbeg))
                      printf("  Fail\n");
                    else
                      { printf("  Local %d-%d vs %d-%d (%d)\n",
                               path->abpos,path->aepos,path->bbpos,path->bepos,path->diffs);
                        printf("  Beyond tips: %d< >%d\n",
                               (cbeg - (path->abpos+path->bbpos-2*off))/2,
                               ((path->aepos+path->bepos-2*off) - cend)/2);
                        Compute_Trace_PTS(align,work,100,GREEDIEST,dave-low,dave-hgh);
                        Print_Alignment(stdout,align,work,4,100,10,0,8);
                      }
                    fflush(stdout);
                  }
                cbeg = a;
                cov  = 16;
                pure = 1;
                dave = d;
                dnum = 1;
#ifdef CHAIN
                s1 = p1;
                s2 = p2;
#endif
              }
            else
              { if (a < cend)
                  cov += (a+16) - cend;
                else
                  cov += 16;
                dave += d;
                dnum += 1;
              }
            cend = a + 16;

            if (p2 >= t2 || a1 < a2) 
              { a1 = hits[++p1].anti;
                pure = 0;
              }
            else
              a2 = hits[++p2].anti;
          }
      }
  }

#endif

#ifdef SHOW_SEEDS
  p = 0;
  for (i = 0; i < sm; i++)
    { printf("Band %d\n",i);
      for (e = p; e < diags[i]; e++)
        { d = hits[e].diag;
          a = hits[e].anti;
          f = ((a+d) >> 1);
          printf("   %4d : %5d  ",d,a);
          Print_Seq(seq+f);
          printf("\n");
        }
      p = diags[i];
    }
#endif
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

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,1);

  Create_GDB(gdb,spath,ftype,1,NULL);

  free(tpath);
  free(spath);

  { uint8 *buffer;
    void  *block;
    int    clen;
    int    i, p;
    Work_Packet pack;
    Work_Data  *work;
    Align_Spec *spec;
    Alignment  _align, *align = &_align;
    Overlap    _ovl, *ovl = &_ovl;

    block  = malloc(0x60000);
    buffer = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;

    work = New_Work_Data();
    spec = New_Align_Spec(.75,100,gdb->freq,0);

    if (block == NULL || buffer == NULL || work == NULL || spec == NULL)
      { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
        exit (1);
      }

    align->aseq  = (char *) buffer;
    align->bseq  = (char *) buffer;
    align->flags = ovl->flags = 0;
    align->path  = &(ovl->path);

    pack.block = block;
    pack.work  = work;
    pack.spec  = spec;
    pack.align = align;

    for (i = 0; i < gdb->ncontig; i++)
      { Get_Contig(gdb,i,NUMERIC,(char *) buffer);
        clen = gdb->contigs[i].clen;
        align->alen = align->blen = clen;
        ovl->aread = ovl->bread = i;
        if (clen < 0x8000)
          filter_block(buffer,0,clen,&pack);
        else
          for (p = 0; p+0x2000 <= clen; p += 0x6000)
            if (p+0x8000 > clen)
              filter_block(buffer+p,p,clen-p,&pack);
            else
              filter_block(buffer+p,p,0x8000,&pack);
      }

    Free_Align_Spec(spec);
    Free_Work_Data(work);
    free(buffer-1);
    free(block);
  }

  exit (0);
}

/*******************************************************************************************
 *
 *  Search an assembly for satellitic repeats
 *
 *  Author:   Gene Myers
 *  Creation: Jan 2024
 *  Last Mod: July 2025
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
#include "alncode.h"

#define TSPACE   100
#define VERSION "0.1"

static char *Usage = "[-v] <source:path>[<fa_extn>|<1_extn>] <target>[.1aln]";

static char dna[4] = { 'a', 'c', 'g', 't' };

static OneFile    *Ofile;
static Work_Data  *Work;
static Align_Spec *Spec;
static Alignment  *Align;
static Overlap    *Over;
static int64      *Trace64;
static int         TraceMax;

static void Print_Seq(uint8 *seq, int len)
{ int j;

  for (j = 0; j < len; j++)
    printf("%c",dna[seq[j]]);
}


/*******************************************************************************************
 *
 *  DEBRUIJN BASED DETECTOR
 *
 ********************************************************************************************/

#undef   SHOW_NODE_SELECTION
#undef   SHOW_EDGE_SELECTION
#undef   SHOW_GRAPH
#undef   SHOW_UGRAPH
#undef   SHOW_CYCLE
#undef   SHOW_ALIGNMENTS
#define  GREEDY

#ifdef GREEDY

#define FTHRESH 10

#else

#define FTHRESH 30

#endif

typedef struct
  { uint16  kmer;
    int16   out;
    int16   status;
    int16   freq;
  } Dnode;

#define KMER(x) node[x].kmer
#define EDGE(x) node[x].out
#define STAT(x) node[x].status
#define FREQ(x) node[x].freq

static char *emer(int x)
{ static char mer[9];
  int i;

  mer[8] = '\0';
  for (i = 7; i >= 0; i--)
    { mer[i] = dna[x&0x3];
      x >>= 2;
    }
  return (mer);
}

static int show_edge(int e, Dnode *node, int16 *edge, int16 *path)
{ int p, x, n;

  p = path[e];
  n = 1;
  for (x = edge[e]; x != p; x = edge[EDGE(x)])
    { n += 1;
#ifdef SHOW_EDGE
      printf("%c",dna[KMER(x)&0x3]);
#endif
    }
#ifdef SHOW_EDGE
  printf("%c",dna[KMER(x)&0x3]);
#endif
  return (n);
}

static void show_segments(int i, Dnode *node, int16 *edge, int16 *path, int16 *pwgt)
{ int e, p, x;
  int clen;

  printf("\n  %d:%s(%d)\n",i,emer(KMER(i)),FREQ(i));
  for (e = EDGE(i); e < EDGE(i+1); e++)
    { p = path[e];
      if (p < 0)
        continue;
      printf("  -> ");
      clen = 1;
      for (x = edge[e]; x != p; x = edge[EDGE(x)])
        { printf("%c",dna[KMER(x)&0x3]);
          clen += 1;
        }
      printf("%c",dna[KMER(x)&0x3]);
      printf("(%d:%d)",clen,pwgt[e]);
      printf(" -> %d:%s(%d)\n",x,emer(KMER(x)),FREQ(x));
    }
}

static Dnode  *Node;
static int16  *Edge;
static int16  *Pdst;
static int16  *Pwgt;

static int     Id;
static int16  *Low;
static int16  *Num;

void processed(int v)
{ int e, u;

  Node[v].status = 3;
  for (e = Node[v].out; e < Node[v+1].out; e++)
    { u = Pdst[e];
      if (u >= 0 && Node[u].status != 3)
        processed(u);
    }
}

void strong_component(int v, int f)
{ int e, u;

  Low[v] = Num[v] = Id++;
  Node[v].status = 2;
  for (e = Node[v].out; e < Node[v+1].out; e++)
    { u = Pdst[e];
      if (u < 0)
        continue;
      if (Num[u] == 0)
        { strong_component(u,e);
          if (Low[u] < Low[v])
            Low[v] = Low[u];
        }
      else if (Node[u].status == 2)
        { if (Num[u] < Low[v])
            Low[v] = Num[u];
        }
      else
        Pdst[e] = -1;
    }

  if (Low[v] == Num[v])
    { processed(v);
      Node[v].status = 4;
      if (f >= 0)
        Pdst[f] = -1;
    }
}

void show_scc(int v, int m)
{ int e, u;

  Node[v].status = m;
  show_segments(v,Node,Edge,Pdst,Pwgt);
  for (e = Node[v].out; e < Node[v+1].out; e++)
    { u = Pdst[e];
      if (u >= 0 && Node[u].status != m)
        show_scc(u,m);
    }
}

typedef struct
  { int next;
    int node;
  } Elem;

static int     Start;
static int16   Free;
static int     Ncycles;
static int16  *List;
static int16  *Blocked;
static int16  *Blist;
static Elem   *Link;

void list_scc(int v, int m)
{ int e, u;

  Node[v].status = m;
  List[Id++] = v;
  for (e = Node[v].out; e < Node[v+1].out; e++)
    { u = Pdst[e];
      if (u >= 0 && Node[u].status != m)
        list_scc(u,m);
    }
}

static inline int16 newElem(int v, int n)
{ int r;

  r = Free;
  Free = Link[Free].next;
  Link[r].next = n;
  Link[r].node = v;
  return (r);
}

static inline int16 freeElem(int r)
{ int n;

  n = Link[r].next;
  Link[r].next = Free;
  Free = r;
  return (n);
}

static void unblock(int u)
{ int c, w;

  Blocked[u] = 0;
  for (c = Blist[u]; c >= 0; c = freeElem(c))
    { w = Link[c].node;
      if (Blocked[w])
        unblock(w);
    }
  Blist[u] = -1;
}

static int circuit(int v)
{ int f, e, w;

  f = 0;
  Blocked[v] = 1;
  for (e = Node[v].out; e < Node[v+1].out; e++)
    { w = Pdst[e];
      if (w < 0)
        continue;
      else if (w == Start)
        { int n, s;

          printf("  Cycle:");
          printf(" %d",Start);
          for (f = 0; f < Id; f++)
            printf(" %d",Pdst[List[f]]);
          printf(" %d\n    ",Start);

          n = 0;
          s = Pwgt[e];
          for (f = 0; f < Id; f++)
            { n += show_edge(List[f],Node,Edge,Pdst);
              if (s > Pwgt[List[f]])
                s = Pwgt[List[f]];
            }
          n += show_edge(e,Node,Edge,Pdst);
          printf("[%d]::%d\n",n,s);
   
          Ncycles += 1;
          f = 1;
        }
      else if (Blocked[w] == 0)
        { List[Id++] = e;
          if (circuit(w))
            f = 1;
          Id -= 1;
        }
      if (Ncycles >= 500)
        { f = 1;
          break;
        }
    }
  if (f)
    unblock(v);
  else
    for (e = Node[v].out; e < Node[v+1].out; e++)
      { int r;

        w = Pdst[e];
        if (w < 0)
          continue;
        for (r = Blist[w]; r >= 0; r = Link[r].next)
          if (Link[r].node == v)
            break;
        if (r < 0)
          Blist[w] = newElem(v,Blist[w]);
      }
  return (f);
}

static int greedy_cycle(int v)
{ int e, best, whch;

  if (Blocked[v] == 0)
    { Blocked[v] = 1;
      best = 0;
      for (e = Node[v].out; e < Node[v+1].out; e++) 
        { if (Pdst[e] >= 0 && Pwgt[e] > best)
            { best = Pwgt[e];
              whch = e;
            }
        }
      Low[Id++] = whch;
      e = greedy_cycle(Pdst[whch]);
      Blocked[v] = 0;
      return (e);
    }
  else
    return (v);
}

/*

typedef struct
  { int lab;
    int elist;
  } POAnode;

typedef struct
  { int enext;
    int node;
  } POAedge;

static void add_to_poa(uint8 *seq, int len)
{
  qstop = qtop;
  while (qbot < qstop)
    { i = queue[qbot].spx;
      v = queue[qbot].vtx;
      qbot += 1;
      for (e = POA[v].elist; e != -1; e = POE[e].enext)
        { w = POE[e].node;
          if (! (i,w) marked)
            { do
                { deg2 = 0;
                  for (f = POA[w].elist; f != -1; f = POE[f].enext)
                    { deg2 += 1;
                      z = POE[f].node;
                      if (POA[z].lab == seq[i])
                        { if (deg2 >= 2 || POE[f].enext >= 0)
                            { queue[qtop].spx = i;
                              queue[qtop].vtx = w;
                              queue[qtop].mark = mark[w];
                              mark[w] = qtop;
                              qtop += 1;
                            }
                          w.mark[i] = 1;
                          w = z;
                          i = i+1;
                          break;
                        }
                    }
                }
              while (f >= 0);
              queue[qtop].spx = i;
              queue[qtop].vtx = w;
              queue[qtop].mark = mark[w];
              mark[w] = qtop;
              qtop += 1;
            }
          i += 1;
          while exit z in A(v) s.t. z.lab == a[i+1])
            { push (i,v) if A(v) > 1 onto level k+1
              v.mark[i] = 1;
              v = z;
              i += 1;
            }
          push (i,v) onto level k+1
    }
}

*/

static void cut_trace(Alignment *align)
{ Path  *path;
  int   *trace, tlen;
  int    p, c, x;
  int    bp, bn, ap, ln;
  uint8 *seq;

  path  = align->path;
  trace = (int *) (path->trace);
  tlen  = path->tlen;
  seq   = (uint8 *) align->aseq;

  p = 0;
  ap = path->abpos;
  bp = path->bbpos;
  ln = ap-bp;
  while (bp < path->bepos)
    { bn = ap;
      printf("%d-%d: ",bp,ap);
      Print_Seq(seq+bp,ap-bp);
      printf("\n");
      while (p < tlen)
        { c = trace[p];
          if (c < 0)
            { c = -c;
              c -= 1;
              x = c-ap; 
              if (bp+(x+1) > bn)
                break;
              bp += x+1;
              ap  = c;
            }
          else
            { c -= 1;
              x = c-bp;
              if (ap+(x+1) > bn+ln)
                break;
              ap += x+1;
              bp = c;
            }
          p += 1;
        }
      ap += (bn-bp);
      bp  = bn;
    }
  printf("%d-%d: ",bn,path->aepos);
  Print_Seq(seq+bn,path->aepos-bn);
  printf("\n");
  fflush(stdout);
}

typedef struct
  { int  start;
    int  cycle;
    int  first;
  } Hit;

static int NSORT(const void *l, const void *r)
{ Hit *x = (Hit *) l;
  Hit *y = (Hit *) r;

  if (x->cycle < y->cycle)
    return (-1);
  else if (x->cycle > y->cycle)
    return (1);
  else
    return (x->first - y->first);
}

static int euler_block(uint8 *seq, int off, int len, void *block)
{ int     nlen, elen;
  double  freq[4];
  int16  *count;  // 64K + 1
  int16  *map;    // 64K
  int16  *index;  // 32K
  Dnode  *node;   // 3.2K + 1
  int16  *edge;   // 3.2K * 4
  int16  *ewgt;   // 3.2K * 4
  int16  *path;   // 3.2K * 4
  int16  *pwgt;   // 3.2K * 4
  int16  *low;    // 3.2K
  int16  *num;    // 3.2K
  int16  *list;   // 3.2K
  int16  *bflag;  // 3.2K
  int16  *blist;  // 3.2K
  Elem   *link;   // 3.2K * 4

  Hit *seed;
  int  nseed, maxseed;

  count = (int16 *) block;         //  64K * 5 + 3.2K * (sizeof(Dnode) + 16)
  index = count + 0x10004;
  node  = (Dnode *) (index + 0x8000);
  edge  = (int16 *) (node + (0x8000/10) + 1);
  ewgt  = (int16 *) (edge + (0x20000/10));
  map   = (int16 *) (ewgt + (0x20000/10));

  path  = map;
  pwgt  = path + (0x20000/10);
  num   = pwgt + (0x20000/10);
  low   = num  + (0x8000/10);
  list  = low  + (0x8000/10);
  bflag = list  + (0x8000/10);
  blist = bflag  + (0x8000/10);
  link  = (Elem *) (blist + (0x8000/10));

  printf("\nPANEL %d-%d\n",off,off+0x8000);

  //  Build index of 8-mers in count and index

  { uint8 *s7;
    int    i, p, x, l7, fq[4];
    uint16 kmer;

    s7 = seq+7;
    l7 = len-7;

    for (i = 0; i < 0x10000; i++)
      count[i] = 0;

    for (i = 0; i < 4; i++)
      fq[i] = 0;

    kmer = seq[0];                //  Count # of each 8-mer
    freq[seq[0]] += 1;
    for (i = 1; i < 7; i++)
      { kmer = (kmer << 2) | seq[i];
        fq[seq[i]] += 1;
      }
    for (i = 0; i < l7; i++)
      { x = s7[i];
        kmer = (kmer << 2) | x;
        count[kmer] += 1;
        fq[x] += 1;
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

    for (i = 0x10000; i > 0; i--)
      count[i] = count[i-1];
    count[0] = 0;

    for (i = 0; i < 4; i++)
      freq[i] = (1.*fq[i])/0x8000;
  }

  //  Make nodes for kmers with count >= FTHRESH(>= 10) and >= random flow in both directions

  { int i, c, p;
    int x, y, n;

#ifdef SHOW_NODE_SELECTION
    printf("\nNODE SELECTION:\n");
    for (i = 0; i < 4; i++)
      printf("  %c: %.4f\n",dna[i],freq[i]);
#endif
    nlen = 0;
    for (i = 0; i < 0x10000; i++)
      { map[i] = -1;
        if ((n = count[i+1]-count[i]) >= FTHRESH)
          { x = 0;
            for (p = 0; p < 4; p++)
              { c = (i>>2) | (p*0x4000);
#ifdef SHOW_NODE_SELECTION
                printf("     < %s : %5d (%5d)\n",emer(c),count[c+1]-count[c],
                                                 (int) ((count[c+1]-count[c])*freq[p]));
#endif
                x += (count[c+1]-count[c]) * freq[p];
              }
            y = 0;
            for (p = 0; p < 4; p++)
              { c = ((i<<2) | p) & 0xffff;
#ifdef SHOW_NODE_SELECTION
                printf("     > %s : %5d (%5d)\n",emer(c),count[c+1]-count[c],
                                                 (int) ((count[c+1]-count[c])*freq[p]));
#endif
                y += (count[c+1]-count[c]) * freq[p];
              }

            if (x <= n && y <= n)
              { KMER(nlen) = i;
                FREQ(nlen) = n;
                STAT(nlen) = 0;
                map[i] = nlen++;
#ifdef SHOW_NODE_SELECTION
                printf("  %s : %5d\n",emer(i),n);
#endif
              }
#ifdef SHOW_NODE_SELECTION
            else
              printf("  %s : %5d XXX\n",emer(i),n);
#endif
          }
      }
  }

  //  Determine positional deBuijn edges amongst nodes
  //    and recognize junction nodes (STAT = 1)

  { int    i, e, c;
    uint16 n, m;
    int    p, q, s, t, h;

#ifdef SHOW_EDGE_SELECTION
    printf("\nEDGE SELECTION:\n");
#endif
    elen = 0;
    for (i = 0; i < nlen; i++)
      { n = KMER(i);
        EDGE(i) = elen;
        for (e = 0; e < 4; e++)
          { m = (n<<2) | e;
            c = map[m];
            if (c >= 0)
              { p = count[n];
                s = count[n+1];
                q = count[m];
                t = count[m+1];
#ifdef SHOW_EDGE_SELECTION
                printf(" %s ->",emer(n));
                printf("%s : ",emer(m));
#endif
                h = 0;
                while (p < s && q < t)
                  { if (index[p] >= index[q])
                      q += 1;
                    else if (index[p] == index[q]-1)
                      { h += 1;
                        p += 1;
                        q += 1;
                      }
                    else
                      p += 1;
                  }

#ifdef SHOW_EDGE_SELECTION
                printf(" %d of %d + %d",h,s-count[n],t-count[m]); 
#endif

                if (h > .3*FREQ(i) || h > .3*FREQ(c))
                  { edge[elen] = c;
                    ewgt[elen] = h;
                    elen += 1;
                    STAT(c) += 1;
                  }
#ifdef SHOW_EDGE_SELECTION
                else
                  printf(" XXX");
                printf("\n");
#endif
              }
          }
      }
    EDGE(nlen) = elen;

    for (i = 0; i < nlen; i++)
      if (STAT(i) != 1 || EDGE(i+1)-EDGE(i) > 1)
        STAT(i) = 1;
      else
        STAT(i) = 0;

#ifdef SHOW_GRAPH
    printf("\nRAW GRAPH:\n");
    for (i = 0; i < nlen; i++)
      { n = KMER(i);
        printf("   %5d : %s, %5d, %d",i,emer(n),FREQ(i),STAT(i));
        for (p = EDGE(i); p < EDGE(i+1); p++)
          printf(" ->%d(%d)",edge[p],ewgt[p]); 
        printf("\n");
      }
#endif
  }

  printf("\nThe uncollapsed graph has %d nodes and %d edges\n",nlen,elen/2);

  //  Collapse unique chains, identify unique cycles (if any)
  //    node[i].status = 1 if core with outedges
  //                     2 if collapsible (in = out = 1) or core with no outedges

  { int i, x, e, s;
    int ncore, nedge;

    for (i = 0; i < nlen; i++)
      if (STAT(i) == 1)
        { for (e = EDGE(i); e < EDGE(i+1); e++)
            { x = edge[e];
              s = ewgt[e];
              while (EDGE(x+1) - EDGE(x) == 1 && STAT(x) == 0)
                { STAT(x) = 2;
                  if (ewgt[EDGE(x)] < s)
                    s = ewgt[EDGE(x)];
                  x = edge[EDGE(x)];
                }
              if (EDGE(x+1) == EDGE(x))
                STAT(x) = 2;

              path[e] = x;
              pwgt[e] = s;
            }
        }

    //  All unmarked nodes are in unique cycles, mark them

    for (i = 0; i < nlen; i++)
      if (STAT(i) == 0)
        { STAT(i) = 1;
          s = ewgt[EDGE(i)];
          for (x = edge[EDGE(i)]; x != i; x = edge[EDGE(x)])
            { STAT(x) = 2;
              if (ewgt[EDGE(x)] < s)
                s = ewgt[EDGE(x)];
            }
          path[EDGE(i)] = i;
          pwgt[EDGE(i)] = s;
        }

      nedge = ncore = 0;
      for (i = 0; i < nlen; i++)
        if (node[i].status != 2)
          { ncore += 1;
            nedge += EDGE(i+1)-EDGE(i);
          }
      printf("\nCollapsed graph has %d nodes, %d edges\n\n",ncore,nedge);
  }

  //  Show the collapsed graph

#ifdef SHOW_UGRAPH
  { int i;

    printf("\nCOLLAPSED GRAPH:\n");
    for (i = 0; i < nlen; i++)
      if (STAT(i) == 1)
        show_segments(i,node,edge,path,pwgt);
  }
#endif

  //  Tarjan's SCC and Johnson's simple cycle algorithms require global access

  Node = node;
  Edge = edge;
  Pdst = path;
  Pwgt = pwgt;
  Low  = low;
  Num  = num;

  //  Find strongly connected components, 1st time
  //    node[i].status = 4 if DFS root of an SCC
  //                     0 if an SCC with no edges || internal to a path
  //                     2 if core and in SCC with edges but not root
  //  

  { int i;

    for (i = 0; i < nlen; i++)
      { num[i] = 0;
        if (STAT(i) != 1)
          STAT(i) = 0;
      }

    Id = 1;
    for (i = 0; i < nlen; i++)
      if (STAT(i) == 1)
        strong_component(i,-1);
  }

  //   Johnson's algorithm for listing all simple cycles
  //     or Greedy algorithm for heaviest cycle.

  { int i, p, e, w, s;
    int loop, mems;

    Blocked = bflag;
    Blist   = blist;
    List    = list;
    Link    = link;

    for (i = 0; i < nlen; i++)
      { Blocked[i] = 0;
        Blist[i]   = -1;
      }

    for (i = 1; i < elen; i++)
      Link[i].next = i-1;
    Link[0].next = -1;
    Free = elen-1;

    nseed   = 0;
    maxseed = 100;
    seed    = malloc(sizeof(Hit)*maxseed);

    for (i = 0; i < nlen; i++)
      if (STAT(i) == 4)
        { loop = -2;
          for (e = EDGE(i); e < EDGE(i+1); e++)
            if (Pdst[e] >= 0)
              { if (loop == -2)
                  loop = e;
                else
                  loop = -1;
              }
          if (loop <= -2)
            { STAT(i) = 0;
              // printf("  %d:%s(%d)   Eliminated\n",i,emer(KMER(i)),FREQ(i));
              continue;
            }

#ifdef SHOW_CYCLE
          printf("\nSCC:\n");
          show_scc(i,2);
          STAT(i) = 4;
          printf("\n");
#endif

          if (loop >= 0 && Pdst[loop] == i)
            { int n;

#ifdef SHOW_CYCLES
              printf("  Cycle: %d %d\n    ",i,i);
#endif
              n = show_edge(loop,node,edge,path);
#ifdef SHOW_CYCLES
              printf("[%d]::%d\n\n",n,pwgt[EDGE(i)]);
#endif

if (nseed >= maxseed)
  printf("Seed Overflow\n");
fflush(stdout);
              seed[nseed].start = i;
              seed[nseed].cycle = n;
              seed[nseed].first = index[count[KMER(i)]];
              nseed += 1;
            }
          else
#ifdef GREEDY
            { Id = 0;
              list_scc(i,1);   //  list[0..Id-1] = all core vertices in SCC(i)
              mems = Id;

              Id = 0;
              Start = greedy_cycle(i);

              for (p = EDGE(Start); p < EDGE(Start+1); p++)
                if (p == Low[0])
                  break;
              if (p >= EDGE(Start+1))
                { for (w = 0; w < Id; w++)
                    if (Pdst[Low[w]] == Start)
                      break;
                  w += 1;
                }
	      else
                w = 0;

#ifdef SHOW_CYCLE
              printf("  Cycle:");
              printf(" %d",Start);
              for (p = w; p < Id; p++)
                printf(" %d",Pdst[Low[p]]);
              printf("\n    ");
#endif
 
              e = 0;
              s = 0x10000;
              for (p = w; p < Id; p++)
                { e += show_edge(Low[p],Node,Edge,Pdst);
                  if (s > Pwgt[Low[p]])
                    s = Pwgt[Low[p]];
                }
#ifdef SHOW_CYCLE
              printf("[%d]::%d\n\n",e,s);
#endif

if (nseed >= maxseed)
  printf("Seed Overflow\n");
fflush(stdout);
              seed[nseed].start = Start;
              seed[nseed].cycle = e;
              seed[nseed].first = index[count[KMER(Start)]];
              nseed += 1;

              for (p = w; p < Id; p++)
                Pdst[Low[p]] = -1;

              Id = 1;
              for (p = 0; p < mems; p++)
                { w = List[p];
                  if (STAT(w) == 1)
                    strong_component(w,-1);
                }
            }
#else
            { Id = 0;
              Start = i;
              Ncycles = 0;
              if (circuit(i) == 0)
                 printf("CIRCUIT FAIL XXX\n");    //  Find cycles of SCC rooted at i
              printf("\n");
              if (Ncycles >= 500)
                printf("TOO MANY CIRCUITS XXX\n");

              Id = 0;
              list_scc(i,1);   //  list[0..Id-1] = all core vertices in SCC(i)
              mems = Id;

              //  Remove edges into i and STAT(i) = 0

              for (p = 0; p < mems; p++)
                { w = list[p];
                  for (e = EDGE(w); e < EDGE(w+1); e++)
                    if (path[e] == i)
                      path[e] = -1;
                  num[w] = 0;
                }
              STAT(i) = 0;

              //  All SCC's of SCC(i)-i are in DFS trees rooted at w st. i->w
              //  Find all SCC's of SCC(i)-i and remove edges out of i

              Id = 1;
              for (e = EDGE(i); e < EDGE(i+1); e++)
                { w = path[e];
                  if (w >= 0)
                    { path[e] = -1;
                      if (STAT(w) == 1)
                        strong_component(w,-1);
                    }
                }
            }
#endif
       }
  }

  { typedef struct
      { int beg;
        int end;
        int bc;
        int ec;
      } SatHit;

    int   len, str;
    int   last, wide, anti;
    int   outhit;
    Path *bpath;
    int   n, e;
    int   nsat, maxsat;
    SatHit *sat;

    qsort(seed,nseed,sizeof(Hit),NSORT);

    nsat   = 0;
    maxsat = 1000;
    sat    = malloc(sizeof(SatHit)*maxsat);

    outhit = 0;
    for (n = 0; n < nseed; n++)
      { len = seed[n].cycle;
        str = seed[n].start;
        if (len == 1)
          continue;

        wide = .2 * len;
        if (wide == 0)
          wide = 1;
#ifdef SHOW_ALIGNMENTS
        printf("  %s %d-%d-%d %d\n",emer(KMER(str)),len-wide,len,len+wide,seed[n].first);
#endif
        last = -1;
        for (e = count[KMER(str)]; e < count[KMER(str)+1]; e++)
          { if (index[e] <= last)
              continue;
            anti = 2*(off + index[e]) + len;
            bpath = Local_Alignment(Align,Work,Spec,len,len,anti,wide,wide);
            if (bpath != NULL && bpath->bepos-bpath->abpos > 20)
              { int beg, end;
                int f, g, s;

                if (Over->path.tlen > TraceMax)
                  { TraceMax = 1.2*Over->path.tlen + 1000;
                    Trace64  = realloc(Trace64,sizeof(int64)*TraceMax);
                  }
                Write_Aln_Overlap(Ofile,Over);
                Compress_TraceTo8(Over,0);
                Write_Aln_Trace(Ofile,Over->path.trace,Over->path.tlen,Trace64);
                Decompress_TraceTo16(Over);

#ifdef SHOW_ALIGNMENTS
                printf(" Hit spans %d-%d\n",bpath->abpos,bpath->bepos);
                Compute_Trace_PTS(Align,Work,100,GREEDIEST,len-wide,len+wide);
                Print_Alignment(stdout,Align,Work,8,100,10,0,10);
#endif

                // cut_trace(Align);

                beg = bpath->abpos - off;
                end = bpath->bepos - off;
                for (f = n+1; f < nseed; f++)
                  { s = KMER(seed[f].start);
                    for (g = count[s]; g < count[s+1]; g++)
                      if (beg <= index[g] && index[g] <= end) 
                        index[g] = -1;
                  }

                last = end;

                sat[nsat].beg = beg;
                sat[nsat].end = end;
                sat[nsat].bc  = seq[beg];
                sat[nsat].ec  = seq[end];
                seq[beg] = 4;
                seq[end] = 4;
                nsat += 1;

                if (bpath->bepos > outhit)
                  outhit = bpath->bepos;
              }
          }
      }

    for (n = 0; n < nsat; n++)
      { seq[sat[n].beg] = sat[n].bc;
        seq[sat[n].end] = sat[n].ec;
      }

    free(sat);
    free(seed);
    return (outhit);
  }
}


/*******************************************************************************************
 *
 *  SPECTRUM CHAIN DETECTOR
 *
 ********************************************************************************************/

typedef struct
  { int d;
    int c;
  } Chord;

static int CSORT(const void *l, const void *r)
{ Chord *x = (Chord *) l;
  Chord *y = (Chord *) r;

  return (y->c - x->c);
}

static int chord_block(uint8 *seq, int off, int len, void *block)
{ double  freq[4];
  int16  *count;  // 64K + 1
  int16  *map;    // 64K
  int16  *index;  // 32K
  Dnode  *node;   // 3.2K + 1
  int16  *edge;   // 3.2K * 4
  int16  *ewgt;   // 3.2K * 4
  int16  *path;   // 3.2K * 4
  int16  *pwgt;   // 3.2K * 4
  int16  *low;    // 3.2K
  int16  *num;    // 3.2K
  int16  *list;   // 3.2K
  int16  *bflag;  // 3.2K
  int16  *blist;  // 3.2K
  Elem   *link;   // 3.2K * 4

  Hit *seed;
  int  nseed, maxseed;

  count = (int16 *) block;         //  64K * 5 + 3.2K * (sizeof(Dnode) + 16)
  index = count + 0x10004;
  node  = (Dnode *) (index + 0x8000);
  edge  = (int16 *) (node + (0x8000/10) + 1);
  ewgt  = (int16 *) (edge + (0x20000/10));
  map   = (int16 *) (ewgt + (0x20000/10));

  path  = map;
  pwgt  = path + (0x20000/10);
  num   = pwgt + (0x20000/10);
  low   = num  + (0x8000/10);
  list  = low  + (0x8000/10);
  bflag = list  + (0x8000/10);
  blist = bflag  + (0x8000/10);
  link  = (Elem *) (blist + (0x8000/10));

  printf("\nPANEL %d-%d\n",off,off+0x8000);

  //  Build index of 8-mers in count and index

  { uint8 *s7;
    int    i, p, x, l7, fq[4];
    uint16 kmer;

    s7 = seq+7;
    l7 = len-7;

    for (i = 0; i < 0x10000; i++)
      count[i] = 0;

    for (i = 0; i < 4; i++)
      fq[i] = 0;

    kmer = seq[0];                //  Count # of each 8-mer
    freq[seq[0]] += 1;
    for (i = 1; i < 7; i++)
      { kmer = (kmer << 2) | seq[i];
        fq[seq[i]] += 1;
      }
    for (i = 0; i < l7; i++)
      { x = s7[i];
        kmer = (kmer << 2) | x;
        count[kmer] += 1;
        fq[x] += 1;
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

    for (i = 0x10000; i > 0; i--)
      count[i] = count[i-1];
    count[0] = 0;

    for (i = 0; i < 4; i++)
      freq[i] = (1.*fq[i])/0x8000;
  }

  { int i;
    int b, e;
    int k, h;
    int d;

    int   ncnt;
    Chord hist[500];

    for (d = 0; d < 500; d++)
      { hist[d].d = d;
        hist[d].c = 0;
      }

    e = 0;
    for (i = 0; i < 0x10000; i++)
      { b = e;
        e = count[i+1];
        h = index[b++];
        while (b < e)
          { k = index[b++];
            d = k-h;
            h = k;
            if (d < 500)
              hist[d].c += 1;
          }
      }

    ncnt = 0;
    for (d = 1; d < 500; d++)
      if (hist[d].c > (d>>2))
        hist[ncnt++] = hist[d];

    qsort(hist,ncnt,sizeof(Chord),CSORT);

    printf("Histo: %d\n",ncnt);
    for (d = 0; d < ncnt; d++)
      printf(" %4d: %5d\n",hist[d].d,hist[d].c);
  }

  return (0);
}
            


/*******************************************************************************************
 *
 *  SEED CHAIN DETECTOR
 *
 ********************************************************************************************/

#define DIAG_MAX 8000
#define BAND_MAX  375  //  = ((DIAG_MAX-1)>>3)+1

#define LAYERS     3   //  Must be >= 1

#undef  SHOW_SEEDS
#undef  SHOW_CHAINS

static void Print_Hit(uint8 *seq, int beg, int diag, int len)
{ static uint16 mask[9] = { 0, 0x3, 0xf, 0x3f, 0xff, 0x3ff, 0xfff, 0x3fff, 0xffff };
  uint16 v, w;
  int    i, p;

  printf("    %5d-%5d: (",beg,beg+len);
  v = 0;
  for (i = 0; i < diag; i++)
    v = v << 2 | seq[beg+i];
  w = v;
  p = 0;
  for (i = 0; i < diag; i++)
    { v = (v << 2 | seq[beg+i]) & mask[diag];
      if (v < w)
        { w = v;
          p = i+1;
        }
    }
  for (i = 0; i < diag; i++)
    printf("%c",dna[seq[beg+i+p]]);
  printf(")^%d",len/diag);
  if (len%diag > 0)
    printf(".%d",len%diag);
  if (p != 0)
    printf("^%d",diag-p);
  printf("\n");
}

typedef struct
  { uint16  diag;
    uint16  ibeg;
  } Seed;

static int spectrum_block(uint8 *seq, int off, int len, void *block)
{ int    i, p, x, c;
  int    d, e, f;
  uint16 kmer;
  int    *count;  // 0x10000
  int    *diags;  // DIAG_MAX
  uint16 *index;  // 0x08000
  Seed   *post;   // 0x08000
  Seed   *hits;   // 0x08000
  uint8 *s7;
  int    l7;

  index = (uint16 *) block;
  post  = (Seed *) (index +  0x8000);
  hits  = post + 0x8000;
  count = (int *) (hits + 0x8000);
  diags = count +  0x10000; 

  (void) off;

  printf("\nPANEL %d-%d\n",off,off+0x8000);
  fflush(stdout);

  for (i = 0; i < 0x10000; i++)   //  Init counters
    count[i] = 0;

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
  index[l7] = 0;

  index[0] |= 0x8000;           //  mark bucket ends and reset count
  for (i = 0; i < 0xffff; i++)
    { index[count[i]] |= 0x8000;
      count[i] = 0;
    }
  count[0xffff] = 0;

#ifdef SORT1
  for (i = 0; i < l7; i++)
    { f = index[i];
      p = index[i] & 0x7fff;
      printf("%c %5d: ",p==f?' ':'+',p);
      Print_Seq(seq+p,8);
      printf("\n");
    }
#endif

  e = index[0] & 0x7fff;
  for (i = 1; i < l7; i++)     //  count ibeg's of all same-kmer position pairs within LAYERS
    { f = index[i];            //    of each other in the index and with diag < DIAG_MAX (3Kbp)
      if (f < 0x8000)
        { d = f-e;
          if (d < DIAG_MAX)
            count[e] += 1;
          e = f;
        }
      else
        e = f & 0x7fff;
    }

  p = 0;                          //  turn counts into ptrs
  for (i = 0; i < 0x08000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  for (i = 0; i < DIAG_MAX; i++)        //  init diagonal tube counters
    diags[i] = 0;

  e = index[0] & 0x7fff;
  for (i = 1; i < l7; i++)       //   place seed pairs in post sorted on ibeg using count
    { f = index[i];              //     ptrs.  Also count diagonal tubes for next sort.
      if (f < 0x8000)
        { d = f-e;
          if (d < DIAG_MAX)
            { c = count[e]++;
              post[c].ibeg = e;
              post[c].diag = d;
              diags[d] += 1;
            }
          e = f;
        }
      else
        e = f & 0x7fff;
    }

  p = 0;                        //  turn diag counts into sort ptrs
  for (i = 0; i < DIAG_MAX; i++)
    { x = diags[i];
      diags[i] = p;
      p += x;
    }

#ifdef SORT2
  printf("Sorted on Anti\n");
  for (c = 0; c < p; c++)
    printf(" %5d %5d\n",post[c].diag,post[c].ibeg); 
#endif

  for (i = 0; i < p; i++)       //  place ibeg/diag pairs in hits in order of diag then ibeg
    { c = post[i].diag;
      hits[diags[c]++] = post[i];
    }

#ifdef SHOW_SEEDS
  p = 0;
  for (i = 1; i < DIAG_MAX; i++)
    { f = diags[i];
      if (p >= e)
        continue;
      printf("Diagonal %d : %d\n",i,f);
      for ( ; p < f; p++)
        { d = hits[p].diag;
          e = hits[p].ibeg;
          printf("   %4d : %5d  ",d,e);
          Print_Seq(seq+e,8);
          printf("\n");
        }
    }
#endif

  { int   ncnt;
    int   outhit;
    Chord hist[DIAG_MAX];
    int   wide, anti, last;
    Path *bpath;

    ncnt = 0;
    p = diags[1];
    for (i = 2; i < DIAG_MAX; i++)
      { f = diags[i];
        if (f-p > 1 && f-p > (i>>6))
          { hist[ncnt].c = f-p;
            hist[ncnt].d = i;
            ncnt += 1;
          }
        p = f;
      }

    qsort(hist,ncnt,sizeof(Chord),CSORT);

    outhit = 0;
    printf("Histo: %d\n",ncnt);
    for (i = 0; i < ncnt; i++)
      { d = hist[i].d;
        last = -1;
// printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
        for (x = diags[d-1]+1; x < diags[d]; x++)
          { p = hits[x].ibeg;
// printf("  p = %d\n",p);
            if (p < last || p - hits[x-1].ibeg > d || seq[p] >= 4)
              continue;
            wide = .2*d;
            if (wide < 1)
              wide = 1;
            anti = 2*(off + p) + d;
            bpath = Local_Alignment(Align,Work,Spec,d,d,anti,wide,wide);
// if (bpath == NULL)
  // printf("    NULL\n");
// else
  // printf("    %d (%d)\n",bpath->bepos-bpath->abpos,2*d);
            if (bpath->bepos - off > last)
              last = bpath->bepos - off;
            if (bpath != NULL && bpath->bepos-bpath->abpos >= 1.95*d)
              { if (Over->path.tlen > TraceMax)
                  { TraceMax = 1.2*Over->path.tlen + 1000;
                    Trace64  = realloc(Trace64,sizeof(int64)*TraceMax);
                  }
                Write_Aln_Overlap(Ofile,Over);
                Compress_TraceTo8(Over,0);
                Write_Aln_Trace(Ofile,Over->path.trace,Over->path.tlen,Trace64);
                Decompress_TraceTo16(Over);

#ifdef SHOW_ALIGNMENTS
                if (last < 0)
                  printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
                printf(" Hit spans %d-%d\n",bpath->abpos,bpath->bepos);
                Compute_Trace_PTS(Align,Work,100,GREEDIEST,d-wide,d+wide);
                Print_Alignment(stdout,Align,Work,8,100,10,0,10);
#endif
                for (f = bpath->abpos - off; f < last; f++)
                  seq[f] = 4;

                if (bpath->bepos > outhit)
                  outhit = bpath->bepos;
              }
          }
      }

    return (outhit);
  }
}

static void filter_block(uint8 *seq, int off, int len, void *block)
{ int    i, p, x, c;
  int    d, e, f;
  uint16 kmer;
  uint32 u;
  int    *count;  // 0x10000
  int    *diags;  // BAND_MAX
  uint16 *index;  // 0x08000
  Seed   *post;   // 0x08000 * LAYERS
  Seed   *hits;   // 0x08000 * LAYERS
  uint8 *s7;
  int    l7;

  index = (uint16 *) block;
  post  = (Seed *) (index +  0x8000);
  hits  = post + 0x8000 * LAYERS;
  count = (int *) (hits + 0x8000 * LAYERS);
  diags = count +  0x10000; 

  (void) off;

#if defined(SHOW_SEED) || defined(SHOW_CHAINS)
  printf("\nPANEL %d-%d\n",off,off+0x8000);
  fflush(stdout);
#endif

  for (i = 0; i < 0x10000; i++)   //  Init counters
    count[i] = 0;

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
  index[l7] = 0;

  index[0] |= 0x8000;           //  mark bucket ends and reset count
  for (i = 0; i < 0xffff; i++)
    { index[count[i]] |= 0x8000;
      count[i] = 0;
    }
  count[0xffff] = 0;

#ifdef SORT1
  for (i = 0; i < l7; i++)
    { f = index[i];
      p = index[i] & 0x7fff;
      printf("%c %5d: ",p==f?' ':'+',p);
      Print_Seq(seq+p,8);
      printf("\n");
    }
#endif

  for (i = 1; i < l7; i++)     //  count ibeg's of all same-kmer position pairs within LAYERS
    { f = index[i];            //    of each other in the index and with diag < DIAG_MAX (3Kbp)
      if (f < 0x8000)
        for (x = 1; x <= LAYERS; x++)
          { p = index[i-x];
            e = p & 0x7fff;
            d = f-e;
            if (d < DIAG_MAX)
              count[e] += 1;
            if (p >= 0x8000)
              break;
          }
    }

  p = 0;                          //  turn counts into ptrs
  for (i = 0; i < 0x08000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  for (i = 0; i < BAND_MAX; i++)        //  init diagonal tube counters
    diags[i] = 0;

  for (i = 1; i < l7; i++)       //   place seed pairs in post sorted on ibeg using count
    { f = index[i];              //     ptrs.  Also count diagonal tubes for next sort.
      if (f < 0x8000)
        for (x = 1; x <= LAYERS; x++)
          { p = index[i-x];
            e = p & 0x7fff;
            d = f-e;
            if (d < DIAG_MAX)
              { c = count[e]++;
                post[c].ibeg = e;
                post[c].diag = d;
                diags[d>>3] += 1;
              }
            if (p >= 0x8000)
              break;
          }
    }

  p = 0;                        //  turn diag counts into sort ptrs
  for (i = 0; i < BAND_MAX; i++)
    { x = diags[i];
      diags[i] = p;
      p += x;
    }

#ifdef SORT2
  printf("Sorted on Anti\n");
  for (c = 0; c < p; c++)
    printf(" %5d %5d\n",post[c].diag,post[c].ibeg); 
#endif

  for (i = 0; i < p; i++)       //  place ibeg/diag pairs in hits in order of diag then ibeg
    { u = post[i].diag;
      hits[diags[u>>3]++] = post[i];
    }

#ifdef SHOW_SEEDS
  p = 0;
  for (i = 0; i < BAND_MAX; i++)
    { e = diags[i];
      if (p >= e)
        continue;
      printf("Band %d\n",i);
      for ( ; p < e; p++)
        { d = hits[p].diag;
          e = hits[p].ibeg;
          printf("   %4d : %5d  ",d,e);
          Print_Seq(seq+e,8);
          printf("\n");
        }
    }
#endif

#ifdef SHOW_CHAINS

  { int p1, p2;
    int t1, t2;
    int a1, a2;
    int b1, b2;
    int lowd;
    int ep;
    int a, w;
    int d, x;
    int nel, hlen;
    int cbeg, cend, cov;
    int bucket[16];

    lowd = -8;
    t1 = 0;
    t2 = diags[0];
    for (i = 1; i < BAND_MAX; i++)
      { p1 = t1;
        t1 = p2 = t2;
        t2 = diags[i];
        a1 = hits[p1].ibeg;
        a2 = hits[p2].ibeg;
        lowd += 8;
        cend = -129;
        ep   = 0;
        nel  = (t1-p1)+(t2-p2);
        cov  = 0;
        while (1)
          { w = (p2 >= t2 || a1 <= a2);
            if (w)
              a = a1;
            else
              a = a2;

            if (nel <= 0 || a >= cend+128 || a >= ep)
              { if (cend > 0 && cov > .5*(cend - cbeg) && cend-cbeg > 6*i && b1 < p1)
                  { printf("  Hit @%5d len = %4d cov = %4d in bands %d - %d\n",
                           off + cbeg,cend-cbeg,cov,lowd,8*(i+1));
                    if ((p1+p2)-(b1+b2) == 1)
                      { if (hits[b1].diag <= 8)
                          { d = hits[b1].diag;
                            x = hits[b1].ibeg;
                            Print_Hit(seq,x,d,8+d);
                          }
                      }
                    else
                      { for (d = 0; d < 16; d++)
                          bucket[d] = 0;
                        hlen = 0;
                        while (b1+b2 < p1+p2)
                          { if (b2 >= p2 || hits[b1].ibeg <= hits[b2].ibeg)
                              { bucket[hits[b1].diag-lowd] += 1;
                                post[hlen++] = hits[b1++];
                              }
                            else
                              { bucket[hits[b2].diag-lowd] += 1;
                                post[hlen++] = hits[b2++];
                              }
                          }

                        if (lowd == 0)
                          { int n, len, beg, fst;

                            while (1)
                              { for (d = 0; d < 8; d++)
                                  if (bucket[d] > 0)
                                    break;
                                if (d == 8) break;

                                len = 0;
                                fst = -1;
                                beg = -1;
                                for (n = 0; n <= hlen; n++)
                                  if (post[n].diag == d || n == hlen)
                                    { if (n < hlen && post[n].ibeg == beg+len)
                                        len += 1;
                                      else
                                        { if (len > 0)
                                            { len += d+7;
                                              Print_Hit(seq,beg,d,len);

                                              for (x = fst; x < n; x++)
                                                if (post[x].ibeg+d+8 <= beg+len)
                                                  { bucket[post[x].diag] -= 1;
                                                    post[x].diag += DIAG_MAX;
                                                  }
                                            }
                                          beg = post[n].ibeg;
                                          fst = n;
                                          len = 1;
                                        }
                                    }
                                }

                          }
                        for (d = 0; d < 16; d++)
                          if (bucket[d] > 0)
                            printf(" %4d: %4d\n",lowd+d,bucket[d]);

/*
  kmer = seq[0];
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    index[i] = kmer = (kmer << 2) | s7[i];

                        for (int n = 0; n < hlen; n++)
                          { d = post[n].diag;
                            x = post[n].ibeg;
                            if (d > DIAG_MAX)
                              { d -= DIAG_MAX;
                                printf("     %4d : %5d  ",d,x);
                                post[n].diag = d;
                              }
                            else
                              printf("   + %4d : %5d  ",d,x);
                            Print_Seq(seq+x,8);
                            printf(" %d-%d  %d-%d",x,x+8,x+d,x+d+8);
                            printf("\n");
                          }
*/
                      }
                  }
                if (nel <= 0)
                  break;
                cbeg = a;
                cov  = 8;
                b1   = p1;
                b2   = p2;
              }
            else
              { if (a < cend)
                  cov += (a+8) - cend;
                else
                  cov += 8;
              }
            cend = a + 8;

            if (w)
              { ep = a1+hits[p1].diag+8; 
                a1 = hits[++p1].ibeg;
              }
            else
              { ep = a2+hits[p2].diag+8; 
                a2 = hits[++p2].ibeg;
              }
            nel -= 1;
          }
      }
  }

#endif

if (p == 9924)
  { Print_Seq(seq+14220,173);
    printf("\n");
  }

  printf("HITS = %d\n",p); fflush(stdout);
}


int main(int argc, char *argv[])
{ char    *spath, *tpath;
  int      ftype;
  GDB     _gdb, *gdb = &_gdb;
  char    *cpath, *APATH, *AROOT;

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

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,0);

  if (ftype != IS_GDB)
    Create_GDB(gdb,spath,ftype,1,NULL);
  else
    { Read_GDB(gdb,tpath);
      if (gdb->seqs == NULL)
        { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
          exit (1);
        }
    }
  free(tpath);

  APATH = PathTo(argv[2]);
  AROOT = Root(argv[2],".1aln");
  cpath = getcwd(NULL,0);

  Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),1,Prog_Name,VERSION,Command_Line,
                         TSPACE,spath,spath,cpath);
  free(cpath);
  free(AROOT);
  free(APATH);
  free(spath);

  if (VERBOSE)
    { fprintf(stderr," Database built, begin scan\n");
      fflush(stderr);
    }

  { uint8 *buffer;
    void  *block;
    int    clen, last;
    int    i, p, d;
    Alignment _align;
    Overlap   _ovl;

    Work = New_Work_Data();
    Spec = New_Align_Spec(.7,TSPACE,gdb->freq,0);
    Align = &_align;
    Over  = &_ovl;
    Align->path = &(Over->path);

    block  = malloc(0x40000*LAYERS + 0x50000 + 4*BAND_MAX);   // 256KB * LAYERS + 332KB
    buffer = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;

    if (block == NULL || buffer == NULL)
      { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
        exit (1);
      }

    TraceMax = 10000;
    Trace64  = malloc(sizeof(int64)*TraceMax);
    Align->flags = 0;
    Over->flags  = 0;
    for (i = 0; i < gdb->ncontig; i++)
      { Get_Contig(gdb,i,NUMERIC,(char *) buffer);
        printf("CONTIG %d\n",i+1);
        clen = gdb->contigs[i].clen;
        Align->aseq = Align->bseq = (char *) buffer;
        Align->alen = Align->blen = clen;
        Over->aread = i;
        Over->bread = i;
        last = -1;
        if (clen < 0x8000)
          { // chord_block(buffer,0,clen,block);
            // euler_block(buffer,0,clen,block);
            spectrum_block(buffer,0,clen,block);
            // filter_block(buffer,0,clen,block);
          }
        else
          for (p = 0; p+0x2000 <= clen; p += 0x6000)
            { d = buffer[p-1];
              buffer[p-1] = 4;
              if (p+0x8000 > clen)
                { // chord_block(buffer+p,p,clen-p,block);
                  // euler_block(buffer+p,p,clen-p,block);
                  spectrum_block(buffer+p,p,clen-p,block);
                  // filter_block(buffer+p,p,clen-p,block);
                  buffer[p-1] = d;
                }
              else
                { // last = chord_block(buffer+p,p,0x8000,block);
                  // last = euler_block(buffer+p,p,0x8000,block);
                  last = spectrum_block(buffer+p,p,0x8000,block);
                  // filter_block(buffer+p,p,0x8000,block);
                  buffer[p-1] = d;
                  if (last >= p+0x6000)
                    p = last-0x6000;
                }
            }
      }

    free(Trace64);
    free(buffer-1);
    free(block);

    Free_Align_Spec(Spec);
    Free_Work_Data(Work);
  }

  oneFileClose(Ofile);

  Close_GDB(gdb);

  if (VERBOSE)
    { TimeTo(stderr,0,1);
      TimeTo(stderr,1,0);
    }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}

/*
 * Written by Joyce Whang
 * Overlapping Non-exhaustive Graph Clustering
 * Jan. 9, 2014
 *
 */

#include <metis.h>
#include <iostream>
#include <vector>
#include <list>
#include "proto_ovfuncs_mllk.h"
using namespace std;

/*************************************************************************
* This function is the entry point for MLKKM
**************************************************************************/
void MLKKM_PartGraphKway_overlap(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                         idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *chainlength, 
                         int *options, int *edgecut, idxtype *part, int levels, std::vector< std::list<int> >& ClustResults,
			 float alpha, float beta, float sigma)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "MLKKM: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  MLKKM_WPartGraphKway_overlap(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, chainlength,
                       tpwgts, options, edgecut, part, levels, ClustResults, alpha, beta, sigma);

  free(tpwgts);
}


/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void MLKKM_WPartGraphKway_overlap(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                          idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *chainlength,
                          float *tpwgts, int *options, int *edgecut, idxtype *part, int levels, 
			  std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma)
{
  int i, j;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);
  
  ctrl.CType = KMETIS_CTYPE; 
  ctrl.IType = KMETIS_ITYPE;
  ctrl.RType = KMETIS_RTYPE;
  //ctrl.dbglvl = KMETIS_DBGLVL;
  ctrl.dbglvl = 7; // modified by Joyce

  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = levels;
  ctrl.maxvwgt = floor(1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt) : (*nvtxs))/ctrl.CoarsenTo));
  ctrl.maxvwgt *= 100;
  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MLKKMPartitioning_overlap(&ctrl, &graph, *nparts, *chainlength, part, tpwgts, 1.03, ClustResults, alpha, beta, sigma);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  //IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}



/*************************************************************************
* This function takes a graph and produces a k-way partitioning of it
**************************************************************************/
int MLKKMPartitioning_overlap(CtrlType *ctrl, GraphType *graph, int nparts, int chain_length, idxtype *part, float *tpwgts, float ubfactor, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma)
{
  int i, j, nvtxs, tvwgt, tpwgts2[2];
  GraphType *cgraph;
  int wgtflag=3, numflag=0, options[10], edgecut;
  float ncut;
  idxtype *cptr, *cind;
  int numcomponents;
  char *mlwkkm_fname = "coarse.graph";

  cptr = idxmalloc(graph->nvtxs, "MLKKMPartitioning: cptr");
  cind = idxmalloc(graph->nvtxs, "MLKKMPartitioning: cind");

  cgraph = Coarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  AllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1; 
  options[OPTION_CTYPE] = MATCH_SHEMKWAY;
  options[OPTION_ITYPE] = IPART_GGPKL;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;
  
  METIS_WPartGraphRecursive(&cgraph->nvtxs, cgraph->xadj, cgraph->adjncy, cgraph->vwgt, 
                           cgraph->adjwgt, &wgtflag, &numflag, &nparts, tpwgts, options, 
                           &edgecut, cgraph->where);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, printf("Initial %d-way partitioning cut: %d\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));
  
  MLKKMRefine_overlap(ctrl, graph, cgraph, nparts, chain_length, tpwgts, ubfactor, ClustResults, alpha, beta, sigma);
  
  //idxcopy(graph->nvtxs, graph->where, part);
  
  GKfree((void **) &graph->gdata, (void **) &graph->rdata, LTERM);

  //return graph->mincut;
  return -1;

}


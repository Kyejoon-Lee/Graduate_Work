/*
 * Written by Joyce Whang
 * NEO-K-Means Graph Clustering
 * date: June 23, 2015
 *
 */

#include <metis.h>
#include <iostream> 
#include <vector>
#include <list>
#include <float.h>
#include <queue>
#include <assert.h>
#include "proto_ovfuncs_mllk.h"
#include "io.h"
#include "mex.h"

void One_level_overlap_refine(GraphType *graph, int nparts, idxtype *w, std::vector< std::list<int> >& AssVec, float alpha, float beta, float sigma);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
  mwSize *jc, *ir, *jc_init, *ir_init;
  double *pr, *pr_init;
  int edgenum, edgenum_init;
  GraphType graph, graph_init;
  idxtype *xadj, *adjncy, *adjwgt, *xadj_init, *adjncy_init, *adjwgt_init, *part;
  float alpha, beta;
  float sigma = 0;
  int nparts, one_level, wgtflag, levels;
  int options[11];
  int numflag = 0, edgecut, chain_length=0;
  timer TOTALTmr, NEOTmr, IOTmr;
  //char filename[256];

  ir = mxGetIr(prhs[0]);
  jc = mxGetJc(prhs[0]);
  pr = mxGetPr(prhs[0]);
  edgenum = (int) (mxGetScalar(prhs[1]));
  nparts = (int) (mxGetScalar(prhs[2]));
  alpha = (float) (mxGetScalar(prhs[3]));
  beta = (float) (mxGetScalar(prhs[4]));
  one_level = (int) (mxGetScalar(prhs[5]));

  cleartimer(TOTALTmr);
  cleartimer(NEOTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);

  /*---------------------------------------- IO process ----------------------------------------*/
  starttimer(IOTmr);  
  InitGraph(&graph);
  graph.nvtxs = mxGetM(prhs[0]);
  vector< list<int> > ClustResults(graph.nvtxs); // cluster assignment vector
  graph.nedges = edgenum;
  //cout<<"(mexFunction) nvtxs: "<<graph.nvtxs<<", nedges: "<<graph.nedges<<"\n";

  xadj = graph.xadj = idxsmalloc(graph.nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph.adjncy = idxmalloc(graph.nedges, "ReadGraph: adjncy");
  adjwgt = graph.adjwgt = idxsmalloc(graph.nedges, 1, "RefineKWay: adjwgt");
  wgtflag = 1;

  for(int i = 0; i < edgenum; i++){
     adjncy[i] = (idxtype) ir[i];
     adjwgt[i] = (idxtype) pr[i];
  }
  for(int i = 0; i <= graph.nvtxs; i++){
     xadj[i] = (idxtype) jc[i];
     //cout<<"xadj["<<i<<"]: "<<xadj[i]<<"\n";
  }

  if (nparts < 2) {
    printf("The number of clusters should be greater than 1.\n");
    exit(0);
  }

  levels = amax((graph.nvtxs)/(40*log2_metis(nparts)), 20*(nparts));

  printf("\n----------------------------------------------------------------------\n");
  printf("NEO-K-Means Graph Clustering, Copyright 2014, Joyce Whang \n");
  printf("Graph Information:\n");
  printf("#Vertices: %d, #Edges: %d, ", graph.nvtxs, graph.nedges/2);

  part = idxmalloc(graph.nvtxs, "main: part");
  printf("#Clusters: %d\n", nparts);
  options[0] = 0;

  /*
  cout<<"graph.nvtxs: "<<graph.nvtxs<<"\n";
  for (int i = 0; i<graph.nvtxs; i++){
     for (int j=xadj[i]; j<xadj[i+1]; j++){
        int neighbor = adjncy[j];
	cout<<"* i: "<<i<<", neighbor: "<<neighbor<<", adjwgt: "<<adjwgt[neighbor]<<"\n";
     }
     cout<<"\n";
  }
  */

  if(one_level){
     std::cout<<"The initialization is provided... Only one level refinement will be conducted...\n";
     // read initialization
     ir_init = mxGetIr(prhs[6]);
     jc_init = mxGetJc(prhs[6]); 
     edgenum_init = (int) (mxGetScalar(prhs[8]));
     InitGraph(&graph_init);
     graph_init.nvtxs = (int) (mxGetScalar(prhs[7]));
     graph_init.nedges = edgenum_init;
     //cout<<"(mexFunction) graph_init nvtxs: "<<graph_init.nvtxs<<", nedges: "<<graph_init.nedges<<"\n";
     xadj_init = graph_init.xadj = idxsmalloc(graph_init.nvtxs+1, 0, "ReadGraph: xadj_init");
     adjncy_init = graph_init.adjncy = idxmalloc(graph_init.nedges, "ReadGraph: adjncy_init");
     adjwgt_init = graph_init.adjwgt = idxsmalloc(graph_init.nedges, 1, "RefineKWay: adjwgt_init");
     for(int i = 0; i < edgenum_init; i++){
        adjncy_init[i] = (idxtype) ir_init[i];
        adjwgt_init[i] = 1; // input is a binary assignment matrix
	//cout<<"adjncy_init["<<i<<"]: "<<adjncy_init[i]<<", adjwgt_init["<<i<<"]: "<<adjwgt_init[i]<<" \n";
     }
     for(int i = 0; i <= graph_init.nvtxs; i++){
        xadj_init[i] = (idxtype) jc_init[i];
	//cout<<"xadj_init["<<i<<"]: "<<xadj_init[i]<<"\n";
     }
     
     //cout<<"graph_init.nvtxs: "<<graph_init.nvtxs<<"\n";
     for (int i = 0; i<graph_init.nvtxs; i++){
        for (int j=xadj_init[i]; j<xadj_init[i+1]; j++){
           int neighbor = adjncy_init[j];
	   ClustResults[i].push_back(neighbor);
	   //cout<<"* i: "<<i<<", neighbor: "<<neighbor<<", adjwgt_init: "<<adjwgt_init[neighbor]<<"\n";
        }
	//cout<<"\n";
     }
     /*
     // print
     cout << "convert to vector-list form...\n";
     for(int iter = 0; iter < graph_init.nvtxs; iter++){
        cout << "node no.: " << iter << " ... ";
        for(list<int>::iterator i=ClustResults[iter].begin(); i!=ClustResults[iter].end(); i++){
	   cout << *i << " ";
        }
     cout << "\n";
     }
     */

  }

  stoptimer(IOTmr);

  /*---------------------------------------- clustering ----------------------------------------*/
  starttimer(NEOTmr);
  /* Main Module */
  if(one_level){ 
     // one-level refinement
     idxtype *w = idxsmalloc(graph.nvtxs, 0, "pingpong: weight");
     Compute_Weights(&graph, w);
     One_level_overlap_refine(&graph, nparts, w, ClustResults, alpha, beta, sigma);
  }else{ 
     // multilevel refinement
     MLKKM_PartGraphKway_overlap(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
	              &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, part, levels, ClustResults, alpha, beta, sigma);
  }
  stoptimer(NEOTmr);
 
  /*---------------------------------------- IO process ----------------------------------------*/
  starttimer(IOTmr);
  /*
  cout << "\n ----- print ClustResults... neo-k-means \n";
  for(vector< list<int> >::iterator iter = ClustResults.begin(); iter != ClustResults.end(); iter++){
    for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
     cout << *i << " ";    
    }
    cout<<"\n";
  }
  */
  //filename={'t','e','s','t'};
  //WriteOverlapClust(filename, ClustResults, graph.nvtxs, nparts, alpha, beta, sigma); 
  double *out, *out2;
  int totalAssign = graph.nvtxs + round(graph.nvtxs * alpha);
  cout<<"** totalAssign: "<<totalAssign<<"\n";
  plhs[0] = mxCreateDoubleMatrix(totalAssign,1,mxREAL);
  out = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(totalAssign,1,mxREAL);
  out2 = mxGetPr(plhs[1]);
  int nodeID = 0, ind = 0;
  for(vector< list<int> >::iterator iter = ClustResults.begin(); iter != ClustResults.end(); iter++){
    nodeID++;
    for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
     //cout <<"nodeID: "<<nodeID<<", clustID: "<< *i+1 << " \n";
     out[ind] = (double) nodeID;   
     out2[ind] = (double) (*i+1); 
     ind = ind+1;
    }
  }
  stoptimer(IOTmr);
  stoptimer(TOTALTmr);
  

  printf("\nTiming Information:\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Clustering:   \t\t %7.3f   (NEO time)\n", gettimer(NEOTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("----------------------------------------------------------------------\n");


  //GKfree((void **) &graph.xadj, (void **) &graph.adjncy, (void **) &graph.adjwgt, (void **) &part);

}  

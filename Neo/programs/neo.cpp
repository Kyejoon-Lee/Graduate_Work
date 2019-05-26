/*
 * Written by Joyce Whang
 * NEO-K-Means Graph Clustering
 * Written: Jan. 9, 2014
 * Updated: June 23, 2015
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

void One_level_overlap_refine(GraphType *graph, int nparts, idxtype *w, std::vector< std::list<int> >& AssVec, float alpha, float beta, float sigma);

main(int argc, char *argv[])
{
  int i, nparts=1, options[11];
  idxtype *part;
  GraphType graph;
  char filename[256];
  int numflag = 0, wgtflag = 0, edgecut, chain_length=0;
  timer TOTALTmr, NEOTmr, IOTmr;
  int levels = 0;
  float alpha = 0.2;
  float beta = 0;
  float sigma = 0;
  int one_level = 0;
  char init_filename[256];

  if(argc < 3){
    printf("Help ----- \n Please provide at least two inputs:\n neo graph_file number_of_clusters\n");
    exit(1);
  }

  for (argv++; *argv != NULL; argv++){
    if ((*argv)[0] == '-'){
      switch ((*argv)[1])
      {
	case 'a':
          alpha = atof(*(++argv));
	  break;
	case 'b':
          beta = atof(*(++argv));
	  break;
	case 's':
          sigma = atof(*(++argv));
	  break;
	case 'i':
	  strcpy(init_filename, *(++argv));
	  one_level = 1;
	  break;
	default:
	  printf("Invalid option %s\n", *argv);
	  exit(0);
      }
    }else{
      strcpy(filename, *argv);
      nparts = atoi(*(++argv));
    }
  }
 
  if (nparts < 2) {
    printf("The number of partitions should be greater than 1.\n");
    exit(0);
  }

  cleartimer(TOTALTmr);
  cleartimer(NEOTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);
  ReadGraph(&graph, filename, &wgtflag);
  stoptimer(IOTmr);

  levels = amax((graph.nvtxs)/(40*log2_metis(nparts)), 20*(nparts));

  printf("\n----------------------------------------------------------------------\n");
  printf("NEO-K-Means Graph Clustering, Copyright 2014, Joyce Whang \n");
  printf("Graph Information:\n");
  printf("  Name: %s, \n  #Vertices: %d, #Edges: %d, ", filename, graph.nvtxs, graph.nedges/2);

  part = idxmalloc(graph.nvtxs, "main: part");
  printf("#Clusters: %d\n", nparts);
  options[0] = 0;

  // cluster assignment vector
  vector< list<int> > ClustResults(graph.nvtxs);

  starttimer(NEOTmr);
  /* Main Module */
  if(one_level==0){ // multilevel refinement
     MLKKM_PartGraphKway_overlap(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
	              &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, part, levels, ClustResults, alpha, beta, sigma);
  }else{ // one level refinement
     std::cout<<"The initialization is provided... Only one level refinement will be conducted...\n";
     // read initialization
     ReadInitialization(init_filename, graph.nvtxs, ClustResults);
     // perform one level refinement
     idxtype *w = idxsmalloc(graph.nvtxs, 0, "pingpong: weight");
     Compute_Weights(&graph, w);
     One_level_overlap_refine(&graph, nparts, w, ClustResults, alpha, beta, sigma);
  }
  stoptimer(NEOTmr);
 
  /*
  cout << "\n ----- print ClustResults... neo-k-means \n";
  for(vector< list<int> >::iterator iter = ClustResults.begin(); iter != ClustResults.end(); iter++){
    for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
     cout << *i << " ";    
    }
    cout<<"\n";
  }
  */

  starttimer(IOTmr);
  WriteOverlapClust(filename, ClustResults, graph.nvtxs, nparts, alpha, beta, sigma); 
  stoptimer(IOTmr);
  stoptimer(TOTALTmr);
  

  printf("\nTiming Information:\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Clustering:   \t\t %7.3f   (NEO time)\n", gettimer(NEOTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("----------------------------------------------------------------------\n");
  FILE *fp;

  GKfree((void **) &graph.xadj, (void **) &graph.adjncy, (void **) &graph.vwgt, (void **) &graph.adjwgt, (void **) &part, LTERM);
}  



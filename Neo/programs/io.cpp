/*
 * Modified by Joyce Whang
 * Overlapping Non-exhaustive Graph Clustering
 * Jan. 9, 2014
 *
 */

#include <metis.h>
#include <iostream> 
#include <vector>
#include <list>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <stdlib.h> 
using namespace std;

void ReadGraph(GraphType *graph, char *filename, int *wgtflag)
{
  int i, j, k, l, edge;
  idxtype *xadj, *adjncy, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  InitGraph(graph);

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    free(line);
    return;
  }

  sscanf(line, "%d %d", &(graph->nvtxs), &(graph->nedges));

  graph->nedges *=2;
  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  // Joyce modified this part to correctly construct a coarsend graph (edge weight)
  *wgtflag = 1;
  //vwgt = graph->vwgt = idxsmalloc(graph->nvtxs, 1, "RefineKWay: graph->vwgt");
  adjwgt = graph->adjwgt = idxsmalloc(graph->nedges, 1, "RefineKWay: graph->adjwgt");

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<graph->nvtxs; i++) {
    do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    for (;;) {
      edge = (idxtype)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (edge < 0)
        break;

      adjncy[k] = edge;
      k++;
    } 
    xadj[i+1] = k;
  }

  fclose(fpin);

  if (k != graph->nedges) {
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("In the first line of the file, you specified that the graph contained\n%d edges. However, I only found %d edges in the file.\n", graph->nedges/2, k/2);
    if (2*k == graph->nedges) {
      printf("\n *> I detected that you specified twice the number of edges that you have in\n");
      printf("    the file. Remember that the number of edges specified in the first line\n");
      printf("    counts each edge between vertices v and u only once.\n\n");
    }
    printf("Please specify the correct number of edges in the first line of the file.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  free(line);

}


void WriteOverlapClust(char *fname, std::vector< list<int> >& ClustResults, int n, int nparts, float alpha, float beta, float sigma)
{
  FILE *fpout;
  int i;
  char filename[256];
  char extracted[256];

  extractfilename(fname, extracted);
  sprintf(filename,"%s_clust_%d_alpha_%1.2f_beta_%1.4f_sigma_%1.2f", extracted, nparts, alpha, beta, sigma);

  ofstream myfile;
  myfile.open (filename);

  //cout << "\n ----- print ClustResults... WriteOverlapClust\n";
  for(vector< list<int> >::iterator iter = ClustResults.begin(); iter != ClustResults.end(); iter++){
    for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
     //cout << *i << " ";  
     myfile << *i << " ";  
    }
    //cout<<"\n";
    myfile<<"\n";
  }

  myfile.close();
}


void ReadInitialization(char *fname, int numNodes, std::vector< list<int> >& AssVec){

  string line;
  ifstream myfile (fname);

  if (myfile.is_open())
  {
    for(int node=0; node<numNodes; node++){
       getline (myfile,line);
       istringstream buf(line);
       istream_iterator<std::string> beg(buf), end;
       vector<std::string> tokens(beg, end); 
       for(std::vector<std::string>::iterator i=tokens.begin(); i!=tokens.end(); i++){
	  string t = *i;
	  AssVec[node].push_back(atoi(t.c_str())-1); //FIXME (Fixed...)
       }
    }
    myfile.close();

  }else{
    cout << "!!!ERROR: Unable to open file\n"; 
    exit(0);
  }

  // print
  /*
  cout << "convert to vector-list form...\n";
  for(int iter = 0; iter < numNodes; iter++){
    cout << "node no.: " << iter << " ... ";
    for(list<int>::iterator i=AssVec[iter].begin(); i!=AssVec[iter].end(); i++){
	cout << *i << " ";
    }
    cout << "\n";
  }
  */
}

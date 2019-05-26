/*
 * Written by Joyce Whang
 * Overlapping Non-exhaustive Graph Clustering
 * Jan. 9, 2014
 *
 */

#include <metis.h>
#include <float.h>
#include <iostream> 
#include <vector>
#include <list>
#include <queue>
#include <assert.h>
#include "proto_ovfuncs_wkkm.h"

using namespace std;

//void Compute_Weights(CtrlType *ctrl, GraphType *graph, idxtype *w)
void Compute_Weights(GraphType *graph, idxtype *w)
{
  int nvtxs, i, j;
  idxtype *xadj, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjwgt = graph->adjwgt;

  //if ((cutType == RASSO) || (cutType == RCUT))
    //for (i=0; i<nvtxs; i++)
      //w[i] = 1;
  //else
    if (adjwgt == NULL)
      for (i=0; i<nvtxs; i++)
	for (j=xadj[i]; j<xadj[i+1]; j++) 
	  w[i] ++;
    else
      for (i=0; i<nvtxs; i++)
	for (j=xadj[i]; j<xadj[i+1]; j++) 
	  w[i] += adjwgt[j];
}



void MLKKMRefine_overlap(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int nparts, int chain_length, float *tpwgts, float ubfactor, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma)
{
  int nlevels, mustfree=0;
  GraphType *ptr;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->UncoarsenTmr));

  // Determine how many levels are there 
  for (ptr=graph, nlevels=0; ptr!=orggraph; ptr=ptr->finer, nlevels++); 
  //cout << "Number of levels is "<< nlevels <<"\n";

  // to remember the previous level assignments
  vector< list<int> > coarser_AssVec;

  for (int level=0; ;level++) {
    timer tmr;
    float result;

    cleartimer(tmr);
    starttimer(tmr);

    // cluster assignment vector
    vector< list<int> > AssVec(graph->nvtxs);
    //cout << "No. of vertices in level " << level << ": " << graph->nvtxs << "\n";

    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));
    // for the coarsest level, convert where to AssVec 
    if(level==0){
      idxtype *where = graph->where;
      for (int it=0; it<graph->nvtxs; it++){
	AssVec[it].push_back(where[it]);
      }
      // print initial partitioning result
      //cout << "\n ----- print where...\n";
      for (int a=0; a<graph->nvtxs; a++) {
        int k = where[a];
        //cout << "where[" << a << "]: " << k << "\n";
      }
    }else{ // Projection
      // print cmap
      //cout << "\n print cmap...\n";
      idxtype *cmap = graph->cmap;
      // initialization of the current level assignment
      for (int p=0; p<graph->nvtxs; p++){
	int k = cmap[p]; // k: previous level node, p: current level node
	//cout << "\n cmap[" << p << "]: " << k << "\n";
  	//cout << "coarser_AssVec[" << k << "] \n";
	for(list<int>::iterator i = coarser_AssVec[k].begin(); i!=coarser_AssVec[k].end(); i++){
		AssVec[p].push_back(*i);
		//cout<<" "<<*i<<" ";		
	}
      }
    }
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));

    /*
    cout << "\n ----- print AssVec... level: "<<level<<"\n";
    for(vector< list<int> >::iterator iter = AssVec.begin(); iter != AssVec.end(); iter++){
	for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
		cout << *i << " ";    
        }
    cout<<"\n";
    }
    */

    idxtype *w = idxsmalloc(graph->nvtxs, 0, "pingpong: weight");
    Compute_Weights(graph, w);

    if (graph == orggraph){
      One_level_overlap_refine(graph, nparts, w, AssVec, alpha, beta, sigma);
      ClustResults = AssVec;
      /*
      cout << "\n ----- print ClustResults... MLKKMRefine_overlap\n";
      for(vector< list<int> >::iterator iter = ClustResults.begin(); iter != ClustResults.end(); iter++){
        for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
	  cout << *i << " ";    
        }
        cout<<"\n";
      }
      */
      break;
    }
    else{
      One_level_overlap_refine(graph, nparts, w, AssVec, alpha, beta, sigma);
    }

    GKfree((void **) &graph->gdata, LTERM);  // Deallocate the graph related arrays 

    // for projection, remember the previous level results
    coarser_AssVec = AssVec;
    
    graph = graph->finer;

    if (graph->vwgt == NULL) {
      //cout<< i <<" ((((((((graph->vwgt == NULL\n";
      graph->vwgt = idxsmalloc(graph->nvtxs, 1, "RefineKWay: graph->vwgt");
      graph->adjwgt = idxsmalloc(graph->nedges, 1, "RefineKWay: graph->adjwgt");
      mustfree = 1;
    }//else
      //cout<< i <<" ((((((((graph->vwgt != NULL\n";
  }
  
  if (mustfree) 
    GKfree((void **) &graph->vwgt, (void **) &graph->adjwgt, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->UncoarsenTmr));
  //cout<<"End of MLKKMRefine\n";
}



void One_level_overlap_refine(GraphType *graph, int nparts, idxtype *w, std::vector< list<int> >& AssVec, float alpha, float beta, float sigma){
// w is the weights

  int nvtxs = graph->nvtxs;
  idxtype * xadj = graph->xadj;
  idxtype * adjncy = graph->adjncy;
  idxtype * adjwgt = graph->adjwgt;

  /*
  std::cout << "One level overlap refine...\n";
  int tmp=0;
  for(vector< list<int> >::iterator iter = AssVec.begin(); iter != AssVec.end(); iter++){
	tmp++;
	for(list<int>::iterator i= (*iter).begin(); i!=(*iter).end(); i++){
		cout <<"node: "<<tmp<<", clust: "<< *i << " ";
	}
	cout<<"\n";
  }
  */
 
  idxtype * clust_sum = idxsmalloc(nparts,0, "One_level_overlap_refine: weight sum");
  float * clust_inv_sum = fmalloc(nparts, "One_level_overlap_refine: sum inverse"); 
  float * clust_squared_inv_sum = fmalloc(nparts, "One_level_overlap_refine: squared sum inverse"); 
  idxtype * clust_squared_sum = idxsmalloc(nparts,0, "One_level_overlap_refine: weight squared sum");
 
  // while loop conditions
  float obj, old_obj;
  old_obj = 0;
  obj = FLT_MAX;
  float epsilon = 0.00001;
  int currit=0;
  //float sigma = 1;
  //float alpha = 0.1; // 0.1
  //float beta = 0.02; // 0.02
  int alphaN = round(nvtxs * alpha);
  int betaN = round(nvtxs * beta);
  std::cout << "sigma: " << sigma << " nvtxs: " << nvtxs << " alpha: " << alpha << " alphaN: " << alphaN << " beta: " << beta << " betaN: " << betaN <<"\n";
  //std::cout << "nvtxs: " << nvtxs << " beta: " << beta << " betaN: " << betaN <<"\n";
  idxtype * clust_linearTerm = idxmalloc(nparts, "One_level_overlap_refine: linearTerm");
  idxtype * clust_in = idxmalloc(nparts, "One_level_overlap_refine: clust_in");
  vector< list<int> > backup_AssVec(graph->nvtxs);
  // --------------------------------

  while(fabs(old_obj-obj) > epsilon && currit < MAXITERATIONS){

    old_obj = obj;
    currit++;
    priority_queue<Distance, vector<Distance>, CompareDistance> ex_queue; // exhaustive
    priority_queue<Distance, vector<Distance>, CompareDistance> ov_queue; // overlapping
    myQueue mq;
    
    //cout<< "----- currit: " << currit << "\n";

    // Compute cluster constant
    for (int i=0; i<nparts; i++){
      clust_sum[i] = 0;
      clust_squared_sum[i] = 0;
    }
    for (int i=0; i<nvtxs; i++){
      for(list<int>::iterator c = AssVec[i].begin(); c!=AssVec[i].end(); c++){
	  clust_sum[*c] += w[i];
      }
    }
    for (int i=0; i<nparts; i++){
      if(clust_sum[i] >0){
        clust_inv_sum[i] = 1.0/clust_sum[i];
        clust_squared_inv_sum[i] = clust_inv_sum[i]*clust_inv_sum[i];
      }
      else{
        clust_inv_sum[i] = clust_squared_inv_sum[i] = 0;
      }
    }
    for (int i=0; i<nvtxs; i++){
      //cout << "node no.: " << i << " ... ";
      for(list<int>::iterator c = AssVec[i].begin(); c!=AssVec[i].end(); c++){
        int me = *c;
	//cout << me << " ";
        for (int j=xadj[i]; j<xadj[i+1]; j++){
	  int neighbor = adjncy[j];
          for(list<int>::iterator nc = AssVec[neighbor].begin(); nc!=AssVec[neighbor].end(); nc++){
	    if(*nc == me){
	      clust_squared_sum[me] += adjwgt[j];
            } 
          } 
        } 
      } 
      //cout << "\n";
    } 

     /*
    // print
    for (int i=0; i<nparts; i++){
      if (clust_sum[i]==0){
        cout<<"clust_sum["<<i<<"]: "<<clust_sum[i]<<"\n";
      }
    }
    */

    /*
    // print
    for(int iter = 0; iter < nvtxs; iter++){
	cout << "node no.: " << iter << " ... ";
	for(list<int>::iterator i=AssVec[iter].begin(); i!=AssVec[iter].end(); i++){
		cout << *i << " ";
	}
	cout << "\n";
    }
    cout << "sigma: " << sigma <<"\n";
    for (int i=0; i<nparts; i++){
      cout << "clust_squared_sum["<<i<<"]: "<< clust_squared_sum[i]<<", clust_squared_inv_sum["<<i<<"]: "<<clust_squared_inv_sum[i]<<", clust_inv_sum["<<i<<"]: "<<clust_inv_sum[i]<< "...... clust_sum["<<i<<"]: "<< clust_sum[i] <<", clust_inv_sum["<<i<<"]: "<<clust_inv_sum[i]<<"\n";
    }
    cout<<"dist = clust_squared_sum[k]*clust_squared_inv_sum[k] - 2*inv_wi*clust_linearTerm[k]*clust_inv_sum[k] + sigma*inv_wi +- sigma*clust_inv_sum[k];\n";
    */
    // ----------------------------------------
    

    // Consider all the nodes
    for (int i=0; i<nvtxs; i++){

      if(w[i]>0){
	float inv_wi=1.0/w[i];
        for (int k=0; k<nparts; k++){
	  clust_linearTerm[k] = 0;
	  clust_in[k] = 0;
	}
        for (int j=xadj[i]; j<xadj[i+1]; j++){
	  int neighbor = adjncy[j];
          for(list<int>::iterator nc = AssVec[neighbor].begin(); nc!=AssVec[neighbor].end(); nc++)
	    clust_linearTerm[*nc] += adjwgt[j];
        }	
	for (int k=0; k<nparts; k++){
	  for(list<int>::iterator c = AssVec[i].begin(); c!=AssVec[i].end(); c++){
	    if (k == *c){
	       clust_in[k]=1;
	    }
	  }
	}
	float clust_dist[nparts];
	Distance d[nparts];
	int min_cind = 0;
	float min_cdist = FLT_MAX;
	for (int k=0; k<nparts; k++){
	  float dist;
	  if (clust_in[k]==1){
	     //cout<< "i: " << i << ", k: " << k << ", clust in \n";
	     dist = clust_squared_sum[k]*clust_squared_inv_sum[k] - 2*inv_wi*clust_linearTerm[k]*clust_inv_sum[k] + sigma*inv_wi - sigma*clust_inv_sum[k]; 
	  }else{
	     //cout<< "i: " << i << ", k: " << k << ", clust out \n";
	     dist = clust_squared_sum[k]*clust_squared_inv_sum[k] - 2*inv_wi*clust_linearTerm[k]*clust_inv_sum[k] + sigma*inv_wi + sigma*clust_inv_sum[k];
	  }
	  clust_dist[k] = dist;
	  if(dist < min_cdist){
	    min_cdist = dist;
	    min_cind = k;
	  }
	}
	for (int k=0; k<nparts; k++){
	  d[k].dist = clust_dist[k];
	  d[k].nodeNum = i;
	  d[k].clustNum = k;
	  if(k==min_cind){
	    ex_queue.push(d[k]); 
	  }else{ 
	    mq.insert_dist(d[k],ov_queue,betaN+alphaN);
	  }
	  //cout << "i: " << i << ", k: " << k << ", clust_linearTerm[k]: " << clust_linearTerm[k] <<", inv_wi: "<< inv_wi << ", dist: " << d[k].dist << "\n";
	}
      } // endif w[i] >0
    } // endfor
    // -------------------------------
    
    // adjust two queues
    for(int p=0; p<betaN; p++){
      Distance d = ex_queue.top();
      ex_queue.pop();
      mq.insert_dist(d,ov_queue,betaN+alphaN);
    }
   
    // Update assignments
    VectorList vl;
    backup_AssVec = AssVec;
    obj=0;
    vl.convert2VecList(AssVec, nvtxs, ex_queue, ov_queue, obj);
    /*
    std::cout << "convert to vector-list form...\n";
    for(int iter = 0; iter < nvtxs; iter++){
	cout << "node no.: " << iter << " ... ";
	for(list<int>::iterator i=AssVec[iter].begin(); i!=AssVec[iter].end(); i++){
		cout << *i << " ";
	}
	cout << "\n";
    }
    */

    cout << "- iteration: "<<currit<<", obj: "<<obj<<"\n";
    // if matrix is not positive definite
    //if (obj > old_obj && currit>2){
    if (obj > old_obj){
        cout << "*** roll back assignments *** currit: "<<currit<<"\n";
	AssVec = backup_AssVec;
        break;
    }
    
    /*
    printf("--- obj: %6.6f, old_obj: %6.6f, abs(obj-old_obj) = %6.6f\n", obj, old_obj, fabs(obj-old_obj));
    cout << "--- currit: " << currit << " MAXITERATIONS: " << MAXITERATIONS << "\n";
    cout << "--- fabs(old_obj-obj): " << fabs(old_obj-obj) << " epsilon: " << epsilon << "\n"; 
    */

  } // endwhile
  //cout << "done with "<<currit<<" iterations... obj: "<<obj<<"\n";

   /*
   std::cout << "!!! final assignment...\n";
   for(int iter = 0; iter < nvtxs; iter++){
     cout << "node no.: " << iter << " ... ";
     for(list<int>::iterator i=AssVec[iter].begin(); i!=AssVec[iter].end(); i++){
	cout << *i << " ";
     }
     cout << "\n";
   }
   */

  free(clust_sum); free(clust_squared_sum); free(clust_linearTerm); free(clust_inv_sum); free(clust_squared_inv_sum);
}



using namespace std;

void One_level_overlap_refine(GraphType *graph, int nparts, idxtype *w, std::vector< std::list<int> >& AssVec, float alpha, float beta, float sigma); 

struct Distance {
    float dist; 
    int nodeNum; 
    int clustNum; 
};

class CompareDistance {
    public:
      bool operator()(Distance& d1, Distance& d2)
      {
        if (d1.dist < d2.dist) return true;
          return false;
      }
};

class myQueue {
    public:
      void insert_dist (Distance& d, priority_queue<Distance, vector<Distance>, CompareDistance>& pq, int queue_size);
};

class VectorList {
    public:
      void convert2VecList (vector< list<int> >& AssVec, int nvtxs, priority_queue<Distance, vector<Distance>, CompareDistance>& ex_queue, priority_queue<Distance, vector<Distance>, CompareDistance>& ov_queue, float &obj);
};

void myQueue::insert_dist (Distance& d, priority_queue<Distance, vector<Distance>, CompareDistance>& pq, int queue_size){

    int currSize = pq.size();
    //cout<< " currSize: " << currSize <<"\n";
    if (currSize < queue_size){
       pq.push(d);
       //cout<< "myQueue::insert_dist inserted..." << d.dist <<" "<< d.nodeNum << " "<< d.clustNum <<"\n";
    }else{
       if (d.dist < pq.top().dist){ 
	//cout<< "myQueue::insert_dist removed..." << pq.top().dist <<" "<< pq.top().nodeNum << " "<< pq.top().clustNum <<"\n";
	pq.pop();
	pq.push(d);
	//cout<< "myQueue::insert_dist inserted..." << d.dist <<" "<< d.nodeNum << " "<< d.clustNum <<"\n";	
       }
    }   
}

void VectorList::convert2VecList (vector< list<int> >& AssVec, int nvtxs, priority_queue<Distance, vector<Distance>, CompareDistance>& ex_queue, priority_queue<Distance, vector<Distance>, CompareDistance>& ov_queue, float &obj){

    assert (AssVec.size()==nvtxs);

    for(int i=0; i<nvtxs; i++){
	AssVec[i].clear();
    }

    obj = 0;
    //cout<< "==== ex_queue.size(): " << ex_queue.size() << "\n";
    while (!ex_queue.empty()) {
       Distance d = ex_queue.top();
       //cout << d.dist << " " << d.nodeNum << " " << d.clustNum << endl;
       AssVec[d.nodeNum].push_back(d.clustNum);
       obj = obj + d.dist;
       //cout<< "==== obj: " << obj << "\n";
       ex_queue.pop();
    }
    
    //cout<< "==== ov_queue.size(): " << ov_queue.size() << "\n";
    while (!ov_queue.empty()) {
       Distance d = ov_queue.top();
       //cout << d.dist << " " << d.nodeNum << " " << d.clustNum << endl;
       AssVec[d.nodeNum].push_back(d.clustNum);
       obj = obj + d.dist;
       //cout<< "==== obj: " << obj << "\n";
       ov_queue.pop();
    }
    //cout<< "==== VectorList::convert2VecList... obj: " << obj << "\n";
}

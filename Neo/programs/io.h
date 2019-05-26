using namespace std;
void ReadGraph(GraphType *graph, char *filename, int *wgtflag);
void WriteOverlapClust(char *fname, std::vector< list<int> >& AssVec, int n, int nparts, float alpha, float beta, float sigma);
void ReadInitialization(char *fname, int numNodes, std::vector< list<int> >& AssVec);

using namespace std;

void MLKKM_PartGraphKway_overlap(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *chainlength,  int *options, int *edgecut, idxtype *part, int levels, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma);

void MLKKM_WPartGraphKway_overlap(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *chainlength, float *tpwgts, int *options, int *edgecut, idxtype *part, int levels, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma);

int MLKKMPartitioning_overlap(CtrlType *ctrl, GraphType *graph, int nparts, int chain_length, idxtype *part, float *tpwgts, float ubfactor, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma);

void MLKKMRefine_overlap(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int nparts, int chain_length, float *tpwgts, float ubfactor, std::vector< std::list<int> >& ClustResults, float alpha, float beta, float sigma);

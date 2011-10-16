#include <vector>
#include "kodtree.h"
using namespace std;

class PointInPolyhedron{

public:
	int isPinPolyhedron(double p[3]); //-1,0,1
	PointInPolyhedron(double (*vti)[3], int numvi,int (*tris)[3],int numtris);//,double epsi);
	double absoluteClosestSqDistance(Point p){
		CellNode3D *cnode=polytree->findaLeafCellContainingPoint(polytree->getRoot(),p);
		int ntri; double dist;
		getAbsoluteClosestTriForPointInGCell(p,cnode,ntri,dist);
		return dist;
	}
	~PointInPolyhedron();
private:
	int testPinPolyhedronForPinGcell(double p[3],CellNode3D *cnode);
	void getRelativeClosestEntityForPointInGCell( double p[3],CellNode3D *cnode,int &id,int &nentity, int &ntri,double &dist);
	void getAbsoluteClosestTriForPointInGCell(double p[3],CellNode3D *cnode, int &tri, double &dist);
	void getTheClosestTriNonLeaf(double p[3],double dist0,CellNode3D *pcell,int &tri,double &dist);
	void getTheClosestTriAmongCell(double p[3],CellNode3D *pcell, double &dist,int &ntri);
	void getEndPointOfTri(int tri, double p0[3],double p1[3],double p2[3]);
	void setGCellAttribOfSubTree(CellNode3D *pcell);
	void getCellSeqWithUnknownAttribFromaCell(CellNode3D *cnode,vector<CellNode3D *> *&pcellseq,
					CellNode3D * &pcellm,int &ia,double pm[3]);
	int classifyVert(double p[3],int vert);
	int  classifyEdge(int tria,int trib);
	void  getVertsAroundaVert(int v, int * &nbverts,int &numnbv);
	int  nextTriOfVert(int v, int ctri);
	int  nextVertOfTri(int tri,int v);
	void  getThePointFormingLeastAngleWith2Points(double p[3],int v, int *nbverts,
		int numnbv,double &maxcosa,int &vridge);
	void  get2TriCom2Vert(int va, int vb, int &ta, int &tb);
	void  formNeighbAndTriOfNode();
	int indexOfNeighbTriToTri(int tria,int trib);
	int indexOfVertAtTri(int v, int ctri);
	int neighbOfTri(int tri,int ind);
	void wrapPointsUpasVerts(void  ** &vti);
	void recNeighbOfTrips(int *numtriofnode,int *tripositionofnode,int *trilist);
	void  getNeighbFromTrilist(int tri,int ind,int &tnb,int &nbindex
								,int *numtriofnode,int *tripositionofnode,int *trilist);
	void  getEdgeOfTri(int np[3], int index, int &a, int &b);
//	void recoverTriused(CellNode3D *pcell);
	static void pofvforcoordnodes3(double p[3],void *pv);
	static bool ifexinfooverlapbox(void *info,int infotype,const Box &bd,double eps);
	static bool ifexinfoshouldbeincell(void *info,int infotype,CellNode3D *cnode);
private:
	static const double epsilonon;
	static const double epsoverlap;
	static const double epscoplanar;
	//double epscell;
	Kodtree *polytree;
	static double (*vertcoord)[3];
	static int numvert;
	static int (*trips)[3],numtri;
	int (*tneighb)[3];
	int *triofnode;
	int *vertattrib;
//	int *triused;
};

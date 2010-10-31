#include <vector>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <math.h>

//----------------------------------------------------------------
// Include files for R functions such as Rprintf()

#include "R.h" // R functions
#include "Rmath.h" // R math

//----------------------------------------------------------------

#include <float.h> // For DBL_MAX

#include "pinpolyhedron.h"
extern  void jf_error(char *);
using namespace std;

const double PointInPolyhedron::epsilonon=0.00000000000001;
const double PointInPolyhedron::epsoverlap=0.000001;
const double PointInPolyhedron::epscoplanar=0.000001;
double (*PointInPolyhedron::vertcoord)[3];
int PointInPolyhedron:: numvert;
int (*PointInPolyhedron::trips)[3];
int PointInPolyhedron::numtri;
int absolute;
int *startaddress=(int *)1;
extern int positionOfPointProjectToTri(double p[3],double p0[3],double p1[3],double p2[3]);
extern double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]);
extern double sqDistPointToSeg3D(double p[3],double p0[3],double p1[3]);
int triIndexFromPt(void *ptri){
	int *ptr=(int *)ptri;
	return ptr-startaddress;
}
int vertIndexFromPt(void *pv){
	int *pt=(int *)pv;
	return pt-startaddress;
}

void  PointInPolyhedron::pofvforcoordnodes3(double p[3],void *pv){

//	extern static int *startaddress;

	int nd=vertIndexFromPt(pv);
	p[0]=vertcoord[nd][0];
	p[1]=vertcoord[nd][1];
	p[2]=vertcoord[nd][2];
}
void PointInPolyhedron::wrapPointsUpasVerts(void  ** &vti){

	vti=new void *[numvert];
	for(int i=0; i<numvert; i++)
		vti[i]=startaddress+i;
}

bool   PointInPolyhedron::ifexinfooverlapbox(void *info,int infotype,const Box &bd,double eps){
	  if(infotype==1){
		  int tri =triIndexFromPt(info);
		  return isTriangleBoxOver(vertcoord[trips[tri][0]],vertcoord[trips[tri][1]],vertcoord[trips[tri][2]],bd,eps);
	  }
	  return false;
  }

bool   PointInPolyhedron::ifexinfoshouldbeincell(void *info,int infotype,CellNode *cnode){
	  if(infotype==1){
		  int tri=triIndexFromPt(info);
		  for(int i=0; i<cnode->numvert; i++){
			  int v=vertIndexFromPt(cnode->vert[i]->vt);
			  if(v==trips[tri][0]||v==trips[tri][1]||v==trips[tri][2])
				  return false;
		  }
	  }
	  return true;
}
int PointInPolyhedron::isPinPolyhedron( double p[3]){

//	int rt;
	CellNode *pcell;
//	vector<CellNode *> *pcellseq;

	pcell=polytree->findaLeafCellContainingPoint(polytree->getRoot(),p);
	if(pcell==0) 
		return -1;
	if(pcell->inoutattrib==-1||pcell->inoutattrib==1)
		return pcell->inoutattrib;
	else if(pcell->inoutattrib ==0)
		return testPinPolyhedronForPinGcell(p,pcell);
	//else if(pcell->inoutattrib==-2)
	if((pcell->inoutattrib=testPinPolyhedronForPinGcell(p,pcell))==0)
		jf_error("err ispointin");
	else return pcell->inoutattrib;
/*	double pm[3];
	CellNode *pcellm=0;
	getCellSeqWithUnknownAttribFromaCell(pcell,pcellseq,pcellm,rt,pm);
	if(rt==0)
		rt=testPinPolyhedronForPinGcell(pm,pcellm);
	if(rt==0)
		jf_error("ispinoPolyhedron");
	if(pcellseq!=0) // what's the meaning, will it be colored when only pcellm alone?
		for(unsigned i=0; i<pcellseq->size(); i++)
			(*pcellseq)[i]->inoutattrib=rt;
	delete pcellseq;
	return rt;
*/
}

int PointInPolyhedron::testPinPolyhedronForPinGcell(double p[3],CellNode *cnode){

	int id,nentity,tri,rt;
	double dist,p0[3],p1[3],p2[3];

	getRelativeClosestEntityForPointInGCell(p,cnode,id,nentity,tri,dist);
	if(dist<=epsilonon)
		return 0;
	if(id==0){
		if(vertattrib[nentity]==-1||vertattrib[nentity]==1)
			return vertattrib[nentity];
		else if(vertattrib[nentity]==2) 
			jf_error("err testpinpolyh");
		else if(vertattrib[nentity]==-2){
			int rt=classifyVert(p,nentity);
			if(rt==-1||rt==1) return rt;
			if(rt==2)
				jf_error("err testpinpolyh1");
		}
	//	nentity=tri; //attrib==0||attrib==-2&&rt==0;
	}
	if(id==1){
		rt=classifyEdge(nentity,tri);
		if(rt==-1||rt==1)
			return rt;
	//	nentity=tri;
	}
	if(id!=1&&id!=0&&id!=2)
		jf_error("err ispoinPolyhedron");
	getEndPointOfTri(tri,p0,p1,p2); //id==2||id==1&&coplanar at edge||id==0&&coplanar at vertex
	if(VolumOf4p(p0,p1,p2,p)<0)
		return 1;
	else 
		return -1;
}


void PointInPolyhedron::getRelativeClosestEntityForPointInGCell( double p[3],CellNode *cnode,int &id,
														   int &nentity, int &ntri,double &dist){

	int ip;
	double p0[3],p1[3],p2[3];

//	if(absolute==0)
//		getRelativeClosestTriForPointInGCell(p,cnode,ntri,dist);
//	else
		getAbsoluteClosestTriForPointInGCell(p,cnode,ntri,dist);
	if(dist==numeric_limits<double>::max())
		jf_error("err getrelativeclosetentityforpingcell");
	getEndPointOfTri(ntri,p0,p1,p2);
	if((ip=positionOfPointProjectToTri(p,p0,p1,p2))==6){
		nentity=ntri;
		id=2;
	}else if(ip<3){
		nentity=trips[ntri][ip];
		id=0;
	}else{
		nentity=neighbOfTri(ntri,ip-3);
		if(nentity<0) jf_error("getrealvie");
		id=1;
	}
}

void PointInPolyhedron::getAbsoluteClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist){

	CellNode *pcell0=0;
	CellNode *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosesttri");
	dist=numeric_limits<double>::max();
	tri=-1;
	
	while(pcell){
		int trin;
		double distn;
		getTheClosestTriNonLeaf(p,dist,pcell->anotherChild(pcell0),trin,distn);
		if(distn<dist){
			dist=distn; tri=trin;
		}
		if(sqdistInnerPointToBoxBound(p,pcell->bound)>=dist)
			return;
		pcell0=pcell;
		pcell=pcell->parent;
	}
}


void PointInPolyhedron::getTheClosestTriNonLeaf(double p[3],double dist0,
										   CellNode *pcell,int &tri,double &dist){
	
	double distn;
	int trin;

	dist=dist0, tri=-1;
	if(sqdistPointToBox(p,pcell->bound)>=dist0) return;
	if(pcell->isLeaf()){
		getTheClosestTriAmongCell(p,pcell,distn,trin);
		if(distn<dist){	dist=distn;	tri=trin; return;}
	}else{
		CellNode *sortsub[2]={pcell->child[0],pcell->child[1]};
		if(sqdistPointToBox(p,pcell->child[0]->bound)>sqdistPointToBox(p,pcell->child[1]->bound)){
			sortsub[0]=pcell->child[1];
			sortsub[1]=pcell->child[0];
		}
		for(int i=0; i<2; i++){
			getTheClosestTriNonLeaf(p,dist,sortsub[i],trin,distn);
			if(distn<dist){ dist=distn; tri=trin; }
		}
	}
}
/*void PointInPolyhedron::recoverTriused(CellNode *pcell){
	
	int tri;
	if(pcell->lpwpinfo!=0)
		for(std::list<WpInfo *>::iterator ite=pcell->lpwpinfo->begin();ite!=pcell->lpwpinfo->end(); ite++){
			if((*ite)->infotype!=1)	continue;
			tri=triIndexFromPt((*ite)->info);
			triused[tri]=0;
		}
//	if(pcell->numvert<=0)
//		return;
	for(int i=0; i<pcell->numvert; i++){
		int v=vertIndexFromPt(pcell->vert[i]->vt);
		int tri0,tri;
		tri0=tri=triofnode[v];
		do{
			triused[tri]=0;
		}while((tri=nextTriOfVert(v,tri))!=tri0);
	}
}*/
void PointInPolyhedron::getTheClosestTriAmongCell(double p[3],CellNode *pcell, double &dist,int &ntri){

	int tri;
	double distemp,p0[3],p1[3],p2[3];

	dist=numeric_limits<double>::max();
	if(!pcell||!pcell->isLeaf())
		jf_error("error gettheclosettriamongcell");
	if(pcell->lpwpinfo!=0)
		for(std::list<WpInfo *>::iterator ite=pcell->lpwpinfo->begin();ite!=pcell->lpwpinfo->end(); ite++){
			if((*ite)->infotype!=1)	continue;
			tri=triIndexFromPt((*ite)->info);
//			triused[tri]=1;
			getEndPointOfTri(tri,p0,p1,p2);
			if((distemp=sqDistPointToTri(p,p0,p1,p2))<dist){
				dist=distemp;
				ntri=tri;
			}
		}
//	if(pcell->numvert<=0)
//		return;
	for(int i=0; i<pcell->numvert; i++){
		int v=vertIndexFromPt(pcell->vert[i]->vt);
		int tri0,tri;
		tri0=tri=triofnode[v];
		do{
//			if(triused[tri]==1) continue;
//			else triused[tri]=1;
			getEndPointOfTri(tri,p0,p1,p2);
			if((distemp=sqDistPointToTri(p,p0,p1,p2))<dist){
				dist=distemp;
				ntri=tri;
			}
		}while((tri=nextTriOfVert(v,tri))!=tri0);
	}
//	recoverTriused(pcell);
}


void PointInPolyhedron::getEndPointOfTri(int tri, double p0[3],double p1[3],double p2[3]){

	p0[0]=vertcoord[trips[tri][0]][0];
	p0[1]=vertcoord[trips[tri][0]][1];
	p0[2]=vertcoord[trips[tri][0]][2];
	p1[0]=vertcoord[trips[tri][1]][0];
	p1[1]=vertcoord[trips[tri][1]][1];
	p1[2]=vertcoord[trips[tri][1]][2];
	p2[0]=vertcoord[trips[tri][2]][0];
	p2[1]=vertcoord[trips[tri][2]][1];
	p2[2]=vertcoord[trips[tri][2]][2];
}
int positionOfPointProjectToSeg3D(double p[3],double p0[3],double p1[3]){

	double vp0p[3],vp0p1[3],vp1p[3];

	vec_2p(p0,p,vp0p);

	vec_2p(p0,p1,vp0p1);
	if(vec_dotp(vp0p,vp0p1)<=0)
		return -1;
	vec_2p(p1,p,vp1p);
	if(vec_dotp(vp1p,vp0p1)>=0)
		return 1;
	return 0;  
}
int positionOfPointProjectToTri(double p[3],double p0[3],double p1[3],double p2[3]){

	double v0p[3],v20[3],v01[3];
	vec_2p(p0,p,v0p);
	vec_2p(p2,p0,v20);
	vec_2p(p0,p1,v01);

	double d0p20=vec_dotp(v0p,v20);
	double d0p01=vec_dotp(v0p,v01);
	if(d0p20>=0&&d0p01<=0) return 0;
	
	double v1p[3],v12[3];
	vec_2p(p1,p,v1p);
	vec_2p(p1,p2,v12);
	double d1p01=vec_dotp(v1p,v01);
	double d1p12=vec_dotp(v1p,v12);
	if(d1p01>=0&&d1p12<=0) return 1;

	double v2p[3];
	vec_2p(p2,p,v2p);
	double d2p12=vec_dotp(v2p,v12);
	double d2p20=vec_dotp(v2p,v20);
	if(d2p12>=0&&d2p20<=0) return 2;

	double nm012[3],nm01p[3],nm12p[3],nm20p[3];
	vec_crop(v20,v01,nm012);

	vec_crop(v01,v0p,nm01p);
	double dt01=vec_dotp(nm012,nm01p);
	if(dt01<=0&&d0p01>=0&&d1p01<=0)
		return 5; // rt==0;
		
	vec_crop(v12,v1p,nm12p);
	double dt12=vec_dotp(nm012,nm12p);
	if(dt12<=0&&d1p12>=0&&d2p12<=0)
		return 3; // rt==0;
	
	vec_crop(v20,v2p,nm20p);
	double dt20=vec_dotp(nm012,nm20p);
	if(dt20<=0&&d2p20>=0&&d0p20<=0)
		return 4; // rt==0;

	if(dt01>0&&dt12>0&&dt20>0) return 6;
	else jf_error("asdf posiotin");
}
extern void sortTrianglesOuterNormAndRecNeighb(double (*vertcoord)[3],int numvert,int (*trips)[3],
										int numtri,int (*tneighb)[3],int *triofnode);
PointInPolyhedron::PointInPolyhedron(double (*vti)[3], int numvi,int (*tris)[3],int numti){//,double epsion=0){

//	epsilonon=epsion;

	numvert=numvi;
	vertcoord=(double (*)[3]) new double[3*numvert];
	memcpy(vertcoord,vti,sizeof(double)*3*numvert);
	numtri=numti;
	trips=(int (*)[3])new int[3*numtri];
	memcpy(trips, tris,sizeof(int)*3*numtri);

	tneighb=(int (*)[3]) new int [3*numtri];
	triofnode=(int *) new int[numvert];
	vertattrib= new int [numvert];
	for(int i=0; i<numvert; i++)
		vertattrib[i]=-2;
//	formNeighbAndTriOfNode();
	sortTrianglesOuterNormAndRecNeighb(vertcoord,numvert,trips,numtri,tneighb,triofnode);
	void **wvti;
	wrapPointsUpasVerts(wvti);
	polytree=new Kodtree(wvti,numvert,pofvforcoordnodes3,3,epsoverlap);
	delete wvti;
    polytree->setFuncExinfoShouldbeInCell(ifexinfoshouldbeincell);
	polytree->setFuncExinfoOverlapBox(ifexinfooverlapbox);
	for(int i=0; i<numtri; i++)
		polytree->insertExinfo(i+startaddress,1);
	setGCellAttribOfSubTree(polytree->getRoot());
	//triused=new int[numtri];
	//for(int i=0; i<numtri; i++)
	//	triused[i]=0;
}

void PointInPolyhedron::setGCellAttribOfSubTree(CellNode *pcell){

	if(!pcell) return;
	if(!pcell->isLeaf())
		for(int i=0; i<2; i++)
			setGCellAttribOfSubTree(pcell->child[i]);
	else if(pcell->lpwpinfo!=0||pcell->numvert!=0)
		pcell->inoutattrib=0;
}

PointInPolyhedron::~PointInPolyhedron(){

	delete []  vertcoord;
	delete [] trips;
	delete [] vertattrib;
	delete [] triofnode;
	delete [] tneighb;
//	delete [] triused;
	polytree->freeSubTree(polytree->getRoot());
}


//void PointInPolyhedron::getCellSeqWithUnknownAttribFromaCell(CellNode *cnode,vector<CellNode *>  * &pcellseq,
					//									CellNode * &pcellm,int &ia,double pm[3]){

//	CellNode *pcellt,*pcell;

//	if(cnode==0) return;
//	pcellseq=new vector<CellNode *>;
//	pcellseq->push_back(cnode);
//	pcellt=cnode;
//	for(;;){
//		pm[0]=pcellt->bound[0]; pm[1]=pcellt->bound[1];pm[2]=pcellt->bound[2];
////		pcell=getTheNeighbOfCellAtSpeciDirectWithRefPoint(pcellt,-1,0,pm);
	//	if(pcell==0){
	//		ia=-1; pcellm=0;
	//		return;
	//	}else if(pcell->inoutattrib!=-2){
	//		ia=pcell->inoutattrib; pcellm=pcell;
	//		return;
	//	}
	//	pcellseq->push_back(pcell);
	//	pcellt=pcell;
	//}
//}

int PointInPolyhedron::classifyVert(double p[3],int vert){

	int numbv, *neighbverts;
	int vertridge;
	double maxcosa;
	int tria,trib;

	getVertsAroundaVert(vert,neighbverts,numbv);
	getThePointFormingLeastAngleWith2Points(p,vert,neighbverts,numbv,maxcosa,vertridge);
	delete neighbverts;
	if(maxcosa>epscoplanar)
	//{
	//	get2TriCom2Vert(vert,vertridge,tria,trib);
	//	double d=sqDistPointToTri(p,vertcoord[trips[tria][0]],vertcoord[trips[tria][1]],vertcoord[trips[tria][2]]);
	//	d=sqDistPointToSeg3D(p,vertcoord[trips[trib][2]],vertcoord[trips[trib][0]]);
	//}
		jf_error("classify");
	get2TriCom2Vert(vert,vertridge,tria,trib);
	int tri0=tria;
	do{
		int rt=classifyEdge(tria,trib);
		if(rt==-1||rt==1){
			vertattrib[vert]=rt;
			return rt;
		}
		tria=trib;
		trib=nextTriOfVert(vert,trib);
	}while(tria!=tri0);
	vertattrib[vert]=0;
	return 0;
}
int PointInPolyhedron::classifyEdge(int tria,int trib){

	int ind=indexOfNeighbTriToTri(tria,trib);
	int vt=trips[tria][ind];
	double dt=VolumOf4p(vertcoord[trips[trib][0]],
					vertcoord[trips[trib][1]],vertcoord[trips[trib][2]],vertcoord[vt]);
	if(fabs(dt)<=epscoplanar)
		return 0;
	else if(dt<0) return -1;
	else return 1;
}
double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]){

	double v0p[3],v20[3],v01[3];
	vec_2p(p0,p,v0p);
	vec_2p(p2,p0,v20);
	vec_2p(p0,p1,v01);

	double d0p20=vec_dotp(v0p,v20);
	double d0p01=vec_dotp(v0p,v01);
	if(d0p20>=0&&d0p01<=0) return SqDistance3D(p,p0);
	
	double v1p[3],v12[3];
	vec_2p(p1,p,v1p);
	vec_2p(p1,p2,v12);
	double d1p01=vec_dotp(v1p,v01);
	double d1p12=vec_dotp(v1p,v12);
	if(d1p01>=0&&d1p12<=0) return SqDistance3D(p,p1);

	double v2p[3];
	vec_2p(p2,p,v2p);
	double d2p12=vec_dotp(v2p,v12);
	double d2p20=vec_dotp(v2p,v20);
	if(d2p12>=0&&d2p20<=0) return SqDistance3D(p,p2);

	double nm012[3],nm01p[3],nm12p[3],nm20p[3];
	vec_crop(v20,v01,nm012);

	vec_crop(v01,v0p,nm01p);
	double dt01=vec_dotp(nm012,nm01p);
	if(dt01<=0&&d0p01>=0&&d1p01<=0)
		return sqDistPointToSeg3D(p,p0,p1); // rt==0;
		
	vec_crop(v12,v1p,nm12p);
	double dt12=vec_dotp(nm012,nm12p);
	if(dt12<=0&&d1p12>=0&&d2p12<=0)
		return sqDistPointToSeg3D(p,p1,p2);
	
	vec_crop(v20,v2p,nm20p);
	double dt20=vec_dotp(nm012,nm20p);
	if(dt20<=0&&d2p20>=0&&d0p20<=0)
		return sqDistPointToSeg3D(p,p2,p0);

//if(dt01>0&&dt12>0&&dt20>0){
  if(dt01>=0&&dt12>=0&&dt20>=0){
		double a=vec_dotp(nm012,v0p);
		return a*a/vec_sqval(nm012);
	}else
		jf_error("asdf posiotin");
}


void PointInPolyhedron::getVertsAroundaVert(int v, int * &nbverts,int &numnbv){

	int ctri,tri0,count;

	ctri=tri0=triofnode[v];
	count=0;
	do{
		count++;
		ctri=nextTriOfVert(v,ctri);
	}while(ctri!=tri0);
	if(count<=2)
		jf_error("err getvertsarounda");
	nbverts=(int *) new int[count];
	numnbv=count;
	count=0;
	do{
		int vt=nextVertOfTri(ctri,v);
		nbverts[count++]=vt;
		ctri=nextTriOfVert(v,ctri);
	}while(ctri!=tri0);
}
int PointInPolyhedron::indexOfVertAtTri(int v, int ctri){

	if(trips[ctri][0]==v) return 0;
	else if(trips[ctri][1]==v) return 1;
	else if(trips[ctri][2]==v) return 2;
	else jf_error("indexoftri #2\n");
}
int PointInPolyhedron::indexOfNeighbTriToTri(int tria,int trinb){

	if(tneighb[tria][0]==trinb) return 0;
	else if(tneighb[tria][1]==trinb) return 1;
	else if(tneighb[tria][2]==trinb) return 2;
	else jf_error("indexofneighb");
}
int PointInPolyhedron::nextTriOfVert(int v, int ctri){
	int ind=indexOfVertAtTri(v,ctri);
	return tneighb[ctri][(1+ind)%3];
}
int PointInPolyhedron::neighbOfTri(int tri,int ind){
	return tneighb[tri][ind];
}
int PointInPolyhedron::nextVertOfTri(int tri,int v){
	if(v==trips[tri][0]) return trips[tri][1];
	else if(v==trips[tri][1]) return trips[tri][2];
	else if(v==trips[tri][2]) return trips[tri][0];
	else jf_error("nextvoftri");
}
void PointInPolyhedron::getThePointFormingLeastAngleWith2Points(double p[3],int v, int *nbverts,int numnbv,double &maxcosa,int &vridge){

	double pv[3],pvp[3],pvpi[3],dp;
	maxcosa=-1.;
	copy3DPoint(vertcoord[v],pv);
	vec_2p(pv,p,pvp);
	vec_uni(pvp);
	for(int i=0; i<numnbv; i++){
		vec_2p(pv,vertcoord[nbverts[i]],pvpi);
		vec_uni(pvpi);
		if((dp=vec_dotp(pvpi,pvp))>maxcosa){
			if(dp>epscoplanar)
				dp=dp;
			maxcosa=dp;
			vridge=nbverts[i];
		}
	}
	
}
void PointInPolyhedron::get2TriCom2Vert(int va, int vb, int &ta, int &tb){

	int tri0;
	tri0=ta=triofnode[va];
	do{
		tb=nextTriOfVert(va,ta);
		if(nextVertOfTri(tb,va)==vb)
			return;
		ta=tb;
	}while(ta!=tri0);
	jf_error("get2triwith");
}
double sqDistPointToSeg3D(double p[3],double p0[3],double p1[3]){

	double vp0p[3],vp0p1[3],vp1p[3];

	vec_2p(p0,p,vp0p);
	vec_2p(p0,p1,vp0p1);
	if(vec_dotp(vp0p,vp0p1)<=0)
		return SqDistance3D(p0,p);
	vec_2p(p1,p,vp1p);
	double prjp1p=vec_dotp(vp1p,vp0p1);
	double sqdp1p=SqDistance3D(p1,p);
	if(prjp1p>=0)
		return sqdp1p;
	double sqd=sqdp1p-prjp1p*prjp1p/vec_sqval(vp0p1);
	if(sqd<0){
		cout<<sqd<<" less than 0"<<endl;
		sqd=0;
	}
	return sqd;  //?
}


void PointInPolyhedron::formNeighbAndTriOfNode(){

	int i;
	int *numtriofnode=new int[numvert];
	int *tripositionofnode=new int[numvert];
	for(i=0; i<numvert; i++)  //record numbers of triangles around each node
		numtriofnode[i]=0;
	for(i=0; i<numtri; i++){
		for(int j=0; j<3; j++)
			numtriofnode[trips[i][j]]++;
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvert; i++)  //positions of first triangles in the trilist for each node
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	int *trilist=(int *) new int[3*numtri];
	for( i=0; i<numtri; i++){
		for(int j=0; j<3; j++){
			trilist[tripositionofnode[ trips[i][j] ]]=i;
			tripositionofnode[ trips[i][j] ]++;
		}
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvert; i++) //recover the first positions
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	for(int i=0; i<numvert; i++)
		triofnode[i]=trilist[tripositionofnode[i]];
	recNeighbOfTrips(numtriofnode,tripositionofnode,trilist);
	delete [] numtriofnode;
	delete [] tripositionofnode;
	delete [] trilist;
}
void PointInPolyhedron::recNeighbOfTrips(int *numtriofnode,int *tripositionofnode,int *trilist){

	int tnb,nbindex;

	for(int i=0; i<numtri; i++)
		for(int j=0; j<3; j++)
			tneighb[i][j]=-1;
	for(int i=0; i<numtri; i++){
		for(int j=0; j<3; j++){
			if(tneighb[i][j]!=-1) continue;
			getNeighbFromTrilist(i,j,tnb,nbindex,numtriofnode,tripositionofnode,trilist);
			tneighb[i][j]=tnb;
			tneighb[tnb][nbindex]=i;
		}
	}
}
void  PointInPolyhedron::getNeighbFromTrilist(int tri,int ind,int &tnb,int &nbindex
								,int *numtriofnode,int *tripositionofnode,int *trilist){

	int a,b;
	getEdgeOfTri(trips[tri],ind,a,b);
	for(int i=0; i<numtriofnode[a]; i++){
		int ip=tripositionofnode[a]+i;
		int ctri=trilist[ip];
		if(ctri==tri) continue;
		if((trips[ctri][0]==a&&trips[ctri][1]==b)||(trips[ctri][1]==a&&trips[ctri][0]==b)){
			tnb=ctri;
			nbindex=2;
			return;
		}else if((trips[ctri][1]==a&&trips[ctri][2]==b)||(trips[ctri][2]==a&&trips[ctri][1]==b)){
			tnb=ctri;
			nbindex=0;
			return;
		}else if((trips[ctri][2]==a&&trips[ctri][0]==b)||(trips[ctri][0]==a&&trips[ctri][2]==b)){
			tnb=ctri;
			nbindex=1;
			return;
		}
	}
	jf_error("err getneighfromtrl");
}

void  PointInPolyhedron::getEdgeOfTri(int np[3], int index, int &a, int &b){
	if(index==0){
		a=np[1]; b=np[2];
	}else	if(index==1){
		a=np[2]; b=np[0];
	}else	if(index==2){
		a=np[0]; b=np[1];
	}else
		jf_error("error getedgeoftri");
}

void PIP_jianfei_cpp(double *vertices, int *numV,
                     int    *faces,    int *numF,
                     double *query,    int *numQ,
                     int    *result) {
//Rprintf("PIP_jianfei_cpp: entered C++ function.\n");
	// Transfer vertex data from flat double array to double[3] array.
	double (*vert)[3];
	vert = ( double (*)[3]) new double [3*(*numV)];
#if 0
	double minX = DLB_MAX, minY=DBL_MAX, minZ=DBL_MAX;
#else
	double minX = FLT_MAX, minY=FLT_MAX, minZ=FLT_MAX;
#endif
	int i;
    for ( i = 0; i < (*numV); i++ ) {
        vert[i][0] = vertices[i+0*(*numV)];
        vert[i][1] = vertices[i+1*(*numV)];
        vert[i][2] = vertices[i+2*(*numV)];
 
        if ( minX > vert[i][0] ) minX = vert[i][0];
        if ( minY > vert[i][1] ) minY = vert[i][1];
        if ( minZ > vert[i][2] ) minZ = vert[i][2];
    }
 
    // Set minimum coordinates to 0.
    for ( i = 0; i < (*numV); i++ ) {
        vert[i][0] -= minX;
        vert[i][1] -= minY;
        vert[i][2] -= minZ;
    }
 
 	// Transfer face data from flat int array to int[3] array.
 	// Decrement vertex indices by 1, since Jianfei's code starts indexing at 0.
	int (*tris)[3];
    tris = (int (*)[3]) new int [3*(*numF)];
    for( i=0; i<(*numF); i++) {
        tris[i][0] = faces[i+0*(*numF)] - 1;
        tris[i][1] = faces[i+1*(*numF)] - 1;
        tris[i][2] = faces[i+2*(*numF)] - 1;
    }

    // Attempt to construct a object of PointInPolyhedron
	PointInPolyhedron *ptpoly = 0;

    try {
 	    ptpoly = new PointInPolyhedron(vert,(*numV),tris,(*numF));
    }
    catch ( int ptpolyError ) {
        // Fill result vector with the code "-2" to indicate
        // a failed initialization.
        for( i=0; i<(*numQ); i++) {
            result[i] = -2;
        }

        // Revert XYZ coordinates back to original values.
        for ( i = 0; i < (*numV); i++ ) {
            vert[i][0] += minX;
            vert[i][1] += minY;
            vert[i][2] += minZ;
        }
        return;
    }

    // Loop over queries, feed them to the Jianfei method.
    // Don't forget about the minX, minY, and minZ shifts.
    double q[3]={0,0,0};
    for( i=0; i<(*numQ); i++) {
        q[0] = query[i+0*(*numQ)] - minX;
        q[1] = query[i+1*(*numQ)] - minY;
        q[2] = query[i+2*(*numQ)] - minZ;
        result[i] = ptpoly->isPinPolyhedron(q);
    }

    // RELEASE MEMORY!!
    delete [] tris;
    delete [] vert;
    delete ptpoly;
}

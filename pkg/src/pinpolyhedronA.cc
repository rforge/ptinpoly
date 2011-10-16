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

bool   PointInPolyhedron::ifexinfoshouldbeincell(void *info,int infotype,CellNode3D *cnode){
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
	CellNode3D *pcell;
//	vector<CellNode3D *> *pcellseq;

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
	CellNode3D *pcellm=0;
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

int PointInPolyhedron::testPinPolyhedronForPinGcell(double p[3],CellNode3D *cnode){

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


void PointInPolyhedron::getRelativeClosestEntityForPointInGCell( double p[3],CellNode3D *cnode,int &id,
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

void PointInPolyhedron::getAbsoluteClosestTriForPointInGCell(double p[3],CellNode3D *cnode, int &tri, double &dist){

	CellNode3D *pcell0=0;
	CellNode3D *pcell=cnode;
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
										   CellNode3D *pcell,int &tri,double &dist){
	
	double distn;
	int trin;

	dist=dist0, tri=-1;
	if(sqdistPointToBox(p,pcell->bound)>=dist0) return;
	if(pcell->isLeaf()){
		getTheClosestTriAmongCell(p,pcell,distn,trin);
		if(distn<dist){	dist=distn;	tri=trin; return;}
	}else{
		CellNode3D *sortsub[2]={pcell->child[0],pcell->child[1]};
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
/*void PointInPolyhedron::recoverTriused(CellNode3D *pcell){
	
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
void PointInPolyhedron::getTheClosestTriAmongCell(double p[3],CellNode3D *pcell, double &dist,int &ntri){

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

void PointInPolyhedron::setGCellAttribOfSubTree(CellNode3D *pcell){

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


//void PointInPolyhedron::getCellSeqWithUnknownAttribFromaCell(CellNode3D *cnode,vector<CellNode3D *>  * &pcellseq,
					//									CellNode3D * &pcellm,int &ia,double pm[3]){

//	CellNode3D *pcellt,*pcell;

//	if(cnode==0) return;
//	pcellseq=new vector<CellNode3D *>;
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

// Three-Dimensional code.
void PIP3D_jianfei_cpp(double *vertices, int *numV,
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
            result[i] = -ptpolyError;
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
        q[0]      = query[i+0*(*numQ)] - minX;
        q[1]      = query[i+1*(*numQ)] - minY;
        q[2]      = query[i+2*(*numQ)] - minZ;
        result[i] = ptpoly->isPinPolyhedron(q);
    }

    // RELEASE MEMORY!!
    delete [] tris;
    delete [] vert;
    delete ptpoly;
}

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
// Two-Dimensional code.
//

void vec_2p2(double *pointa,double *pointb, double *vector);
double vec_dotp2(double *vector1,double *vector2);
double squareDist2p2( double *v1,double *v2);
double vec_crop2(double *vector1,double *vector2);
double squareDistPointToLine2( double p[2],double ps[2], double pe[2]);
void copy2DPoint(double pfr[2],double pto[2]);
//void jf_error(char *ch); Already declared above with 3D code
bool if3PointRightHandSort(double p0[2],double p1[2],double p2[2]);
bool ifBoxContainPoint(double p[2],double bound[4]);
//bool ifBoxWith3OpenEdgeContainPoint(double p[2],double bound[4],int ie);
bool if2BoxNeighb(double a[4],double b[4]);
int positionOfPointProjectToSeg(double p[2],double p0[2],double p1[2]);
bool ifPointOverlapWithBox2D(double p[2],double bd[4],double eps);
int convexityOf3Point(double p0[2],double p1[2],double p2[2],double eps);
bool ifSegOverlapBox2D(double ps[2],double pe[2],double bd[4],double eps);
bool if4CornerOfBoxAtDifferentSideOfSeg(double ps[2],double pe[2],double bound[4]);
void boxOfPoints( double (*p)[2] , int num ,double box[4] );
double squareDistPointToSeg(double p[2], double p0[2],double p1[2]);
void getTheClosestPointAtSeg(double p[2], double p0[2],double p1[2],double pcha[2]);
double sqdistPointToBox(double p[2],double bd[4]);

class PolyQuadtree;
struct CellNode2D;
#include <vector>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <time.h>
using namespace std;


//inline double min ( double a, double b){
//	return b<a? b : a;
//}
//
//template <class T>
//inline const T& max ( const T& a, const T& b){
//	return a<b? b : a;
//}

class PolyQuadtree{
public:
	int isPinpolygon(double p[2]); //-1,0,1
	PolyQuadtree(double (*vti)[2], int numvi);
	PolyQuadtree(double (*vti)[2], int numvi,int (*seg2ei)[2],int numsi);

	double absoluteClosestSqDistance(double p[2]){
		CellNode2D *cnode=findaLeafCellContainingPoint(root,p);
		int nseg; double dist;
		getAbsoluteClosestSegForPointInGCell(p,cnode,nseg,dist);
		return dist;
	}
	bool isPointInGreyCell(double p[2]);	
	~PolyQuadtree();
private:
	int testPinpolygonForPinGcell(double p[2],CellNode2D *cnode);
	void getRelativeClosestEntityForPointInGCell( double p[2],CellNode2D *cnode,int &id ,int &nentity,
		double &dist);
	void getRelativeClosestSegForPointInGCell(double p[2], CellNode2D *cnode,int &nseg, double &dist);
	void findaCloserConvergentCharacterPoint(double p[2], CellNode2D *pcell0,double dist0,
									 CellNode2D *&pncell,double &dist,int &nseg,double pcha[2]);
	void findTheCharacterPointOfSeg(int nseg,double p[2],CellNode2D *pcell0,double pcha[2],CellNode2D * &pncell);
	void getTheClosestSegAmongCell(double p[2],CellNode2D *pcell, double &dist,int &nseg);
	void getEndPointOfSeg(int seg, double p0[2],double p1[2]);
	CellNode2D *findaLeafCellContainingPoint(CellNode2D *pcell,double p[2]);
	bool if2CellNeighb(CellNode2D *pcell0, CellNode2D *pcell1);
//	CellNode2D *findTheLeafCellWithPointAtSpecifiedEdge(CellNode2D *pcell, double p[2], int ie);
	CellNode2D *getTheNeighbOfCellAtSpeciDirectWithRefPoint(CellNode2D *pcell,int ix,int iy,
													  double pmf[2],double pmto[2]);
	CellNode2D *getaNeighbCellCloserToAnotherCell(CellNode2D *pscell,CellNode2D *pecell,double pmf[2],
		double pt[2]);
	void creV2seg(void);
	void insertVertInSubTree(int v, CellNode2D *cnode);
	void splitNode(CellNode2D *cnode);
	void insertSegInSubTree(int seg,CellNode2D *cnode);
	void compVertattrib(void);
	void getCellSeqWithUnknownAttribFromaCell(CellNode2D *cnode,vector<CellNode2D *> *&pcellseq, 
		CellNode2D * &pcellm,int &ia,double pm[2]);
	void freeSubQuadtree(CellNode2D *pcell);
	void setGCellAttribOfSubTree(CellNode2D *pcell);
	CellNode2D *findTheNearestAncestorContainingPoint(CellNode2D *pcell,double pcha[2]);
	void sortTheDistancesOfChildrenFromPoint(double p[2],CellNode2D *pcell,CellNode2D *sortsub[4]);
	void getTheClosestSegNonLeaf(double p[2],CellNode2D *pcell0,double dist0,
										   CellNode2D *pcell,int &seg,double &dist);
	bool getTheRelativeClosestSegNonLeaf(double p[2],CellNode2D *pcell0,//double dist0,
										   CellNode2D *pcell,int &seg,double &dist,double pcha[2]);
	void getAbsoluteClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist);
	CellNode2D* getNextCell(CellNode2D *cnode,double ps[2],double pe[2]);
//	CellNode2D *findTheNearestAncestorWithPointAtEdge(CellNode2D *pcell,double pt[2],int ie);

private:
	static const double epsilonon;
	static const double epsoverlap;
	static const double epscoplanar;
	double epscell;
	CellNode2D *root;
	double (*vert)[2];
	int numvert;
	int (*v2seg)[2];
	int *vertattrib;
	int (*seg2end)[2],numseg;
};
struct CellNode2D{
	vector<int> *psegar;
	int vertincell;
	int inoutattrib; //-1,0,1
	double bound[4];
	CellNode2D *child[4];
	CellNode2D *parent;
	
	CellNode2D(double bd[4]);
	~CellNode2D();
	bool isEmpty(){ return psegar==0&&vertincell==-1;}
	bool isLeaf(){ return !child[0];}
};
const double PolyQuadtree::epsilonon=0.00000000000001;
const double PolyQuadtree::epsoverlap=0.000001;
const double PolyQuadtree::epscoplanar=0.000001;

// int absolute; Already defined above.

int PolyQuadtree::isPinpolygon( double p[2]){

	int rt;
	CellNode2D *pcell;
	vector<CellNode2D *> *pcellseq;

	pcell=findaLeafCellContainingPoint(root,p);
	if(pcell==0) 
		return -1;
	if(pcell->inoutattrib==-1||pcell->inoutattrib==1)
		return pcell->inoutattrib;
	else if(pcell->inoutattrib ==0)
		return testPinpolygonForPinGcell(p,pcell);
	//else if(pcell->inoutattrib==-2)
	double pm[2];
	CellNode2D *pcellm=0;
	getCellSeqWithUnknownAttribFromaCell(pcell,pcellseq,pcellm,rt,pm);
	if(rt==0)
		rt=testPinpolygonForPinGcell(pm,pcellm);
	if(rt==0)
		jf_error("ispinopolygon");
	if(pcellseq!=0)
		for(unsigned i=0; i<pcellseq->size(); i++)
			(*pcellseq)[i]->inoutattrib=rt;
	delete pcellseq;
	return rt;
}

int PolyQuadtree::testPinpolygonForPinGcell(double p[2],CellNode2D *cnode){

	int id,nentity;
	double dist,p0[2],p1[2];

	getRelativeClosestEntityForPointInGCell(p,cnode,id,nentity,dist);
	if(dist<=epsilonon)
		return 0;
	if(id==0){
		if(vertattrib[nentity]!=0)
			return vertattrib[nentity];
		else
			nentity=v2seg[nentity][0];
	}
	if(id!=1&&id!=0)
		jf_error("err ispoinpolygon");
	getEndPointOfSeg(nentity,p0,p1);
	if(if3PointRightHandSort(p,p0,p1))
		return 1;
	else 
		return -1;
}


void PolyQuadtree::getRelativeClosestEntityForPointInGCell( double p[2],CellNode2D *cnode,int &id,
														   int &nentity, double &dist){

	int nseg,ip;
	double p0[2],p1[2];

//	if(absolute==0)
		getRelativeClosestSegForPointInGCell(p,cnode,nseg,dist);
//	else
//		getAbsoluteClosestSegForPointInGCell(p,cnode,nseg,dist);
	if(dist==numeric_limits<double>::max())
		jf_error("err getrelativeclosetentityforpingcell");
	getEndPointOfSeg(nseg,p0,p1);
	if((ip=positionOfPointProjectToSeg(p,p0,p1))==0){
		nentity=nseg;
		id=1;
	}else if(ip==-1){
		nentity=seg2end[nseg][0];
		id=0;
	}else{
		nentity=seg2end[nseg][1];
		id=0;
	}
}
double sqdistPointToBox(double p[2],double bd[4]){

	double a=0,b=0;

	if(p[0]>bd[2]) a=p[0]-bd[2];
	else if(p[0]<bd[0]) a=bd[0]-p[0];
	else a=0;
	if(p[1]>bd[3]) b=p[1]-bd[3];
	else if(p[1]<bd[1]) b=bd[1]-p[1];
	else b=0;
	return a*a+b*b;
}
double sqdistInnerPointToBoxBound(double p[2],double bd[4]){

	double a=min(p[0]-bd[0],bd[2]-p[0]);
	double b=min(p[1]-bd[1],bd[3]-p[1]);
	double c= min(a,b);
	return c*c;
}
void PolyQuadtree::getAbsoluteClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist){

	CellNode2D *pcell0=0;
	CellNode2D *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosestseg");
	dist=numeric_limits<double>::max();
	seg=-1;
	
	while(pcell){
		int segn;
		double distn;
		getTheClosestSegNonLeaf(p,pcell0,dist,pcell,segn,distn);
		if(distn<dist){
			dist=distn; seg=segn;
		}
		if(sqdistInnerPointToBoxBound(p,pcell->bound)>=dist)
			return;
		pcell0=pcell;
		pcell=pcell->parent;
	}
}
void findOutPointofBox2D(double ps[2],double pe[2],double *pl,double *ph,double eps,double px[2]){

	for(int i=0; i<2; i++){
		if(ps[i]<=pe[i]){
			if(pe[i]<=ph[i]) px[i]=pe[i];
			else px[i]=ph[i]+eps;
		}else{
			if(pe[i]<pl[i]) px[i]=pl[i]-eps;
			else px[i]=pe[i];
		}
	}
}
CellNode2D* PolyQuadtree::getNextCell(CellNode2D *cnode,double ps[2],double pe[2]){

	double px[2];
	CellNode2D *pancestor,*pcell;

	findOutPointofBox2D(ps,pe,cnode->bound,cnode->bound+2,epscell,px);
	if((pancestor=findTheNearestAncestorContainingPoint(cnode,px))==0)
		return 0;
	else pcell=findaLeafCellContainingPoint(pancestor,px);
	if(pcell==cnode)
		jf_error("err epscell, contact the developer please,liujianfei@pku.edu.cn");
	return pcell;
}

void PolyQuadtree::getRelativeClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist){

//	CellNode2D *pcell0=0;
	CellNode2D *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosesttri");
	dist=numeric_limits<double>::max();
	seg=-1;
	double pcha[2]={dist,0},p0[2],p1[2];
	
	for(;;){
		int segn;
		double distn;
		bool update;
		update=false;
		getTheClosestSegAmongCell(p,pcell,distn,segn);
		if(distn<dist){
			dist=distn; seg=segn; update=true;
			getEndPointOfSeg(seg,p0,p1);
			getTheClosestPointAtSeg(p,p0,p1,pcha);
		}
		if(ifBoxContainPoint(pcha,pcell->bound)) break;
		if(update)
			pcell=getNextCell(cnode,p,pcha);
		else
			pcell=getNextCell(pcell,p,pcha);
		if(pcell==0) jf_error(" err getrelative");
	}
}
/*
void PolyQuadtree::getRelativeClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist){

	CellNode2D *pcell0=0;
	CellNode2D *pcell=cnode;
	dist=numeric_limits<double>::max();
	seg=-1;
	double pcha[2]={numeric_limits<double>::max(),0};
	while(pcell){
		int segn;
		double distn;
		bool rt=getTheRelativeClosestSegNonLeaf(p,pcell0,pcell,seg,dist,pcha);//dist,pcell,segn,distn,pcha);
		//getTheClosestSegNonLeaf(p,pcell0,dist,pcell,segn,distn);
	//	if(distn<dist){
	//		dist=distn; seg=segn;
		//	double p0[2],p1[2];
		//	getEndPointOfSeg(seg,p0,p1);
		//	getTheClosestPointAtSeg(p,p0,p1,pcha);
	//	}
		if(pcha[0]==numeric_limits<double>::max()) jf_error("err getrelativeclosestseg");
//		if(!rt&&(sqdistInnerPointToBoxBound(p,pcell->bound)>=dist||
//		if(ifBoxContainPoint(pcha,pcell->bound)) return;
//			jf_error("ewr");
		if(	rt)	return;
		pcell0=pcell;
		pcell=pcell->parent;
	}
	jf_error("jsald;fjlsakdf;ajsdf;");
}

*/
void PolyQuadtree::sortTheDistancesOfChildrenFromPoint(double p[2],CellNode2D *pcell,CellNode2D *sortsub[4]){

	double distpc[4];

	if(pcell==0||pcell->isLeaf())
		jf_error("sortthedist");
	for(int i=0; i<4; i++)
		distpc[i]=sqdistPointToBox(p,pcell->//bound);
											child[i]->bound);
	for(int i=0; i<4; i++){
		int index=0;
		for(int j=0; j<4; j++)
			if(distpc[i]>distpc[j]||i>j&&distpc[i]==distpc[j])
				index++;
		sortsub[index]=pcell->child[i];
	}
}
void PolyQuadtree::getTheClosestSegNonLeaf(double p[2],CellNode2D *pcell0,double dist0,
										   CellNode2D *pcell,int &seg,double &dist){
	
	double distn;
	int segn;
	dist=dist0, seg=-1;
	if(pcell==pcell0||sqdistPointToBox(p,pcell->bound)>=dist0) return;
	if(pcell->isLeaf()){
		getTheClosestSegAmongCell(p,pcell,distn,segn);
		if(distn<dist){	dist=distn;	seg=segn; return;}
	}else{
		CellNode2D *sortsub[4];
		sortTheDistancesOfChildrenFromPoint(p,pcell,sortsub);
		for(int i=0; i<4; i++){
			getTheClosestSegNonLeaf(p,pcell0,dist,sortsub[i],segn,distn);
			if(distn<dist){ dist=distn; seg=segn; }
		}
	}
}

bool PolyQuadtree::getTheRelativeClosestSegNonLeaf(double p[2],CellNode2D *pcell0,//double dist0,
										   CellNode2D *pcell,int &seg,double &dist,double pcha[2]){
	
	double distn;
	int segn;

//	double pcha[2];  

	//dist=dist0, seg=-1;
	if(pcell==pcell0||sqdistPointToBox(p,pcell->bound)>=dist) return false;
	if(pcell->isLeaf()){
		getTheClosestSegAmongCell(p,pcell,distn,segn);
		//if(distn<dist){	dist=distn;	seg=segn; return;}
			if(distn<dist){ 
				dist=distn; seg=segn; 
				double p0[2],p1[2];
				getEndPointOfSeg(seg,p0,p1);
				getTheClosestPointAtSeg(p,p0,p1,pcha);
			}
		if(pcha[0]==numeric_limits<double>::max()) jf_error("err getrelativeclosestseg");
		if(ifBoxContainPoint(pcha,pcell->bound)) //sortsub[i]->bound))
			return true;
	}else{
		CellNode2D *sortsub[4];
		sortTheDistancesOfChildrenFromPoint(p,pcell,sortsub);
		for(int i=0; i<4; i++){
//			if(sortsub[i]==pcell0) continue;
			bool rt=getTheRelativeClosestSegNonLeaf(p,pcell0,sortsub[i],seg,dist,pcha);//segn,distn,pcha) ;
//			if(distn<dist){ dist=distn; seg=segn;} 
			if(rt) return true;
		}
	}
	if(ifBoxContainPoint(pcha,pcell->bound)||pcha[0]==numeric_limits<double>::max()) //sortsub[i]->bound))
		jf_error("jkasldf");
	//	return true;
	return false;
}
/*
void PolyQuadtree::getRelativeClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist){

	CellNode2D *pcellp,*pcell,*pncell;
	double distn;
	double pcha[2],pmto[2],pmove[2];
	int nseg;

	pcellp=pcell=cnode;
	if(!pcell||pcell->isEmpty())
		jf_error("getrelativeclosestseg");
	dist=numeric_limits<double>::max();
	seg=-1;
	for(;;){
		findaCloserConvergentCharacterPoint(p,pcell,dist,pncell,distn,nseg,pcha);
		if(distn<dist){
			seg=nseg; dist=distn;
			pcell=pncell;
			copy2DPoint(pcha,pmove);
		}
		if((pcellp)==(pcell)||if2CellNeighb(pcellp,pcell))
			return;
		pcell=getaNeighbCellCloserToAnotherCell(pcell,pcellp,pmove,pmto);
		copy2DPoint(pmto,pmove);
	}
}*//*
void PolyQuadtree::getRelativeClosestSegForPointInGCell(double p[2],CellNode2D *cnode, int &seg, double &dist){

	CellNode2D *pcellp,*pscell,*pecell=0;
	double distn;
	double pcha[2],pmto[2],pmove[2];
	int nseg;

	pcellp=pscell=cnode;
	if(!pscell||pscell->isEmpty())
		jf_error("getrelativeclosestseg");
	dist=numeric_limits<double>::max();
	seg=-1;
	for(;;){
		getTheClosestSegAmongCell(p,pscell,distn,nseg);
		if(distn<dist){
			seg=nseg; dist=distn;
			findTheCharacterPointOfSeg(nseg,p,pscell,pcha,pecell);
			if(pscell==pecell) return;
			else {pscell=pcellp;copy2DPoint(p,pmove);}
		//	copy2DPoint(pcha,pmto);
		}
		if(pecell==0) jf_error("err getrelativecloseseg");
		if(pscell==pecell) return;
		if(if2CellNeighb(pscell,pecell)) pscell=pecell;
		else{
			pscell=getaNeighbCellCloserToAnotherCell(pscell,pecell,pmove,pmto);
			copy2DPoint(pmto,pmove);
		}
	}
}*/
CellNode2D * PolyQuadtree::findTheNearestAncestorContainingPoint(CellNode2D *pcell0,double pcha[2]){
	
	CellNode2D *pcell=pcell0;
	for(;;){
		if(pcell==0) return 0;
		if(ifBoxContainPoint(pcha,pcell->bound)) return pcell;
		else pcell=pcell->parent;
	}
}
void PolyQuadtree::findTheCharacterPointOfSeg(int nseg,double p[2],CellNode2D *pcell0,double pcha[2],CellNode2D * &pncell){

	double p0[2],p1[2];

	getEndPointOfSeg(nseg,p0,p1);
	getTheClosestPointAtSeg(p,p0,p1,pcha);
	CellNode2D *pancestor=findTheNearestAncestorContainingPoint(pcell0,pcha);
	pncell=findaLeafCellContainingPoint(pancestor,pcha);
	if(pncell==0) jf_error("err findaclosercp");
}
void PolyQuadtree::findaCloserConvergentCharacterPoint(double p[2], CellNode2D *pcell0,double dist0,
									 CellNode2D * &pncell,double &distn,int &nseg,double pcha[2]){

	CellNode2D *pcell,*pancestor;
	double p0[2],p1[2];

	getTheClosestSegAmongCell(p,pcell0,distn,nseg);
	if(distn>=dist0)
		return;
	pcell=pcell0;
	for(;;){
		getEndPointOfSeg(nseg,p0,p1);
		getTheClosestPointAtSeg(p,p0,p1,pcha);
		pancestor=findTheNearestAncestorContainingPoint(pcell,pcha);
		pncell=findaLeafCellContainingPoint(pancestor,pcha);
		if(pncell==0) jf_error("err findaclosercp");
		if(pncell==pcell) return; 
		else
			pcell=pncell;
		getTheClosestSegAmongCell(p,pcell,distn,nseg);
	}
}

void PolyQuadtree::getTheClosestSegAmongCell(double p[2],CellNode2D *pcell, double &dist,int &nseg){

	int seg;
	double distemp,p0[2],p1[2];

	dist=numeric_limits<double>::max();
	if(!pcell||!pcell->isLeaf())
		jf_error("error gettheclosetsegamongcell");
	if(pcell->psegar!=0)
		for(int i=0; i<pcell->psegar->size(); i++){
			seg=(*(pcell->psegar))[i];
			getEndPointOfSeg(seg,p0,p1);
			if((distemp=squareDistPointToSeg(p,p0,p1))<dist){
				dist=distemp;
				nseg=seg;
			}
		}
	if(pcell->vertincell==-1)
		return;
	seg=v2seg[pcell->vertincell][0];
	getEndPointOfSeg(seg,p0,p1);
	if((distemp=squareDistPointToSeg(p,p0,p1))<dist){ 
		dist=distemp;
		nseg=seg;
	}
	seg=v2seg[pcell->vertincell][1];
	getEndPointOfSeg(seg,p0,p1);
	if((distemp=squareDistPointToSeg(p,p0,p1))<dist){
		dist=distemp;
		nseg=seg;
	}
}

double squareDistPointToSeg(double p[2],double p0[2],double p1[2]){

	double vp0p[2],vp0p1[2],vp1p[2];

	vec_2p2(p0,p,vp0p);
	vec_2p2(p0,p1,vp0p1);
	if(vec_dotp2(vp0p,vp0p1)<=0)
		return squareDist2p2(p0,p);
	vec_2p2(p1,p,vp1p);
	if(vec_dotp2(vp1p,vp0p1)>=0)
		return squareDist2p2(p1,p);
	return squareDistPointToLine2(p,p0,p1);
}

void PolyQuadtree::getEndPointOfSeg(int seg, double p0[2],double p1[2]){

	p0[0]=vert[seg2end[seg][0]][0];
	p0[1]=vert[seg2end[seg][0]][1];
	p1[0]=vert[seg2end[seg][1]][0];
	p1[1]=vert[seg2end[seg][1]][1];
}

int positionOfPointProjectToSeg(double p[2],double p0[2],double p1[2]){

	double vp0p[2],vp0p1[2],vp1p[2];

	vec_2p2(p0,p,vp0p);
	vec_2p2(p0,p1,vp0p1);
	if(vec_dotp2(vp0p,vp0p1)<=0)
		return -1;
	vec_2p2(p1,p,vp1p);
	if(vec_dotp2(vp1p,vp0p1)>=0)
		return 1;
	return 0;
}

void getTheClosestPointAtSeg(double p[2],double p0[2],double p1[2],double pcha[2]){

	double ratio,dp1,dp2,vp0p[2],vp0p1[2],vp1p[2];

	vec_2p2(p0,p,vp0p);
	vec_2p2(p0,p1,vp0p1);
	if((dp1=vec_dotp2(vp0p,vp0p1))<=0){
		pcha[0]=p0[0];
		pcha[1]=p0[1];
		return;
	}
	vec_2p2(p1,p,vp1p);
	if((dp2=vec_dotp2(vp1p,vp0p1))>=0){
		pcha[0]=p1[0];
		pcha[1]=p1[1];
		return;
	}
	ratio=dp1/(dp1-dp2);
	pcha[0]=p0[0]+ratio*vp0p1[0];
	pcha[1]=p0[1]+ratio*vp0p1[1];
}

CellNode2D * PolyQuadtree:: findaLeafCellContainingPoint(CellNode2D *pcell,double p[2]){

	CellNode2D *rtpcell;

	if(!pcell||!ifBoxContainPoint(p,pcell->bound))
		return 0;
	if(pcell->isLeaf())
		return pcell;
	for( int i=0; i<4; i++)
		if((rtpcell=findaLeafCellContainingPoint(pcell->child[i],p))!=0)
			return rtpcell;
	jf_error("err findaleafcellcontainp");
}

bool PolyQuadtree:: if2CellNeighb(CellNode2D *pcell0, CellNode2D *pcell1){

	if(!pcell0||!pcell1) 
		jf_error("err is2cellneigh");
	if(if2BoxNeighb(pcell0->bound,pcell1->bound))
		return true;
	else
		return false;
}

bool if2BoxNeighb(double a[4],double b[4]){

	if(a[0]>b[2]||a[1]>b[3]||a[2]<b[0]||a[3]<b[1])
		return false;
	return true;
}

CellNode2D * PolyQuadtree::getaNeighbCellCloserToAnotherCell(CellNode2D *pscell,
									CellNode2D *pecell,double pmf[2],double pt[2]){

	int idx=0,idy=0;

	if(!pscell||!pecell)
		jf_error("err getneighbcellcloser");
//	pt[0]=pmf[0]; pt[1]=pmf[1];
	if(pscell->bound[0]>pecell->bound[2])
		idx=-1;
	else if(pscell->bound[1]>pecell->bound[3])
		idy=-1;
	if(pscell->bound[2]<pecell->bound[0])
		idx=1;
	else if(pscell->bound[3]<pecell->bound[1])
		idy=1;
	if(idx==0&&idy==0)
		jf_error("err getaneigh");
	return getTheNeighbOfCellAtSpeciDirectWithRefPoint(pscell,idx,idy,pmf,pt);
}
CellNode2D * PolyQuadtree::getTheNeighbOfCellAtSpeciDirectWithRefPoint(CellNode2D *pcell,int ix,int iy,
													  double pmf[2],double pmto[2]=0){
	double pm[2];
	copy2DPoint(pmf,pm);
	if(ix==-1)
		pm[0]=pcell->bound[0];
	else if(ix==1)
		pm[0]=pcell->bound[2];
	if(iy==-1)
		pm[1]=pcell->bound[1];
	else if(iy==1)
		pm[1]=pcell->bound[3];
	if(pmto!=0)
		copy2DPoint(pm,pmto);
	if(ix==-1)
		pm[0]-=epscell;
	else if(ix==1)
		pm[0]+=epscell;
	if(iy==-1)
		pm[1]-=epscell;
	else if(iy==1)
		pm[1]+=epscell;
	CellNode2D *pancestor=findTheNearestAncestorContainingPoint(pcell,pm);
	if(pancestor==0)
		return 0;
	return findaLeafCellContainingPoint(pancestor,pm);
}


/*
CellNode2D * PolyQuadtree::findTheNearestAncestorWithPointAtEdge(CellNode2D *pcell0,double pt[2],int ie){

	CellNode2D *pcell=pcell0;
	for(;;){
		if(pcell==0) return 0;
		if(ifBoxWith3OpenEdgeContainPoint(pt,pcell->bound,ie)) return pcell;
		else pcell=pcell->parrent;
	}
}

CellNode2D *PolyQuadtree::findTheLeafCellWithPointAtSpecifiedEdge(CellNode2D *pcell, double p[2], int ie){

	CellNode2D *rtpcell;

	if(!pcell||!ifBoxWith3OpenEdgeContainPoint(p,pcell->bound,ie))
		return 0;
	if(pcell->isLeaf())
		return pcell;
	for( int i=0; i<4; i++)
		if((rtpcell=findTheLeafCellWithPointAtSpecifiedEdge(pcell->child[i],p,ie))!=0)
			return rtpcell;
	jf_error("err findaleafcellcontainp");
}

bool ifBoxWith3OpenEdgeContainPoint(double p[2],double bound[4],int ie){

	if(p[0]<bound[0]||p[1]<bound[1]||p[0]>bound[2]||p[1]>bound[3]) return false;
	if(p[0]>bound[0]&&p[1]>bound[1]&&p[0]<bound[2]&&p[1]<bound[3]) return true;
	if((ie==0&&p[1]==bound[1]&&p[0]!=bound[2])||(ie==1&&p[1]==bound[1]&&p[0]!=bound[0])||
	   (ie==2&&p[0]==bound[2]&&p[1]!=bound[3])||(ie==3&&p[0]==bound[2]&&p[1]!=bound[1])||
	   (ie==4&&p[1]==bound[3]&&p[0]!=bound[2])||(ie==5&&p[1]==bound[3]&&p[0]!=bound[0])||
	   (ie==6&&p[0]==bound[0]&&p[1]!=bound[3])||(ie==7&&p[0]==bound[0]&&p[1]!=bound[1]) ) 
	   return true;
	return false;
}

*/
void vec_2p2(double *pointa,double *pointb, double *vec)
{
  int i;
  for(i=0; i<2; i++)
		vec[i]=pointb[i]-pointa[i];
}

double vec_dotp2(double *vector1,double *vector2)
{
   return vector1[0]*vector2[0]+vector1[1]*vector2[1];
}

double squareDist2p2( double *v1,double *v2){

  return (v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]) ;
}

double vec_crop2(double *vector1,double *vector2)
{
  return vector1[0]*vector2[1]-vector2[0]*vector1[1] ;
}
double squareDistPointToLine2( double p[2],double ps[2], double pe[2] )
{
  	double pse[2],psp[2];

	vec_2p2(ps,pe,pse);
	double sqdse=pse[0]*pse[0]+pse[1]*pse[1];
	if(sqdse<=numeric_limits<double>::epsilon())
		jf_error("too short line found in squredistptol");
	vec_2p2(ps,p,psp);
	double d=vec_crop2(psp,pse);
	return d*d/sqdse;
}

void copy2DPoint(double pfr[2],double pto[2]){

	pto[0]=pfr[0];
	pto[1]=pfr[1];
}

/* Already defined in kodtree.cc
void jf_error(char *ch){

		printf("%s\n",ch);
		exit(1);
}
*/

bool if3PointRightHandSort(double p0[2],double p1[2],double p2[2]){

	double p01[2],p02[2];
	vec_2p2(p0,p1,p01);
	vec_2p2(p0,p2,p02);
	if(vec_crop2(p01,p02)>0) return true;
	else return false;
}

bool ifBoxContainPoint(double p[2],double bound[4]){

	if(p[0]>=bound[0]&&p[1]>=bound[1]&&p[0]<=bound[2]&&p[1]<=bound[3])
		return true;
	else
		return false;
}
PolyQuadtree::PolyQuadtree(double (*vti)[2], int numvi){

	numvert=numvi;
	numseg=numvi;
	vert=( double (*)[2]) new double [2*numvert]; ///??
	v2seg= (int (*)[2]) new int [2*numvert];
	vertattrib= new int [numvert];
	seg2end= (int (*)[2]) new int [2*numseg];
	for(int i=0; i<numvert; i++){
		vert[i][0]=vti[i][0];
		vert[i][1]=vti[i][1];
	}
	for(int i=0; i<numseg; i++){
		seg2end[i][0]=i;
		seg2end[i][1]=(i==numseg-1?0:i+1);
	}
	creV2seg();
	double bd[4];
	boxOfPoints(vert,numvert,bd);
	double lcube=max(bd[2]-bd[0],bd[3]-bd[1]);
	bd[2]=bd[0]+lcube,bd[3]=bd[1]+lcube;
	epscell=numeric_limits<double>::epsilon()*(1+max(bd[2]-bd[0],bd[3]-bd[1]));
	root=new CellNode2D(bd);
	for( int i=0; i<numvert; i++)
		insertVertInSubTree(i,root);
	for(int i=0; i<numseg; i++)
		insertSegInSubTree(i,root);
	compVertattrib();
	setGCellAttribOfSubTree(root);
}

PolyQuadtree::PolyQuadtree(double (*vti)[2], int numvi,int (*seg2ei)[2],int numsi){

	numvert=numvi;
	numseg=numsi;
	vert=( double (*)[2]) new double [2*numvert]; ///??
	v2seg= (int (*)[2]) new int [2*numvert];
	vertattrib= new int [numvert];
	seg2end= (int (*)[2]) new int [2*numseg];
	for(int i=0; i<numvert; i++){
		vert[i][0]=vti[i][0];
		vert[i][1]=vti[i][1];
	}
	for(int i=0; i<numseg; i++){
		seg2end[i][0]=seg2ei[i][0];
		seg2end[i][1]=seg2ei[i][1];
	}
	creV2seg();
	double bd[4];
	boxOfPoints(vert,numvert,bd);
	double lcube=max(bd[2]-bd[0],bd[3]-bd[1]);
	bd[2]=bd[0]+lcube,bd[3]=bd[1]+lcube;
	epscell=numeric_limits<double>::epsilon()*(1+max(bd[2]-bd[0],bd[3]-bd[1]));
	root=new CellNode2D(bd);
	for( int i=0; i<numvert; i++)
		insertVertInSubTree(i,root);
	for(int i=0; i<numseg; i++)
		insertSegInSubTree(i,root);
	compVertattrib();
	setGCellAttribOfSubTree(root);
}

void PolyQuadtree::setGCellAttribOfSubTree(CellNode2D *pcell){

	if(!pcell) return;
	if(!pcell->isLeaf())
		for(int i=0; i<4; i++)
			setGCellAttribOfSubTree(pcell->child[i]);
	else if(pcell->psegar!=0||pcell->vertincell!=-1)
		pcell->inoutattrib=0;
}

PolyQuadtree::~PolyQuadtree(){

	delete []  vert;
	delete [] v2seg;
	delete [] vertattrib;
	delete [] seg2end;
	freeSubQuadtree(root);
}


void PolyQuadtree::creV2seg(void){

	int v0,v1,i;
	for(i=0; i<numseg; i++){
		v0=seg2end[i][0];
		v1=seg2end[i][1];
		if(v0<0||v0>=numvert||v1<0||v1>=numvert)
			jf_error("crev2seg");
		v2seg[v0][1]=i;
		v2seg[v1][0]=i;
	}
}

void PolyQuadtree::insertVertInSubTree(int v, CellNode2D *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(!ifPointOverlapWithBox2D(vert[v],cnode->bound,epsoverlap ))
		return;
	if(!cnode->isLeaf()){
		for(int i=0; i<4; i++)
			insertVertInSubTree(v,cnode->child[i]);
		return;
	}
	if(cnode->vertincell==-1){
		cnode->vertincell=v;
	}else{
		splitNode(cnode);
		for(int i=0; i<4; i++)
			insertVertInSubTree(v,cnode->child[i]);
	}
}
bool ifPointOverlapWithBox2D(double p[2],double bd[4],double eps=0){

	double bound[4];

	double a=bd[2]-bd[0];
	double b=bd[3]-bd[1];
	bound[0]=bd[0]-eps*a;
	bound[1]=bd[1]-eps*b;
	bound[2]=bd[2]+eps*a;
	bound[3]=bd[3]+eps*b;
	if(p[0]>=bound[0]&&p[1]>=bound[1]&&p[0]<=bound[2]&&p[1]<=bound[3])
		return true;
	else
		return false;
}

void PolyQuadtree::splitNode(CellNode2D *cnode){

	for(int i=0; i<4; i++){
		cnode->child[i]=new CellNode2D(cnode->bound);
		cnode->child[i]->parent=cnode;
	}

	double x=(cnode->bound[0]+cnode->bound[2])/2.;
	double y=(cnode->bound[1]+cnode->bound[3])/2.;
	cnode->child[0]->bound[2]=x; cnode->child[0]->bound[1]=y;
	cnode->child[1]->bound[0]=x; cnode->child[1]->bound[1]=y;
	cnode->child[2]->bound[2]=x; cnode->child[2]->bound[3]=y;
	cnode->child[3]->bound[0]=x; cnode->child[3]->bound[3]=y;
	if(cnode->vertincell==-1)
		return;
	for(int i=0; i<4; i++)
		insertVertInSubTree(cnode->vertincell,cnode->child[i]);
	cnode->vertincell=-1;
}
void PolyQuadtree::insertSegInSubTree(int seg,CellNode2D *cnode){

	if(!cnode)
		jf_error("insertseginsubtree");
//	if(cnode->bound[0]<=15.219837109376&&cnode->bound[0]>=15.219837109374)
//		cout<<2;
	if(!ifSegOverlapBox2D(vert[seg2end[seg][0]],vert[seg2end[seg][1]],cnode->bound,epsoverlap))
		return;
	if(!cnode->isLeaf()){
		for(int i=0; i<4; i++)
			insertSegInSubTree(seg,cnode->child[i]);
		return;
	}
	if(cnode->vertincell==seg2end[seg][0]||cnode->vertincell==seg2end[seg][1])
		return;
	if(!cnode->psegar)
		cnode->psegar =new vector<int>;
	cnode->psegar->push_back(seg);
}

void PolyQuadtree::compVertattrib(void){


	for(int i=0 ;i<numvert; i++){
		int vf=seg2end[v2seg[i][0]][0];
		int vb=seg2end[v2seg[i][1]][1];
		if(i==11821) 
			i=i;
		vertattrib[i]=-convexityOf3Point(vert[vf],vert[i],vert[vb],epscoplanar);
	}
}

int convexityOf3Point(double p0[2],double p1[2],double p2[2],double eps=0){

	double d,p10[2],p12[2];

	vec_2p2(p1,p0,p10);
	vec_2p2(p1,p2,p12);
	d=vec_crop2(p10,p12);
	if(fabs(d)<=eps&&vec_dotp2(p10,p12)<0) return 0;
	if(d>0) 
		return -1;
	else  
		return 1;
}
bool ifSegOverlapBox2D(double ps[2],double pe[2],double bd[4],double eps){

	double bound[4];

	double a=bd[2]-bd[0];
	double b=bd[3]-bd[1];
	bound[0]=bd[0]-eps*a;
	bound[1]=bd[1]-eps*b;
	bound[2]=bd[2]+eps*a;
	bound[3]=bd[3]+eps*b;

	if(ps[0]<bound[0]&&pe[0]<bound[0]||ps[0]>bound[2]&&pe[0]>bound[2]||
	   ps[1]<bound[1]&&pe[1]<bound[1]||ps[1]>bound[3]&&pe[1]>bound[3]) 
	   return false;
	else if(ps[0]>=bound[0]&&ps[1]>=bound[1]&&ps[0]<=bound[2]&&ps[1]<=bound[3]
	 ||pe[0]>=bound[0]&&pe[1]>=bound[1]&&pe[0]<=bound[2]&&pe[1]<=bound[3])
		return true;
	else
		return if4CornerOfBoxAtDifferentSideOfSeg(ps,pe,bound);
}
bool if4CornerOfBoxAtDifferentSideOfSeg(double ps[2],double pe[2],double bound[4]){

 	double pse[2],psp[2],p[2],d,dc;

	vec_2p2(ps,pe,pse);
	p[0]=bound[0]; p[1]=bound[1];
	vec_2p2(ps,p,psp);
	d=vec_crop2(psp,pse);

	p[0]=bound[2]; p[1]=bound[1];
	vec_2p2(ps,p,psp);
	dc=vec_crop2(psp,pse);
	if(dc*d<0) return true;

	p[0]=bound[2]; p[1]=bound[3];
	vec_2p2(ps,p,psp);
	dc=vec_crop2(psp,pse);
	if(dc*d<0) return true;

	p[0]=bound[0]; p[1]=bound[3];
	vec_2p2(ps,p,psp);
	dc=vec_crop2(psp,pse);
	if(dc*d<0) return true;
	
	return false;
}
void boxOfPoints( double (*p)[2] , int num ,double box[4] ){

  int i , j ;
  double a ;

  if( num<=0 ) jf_error("boxofP" ) ;
  for(i=0 ; i<2 ; i++ )
    box[i]=box[i+2]=p[0][i] ;
  for(j=1 ; j<num ; j++ ){
    for( i=0 ; i<2 ; i++ ){
	 if( p[j][i]<box[i] ) box[i]=p[j][i] ;
	 if( p[j][i]>box[i+2] ) box[i+2]=p[j][i] ;
    }
  }
  a=max( box[2]-box[0] ,box[3]-box[1] ) ;
  for( i=0 ; i<2 ; i++ ){
    box[i] -= 0.01*a ;
    box[i+2] += 0.01*a ;
  }
}
void PolyQuadtree::getCellSeqWithUnknownAttribFromaCell(CellNode2D *cnode,vector<CellNode2D *>  * &pcellseq,
														CellNode2D * &pcellm,int &ia,double pm[2]){

	CellNode2D *pcellt,*pcell;

	if(cnode==0) return;
	pcellseq=new vector<CellNode2D *>;
	pcellseq->push_back(cnode);
	pcellt=cnode;
	for(;;){
		pm[0]=pcellt->bound[0]; pm[1]=pcellt->bound[1];
		pcell=getTheNeighbOfCellAtSpeciDirectWithRefPoint(pcellt,-1,0,pm);
		if(pcell==0){
			ia=-1; pcellm=0;
			return;
		}else if(pcell->inoutattrib!=-2){
			ia=pcell->inoutattrib; pcellm=pcell;
			return;
		}
		pcellseq->push_back(pcell);
		pcellt=pcell;
	}
}

void PolyQuadtree::freeSubQuadtree(CellNode2D *pcell){

	if(pcell ==0) return;
	for(int i=0; i<4; i++)
		freeSubQuadtree(pcell->child[i]);
	delete pcell;
}

CellNode2D ::CellNode2D(double bd[4]){

	psegar=0;
	vertincell=-1;
	inoutattrib=-2;
	for(int i=0; i<4; i++)
		bound[i]=bd[i];
	child[0]=child[1]=child[2]=child[3]=0;
	parent=0;
}

CellNode2D ::~CellNode2D (){
	delete psegar;
}

void PIP2D_jianfei_cpp(double *vertices, int *numV,
                       double *query,    int *numQ,
                       int    *result) {
//Rprintf("PIP_jianfei_cpp: entered C++ function.\n");
	// Transfer vertex data from flat double array to double[3] array.
	double (*vert)[2];
	vert = ( double (*)[2]) new double [2*(*numV)];
#if 0
	double minX = DLB_MAX, minY=DBL_MAX;
#else
	double minX = FLT_MAX, minY=FLT_MAX;
#endif
	int i;
    for ( i = 0; i < (*numV); i++ ) {
        vert[i][0] = vertices[i+0*(*numV)];
        vert[i][1] = vertices[i+1*(*numV)];
 
        if ( minX > vert[i][0] ) minX = vert[i][0];
        if ( minY > vert[i][1] ) minY = vert[i][1];
    }
 
    // Set minimum coordinates to 0.
    for ( i = 0; i < (*numV); i++ ) {
        vert[i][0] -= minX;
        vert[i][1] -= minY;
    }
 
    // Attempt to construct a object of PointInPolyhedron
    PolyQuadtree *ptpoly = 0;

    try {
 	    ptpoly = new PolyQuadtree(vert,(*numV));
    }
    catch ( int ptpolyError ) {
        // Fill result vector with the code "-2" to indicate
        // a failed initialization.
        for( i=0; i<(*numQ); i++) {
            result[i] = -ptpolyError;
        }

        // Revert XYZ coordinates back to original values.
        for ( i = 0; i < (*numV); i++ ) {
            vert[i][0] += minX;
            vert[i][1] += minY;
        }
        return;
    }

    // Loop over queries, feed them to the Jianfei method.
    // Don't forget about the minX, minY, and minZ shifts.
    double q[2]={0,0};
    for( i=0; i<(*numQ); i++) {
        q[0]      = query[i+0*(*numQ)] - minX;
        q[1]      = query[i+1*(*numQ)] - minY;
        result[i] = ptpoly->isPinpolygon(q);
    }

    // RELEASE MEMORY!!
    delete [] vert;
    delete ptpoly;
}

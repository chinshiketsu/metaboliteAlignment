// align.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
using namespace std;
//#define  INPUT_OVERVIEW	
#define  ALIGN_OUTPUT
struct Compound {
	int cID;
	string cName;
	double Mz, Rt, Area;
	Compound(double mz,double rt, double area):Mz(mz),Rt(rt),Area(area){}
	Compound(int id, double mz, double rt, double area) :cID(id), Mz(mz), Rt(rt), Area(area) {}
	Compound(){}
	string toString() {
		char buf[256] = {0};
		sprintf(buf,"ID=%d  m/z=%.10f  Rt=%.10f  Area=%.10f\n", cID,Mz, Rt, Area);
		return string(buf);
	}
};
struct AreaCmp {
	bool operator() (const Compound& c1, const Compound& c2) const{
		return c1.Area > c2.Area;
	}
};
struct MzCmp {
	bool operator() (const Compound& c1, const Compound& c2) const{
		return c1.Mz < c2.Mz;
	}
};

bool CmpByIncreasingRt(Compound cpd1, Compound cpd2) {
	return cpd1.Rt < cpd2.Rt;
}

typedef multiset<Compound, AreaCmp>  AreaOrderedSet;
typedef multiset<Compound, MzCmp>    MzOrderedSet;
typedef multiset<Compound>::iterator CompoundItr;
//typedef set<Compound, AreaCmp>::iterator AreaOrderedItr;
//typedef set<Compound, MzCmp>::iterator   MzOrderedItr;
typedef tuple<int, vector<Compound> > CpdCluster;

class UnalignedList {
public:
	AreaOrderedSet areaRank;
	MzOrderedSet   mzRank;
public:
	size_t size();
	bool empty();
	void insert(Compound x);
	void remove(const Compound& x);
	void mzSearch(double lb, double ub, CompoundItr& beginItr, CompoundItr& endItr);
	CompoundItr areaSearch(Compound x);
	Compound getTopArea();
};


void readDataset(UnalignedList& UnALst, const char* inputfile) {
	ifstream fin;
	//char inputfile[] = "dataset.txt";
	fin.open(inputfile, ifstream::in);
	double mz, rt, area;
	int id = 0;
	while (fin>>mz>>rt>>area) {
		Compound cpd = Compound(id++, mz, rt, area);
		UnALst.insert(cpd);
	}
	fin.close();
}

void generateTxtReport(char* filename) {
/*
	ofstream fout;
	fout.open(filename, ofstream::out);

	fout << alignList.size() << " aligned compounds in total." << endl;
	for (auto tuple : alignList) {
		fout << "AlignID=" << get<0>(tuple) << endl;
		auto vec = get<1>(tuple);
		for (auto cpd : vec)
			fout << cpd.toString();
		fout << endl;
	}


	fout.close();*/
}

//================== main function ==================
int main()
{	
	double AlignWindowPhase1  = 20;
	double AlignWindowPhase2  = 30;
	double MassTol = 1;

	ofstream fout;
	char outputfile[] = "report.txt";
	fout.open(outputfile, ofstream::out);

	UnalignedList    unalignList;
	vector<CpdCluster> alignList;
	readDataset(unalignList, "three_ds.txt");
#ifdef INPUT_OVERVIEW	
	for (auto cpd : unalignList.areaRank)
		fout << cpd.toString();
	fout << unalignList.size();
#endif	
	int AppearTimes = 2;
	int AlignID = 0;
	while (!unalignList.empty()) {
		Compound cpd = unalignList.getTopArea();
		double AlignRt       = cpd.Rt;
		double AlignMoleMass = cpd.Mz;
		
		vector<Compound> alignGroup;
		CompoundItr itrBegin, itrEnd;
		unalignList.mzSearch(AlignMoleMass - MassTol, AlignMoleMass + MassTol, itrBegin, itrEnd);//find Mass-Valid compounds 

		for (auto itr = itrBegin; itr != itrEnd; ++itr) {//find Rt-Valid compounds within Mass-Valid compounds
			auto cpd = *itr;
			if (AlignRt - AlignWindowPhase1*0.5 < cpd.Rt && cpd.Rt < AlignRt + AlignWindowPhase1*0.5) 
				alignGroup.push_back(cpd);
		}
		for (auto cpd : alignGroup)
			unalignList.remove(cpd);
		//cout << AlignID << endl;
		//for (auto item : alignGroup)
		//	cout << item.toString();
		//cout << "__________________________" << endl;

		int alignSize = alignGroup.size();
		if ((int)alignSize >= AppearTimes) {//only align the compounds when it is contained in all groups 
			sort(alignGroup.begin(), alignGroup.end(), CmpByIncreasingRt);
			AlignRt = alignGroup[alignSize / 2].Rt; //find the compound with median Rt
			AlignMoleMass = alignGroup[alignSize / 2].Mz;

			//enlargement of compounds after the AlignRt and the AlignMoleMass has been revised
			unalignList.mzSearch(AlignMoleMass - MassTol, AlignMoleMass + MassTol, itrBegin, itrEnd);
			for (auto itr = itrBegin; itr != itrEnd; ++itr) {
				auto cpd = *itr;
				if (AlignRt - AlignWindowPhase2*0.5 < cpd.Rt && cpd.Rt < AlignRt + AlignWindowPhase2*0.5)//The Rt bound changes
					alignGroup.push_back(cpd);
			}
		}

		//for (auto item : alignGroup)
		//	cout << item.toString();
		//cout << "_________________________((" << endl;

		//remove selected compounds
		for (auto cpd : alignGroup)
			unalignList.remove(cpd);

		//for (auto cpd : alignGroup)
		//	cout << cpd.toString()<<endl;
		if ((int)alignGroup.size() >= AppearTimes) {
			CpdCluster alignTuple = make_tuple(AlignID++, alignGroup);
			alignList.push_back(alignTuple);
		}
		alignGroup.clear();
		//if (AlignID==18)
		//	break;
	}

#ifdef ALIGN_OUTPUT
	fout << alignList.size() << " aligned compounds in total." << endl;
	for (auto tuple : alignList) {
		fout << "AlignID=" << get<0>(tuple) << endl;
		auto vec = get<1>(tuple);
		for (auto cpd : vec)
			fout << cpd.toString();
		fout << endl;
	}
#endif
	//system("pause");


	fout.close();
	
    return 0;
}

//==================== Unit Test ====================


//================ member funciton impl ==============
size_t UnalignedList::size() {
	return areaRank.size();
}

bool UnalignedList::empty() {
	return areaRank.empty();
}

void UnalignedList::insert(Compound x){
	areaRank.insert(x);
	mzRank.insert(x);
}

void UnalignedList::remove(const Compound & x){
	auto areaItr=areaRank.find(x);
	auto mzItr = mzRank.find(x);
	if (areaItr != areaRank.end())
		areaRank.erase(areaItr);
	if (mzItr != mzRank.end())
		mzRank.erase(mzItr);
}

void UnalignedList::mzSearch(double lb, double ub, CompoundItr & beginItr, CompoundItr & endItr){
	beginItr = mzRank.lower_bound(Compound(lb, 0, 0));
	endItr   = mzRank.upper_bound(Compound(ub, 0, 0));
}

CompoundItr UnalignedList::areaSearch(Compound x){
	return areaRank.find(x);
}

Compound UnalignedList::getTopArea(){
	return *(areaRank.begin());
}



#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include "Base.h"
#include "QCMOEA_PM.h"
#include "rand.h"
using namespace std;
using namespace Base;

QCMOEA_PM::QCMOEA_PM()
{
}

QCMOEA_PM::~QCMOEA_PM()
{
}

void QCMOEA_PM::SetParameters(int AN, int NumberOfExtractedFams, long TimeLimit, int AllJobTotalPTime)
{

	this->m_FamD = NumberOfExtractedFams;
	this->m_TimeLimit = TimeLimit;
	this->m_Temperature = 0.6 * AllJobTotalPTime / (this->m_Factories * this->m_Jobs * this->m_Machines);

	this->m_Popsize = 1000;
	this->m_RefSize = AN;
	this->m_PS1 = AN;
	this->m_PS2 = AN;
	this->thr_zeta = 1.0;

	m_CMOEAPopulation.clear();
	m_CMOEAPopulation.resize(m_RefSize);

	m_nadirpointMS = -1;
	m_nadirpointTEC = -1;

	m_idealpointMS = -1;
	m_idealpointTEC = -1;

	//档案集 Reference set
	m_RefSpanArray.clear();
	m_RefSpanArray.resize(m_RefSize);

	m_RefTECArray.clear();
	m_RefTECArray.resize(m_RefSize);

	m_RefSpeedVector.clear();
	m_RefSpeedVector.resize(m_RefSize);

	m_nRefCriFacArray.clear();
	m_nRefCriFacArray.resize(m_RefSize, 0);

	m_RefFacFamSeqArray.clear();
	m_RefFacFamSeqArray.resize(m_RefSize);

	m_RefJobSeqInFamArray.clear();
	m_RefJobSeqInFamArray.resize(m_RefSize);

	m_RefFacSpanArray.clear();
	m_RefFacSpanArray.resize(m_RefSize);

	m_RefFacECArray.clear();
	m_RefFacECArray.resize(m_RefSize);


	//组种群
	m_SpanArray1.clear();
	m_SpanArray1.resize(m_PS1);

	m_TECArray1.clear();
	m_TECArray1.resize(m_PS1);

	m_Map1.clear();
	m_Map1.resize(m_PS1);

	m_FacFamSeqArray.clear();
	m_FacFamSeqArray.resize(m_PS1);

	//工件种群
	m_SpanArray2.clear();
	m_SpanArray2.resize(m_PS2);

	m_TECArray2.clear();
	m_TECArray2.resize(m_PS2);

	m_Map2.clear();
	m_Map2.resize(m_PS2);

	m_JobSeqInFamArray.clear();
	m_JobSeqInFamArray.resize(m_PS2);
}

double QCMOEA_PM::GetTemperature()
{
	int Sum = 0; 
	for (int j = 0; j < this->m_Jobs; j++)
		for (int m = 0; m < this->m_Machines; m++)
			Sum += this->m_JobOperPTime[j][m];
	double t = 0.6 * Sum / (this->m_Factories * this->m_Jobs * this->m_Machines);
	return t;
}

void QCMOEA_PM::BasedindRandFamInFacTobestPos(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC,
	int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacTEC, int& ObjectMS, float& ObjectTEC)
{
	int CriFac = 0;
	int max = INT_MIN;
	for (int i = 0; i < FacFamSeq.size(); i++)
	{
		if (max < FacSpan[i])
		{
			max = FacSpan[i];
			CriFac = i;
		}
	}
	int BestFac, BestPos;
	float Ind = 0.0;
	vector<vector<int>> TempFacFamSeq = FacFamSeq;

	for (int fam = 0; fam < TempFacFamSeq[CriFac].size(); fam++)
	{
		if (FacFamSeq[CriFac].size() < 2)
			continue;

		int CurFam = TempFacFamSeq[CriFac][fam];

		vector<int>::iterator it = find(FacFamSeq[CriFac].begin(), FacFamSeq[CriFac].end(), CurFam);
		FacFamSeq[CriFac].erase(it);

		BestFac = -1;
		BestPos = -1;
		Ind = this->FindBestPosToInsertFamForAllFac_Ind(FacFamSeq, JobSeqInFam, SpeedMatrix, CurFam, BestFac, BestPos,
			nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacTEC, ObjectMS, ObjectTEC, CMOEAPopulation);
		FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);

	}
}

void QCMOEA_PM::BasedindSwapFam(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int OrgMS = ObjectMS;
	float OrgTEC = ObjectTEC;
	if (OrgMS > nadirpointMS)
		nadirpointMS = OrgMS;
	if (OrgMS < idealpointMS)
		idealpointMS = OrgMS;
	if (OrgTEC > nadirpointTEC)
		nadirpointTEC = OrgTEC;
	if (OrgTEC < idealpointTEC)
		idealpointTEC = OrgTEC;
	float OrgnormalMS = -1;
	float OrgnormalTEC = -1;

	OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	float Orgconvergence_ind = 0;
	int prank = Individual_PM::GetPRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC, OrgMS, OrgTEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	Orgconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);


	for (int i = 0; i < m_Families / m_Factories; i++)
	{
		int crifac = -1;
		int max = -1;
		int optfac = -1;
		int min = INT_MAX;
		for (int i = 0; i < FacFamSeq.size(); i++)
		{
			if (max < FacSpan[i])
			{
				max = FacSpan[i];
				crifac = i;
			}
			if (min > FacSpan[i])
			{
				min = FacSpan[i];
				optfac = i;
			}
		}
		if (crifac == optfac)
		{
			break;
		}

		if (FacFamSeq[crifac].empty())
		{
			break;
		}
		int Pos1 = wyt_rand(FacFamSeq[crifac].size());
		int Fam1 = FacFamSeq[crifac][Pos1];
		FacFamSeq[crifac].erase(FacFamSeq[crifac].begin() + Pos1);

		if (FacFamSeq[optfac].empty())
		{
			break;
		}
		int Pos2 = rand() % FacFamSeq[optfac].size();
		int Fam2 = FacFamSeq[optfac][Pos2];
		FacFamSeq[optfac].erase(FacFamSeq[optfac].begin() + Pos2);

		FacFamSeq[crifac].insert(FacFamSeq[crifac].begin() + Pos1, Fam2);
		FacSpan[crifac] = GetCTimeForPerFac(FacFamSeq[crifac], JobSeqInFam, SpeedMatrix);
		FacEC[crifac] = GetTECForPerFac(FacFamSeq[crifac], JobSeqInFam, SpeedMatrix);

		FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, Fam1);
		FacSpan[optfac] = GetCTimeForPerFac(FacFamSeq[optfac], JobSeqInFam, SpeedMatrix);
		FacEC[optfac] = GetTECForPerFac(FacFamSeq[optfac], JobSeqInFam, SpeedMatrix);

		int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
		float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0);

		if (AfterMakespan > nadirpointMS)
			nadirpointMS = AfterMakespan;
		if (AfterMakespan < idealpointMS)
			idealpointMS = AfterMakespan;

		if (AfterTEC > nadirpointTEC)
			nadirpointTEC = AfterTEC;
		if (AfterTEC < idealpointTEC)
			idealpointTEC = AfterTEC;

		bool flag = true;
		if ((AfterMakespan < OrgMS && AfterTEC <= OrgTEC) || (AfterMakespan <= OrgMS && AfterTEC < OrgTEC))
		{
			for (int i = 0; i < CMOEAPopulation.size(); i++)
			{
				if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				Individual_PM tempIndi;
				tempIndi.m_FacFamSeqArray = FacFamSeq;
				tempIndi.m_JobSeqInFamArray = JobSeqInFam;
				tempIndi.MS = AfterMakespan;
				tempIndi.TEC = AfterTEC;
				tempIndi.m_SpeedVector = SpeedMatrix; 
				CMOEAPopulation.push_back(tempIndi);
			}
		}
	
		float AfternormalMS = -1;
		float AfternormalTEC = -1;

		AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
		AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

		float Afterconvergence_ind = 0.0;
		int prank = Individual_PM::GetPRank(CMOEAPopulation, AfternormalMS, AfternormalTEC, AfterMakespan, AfterTEC);
		pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, AfternormalMS, AfternormalTEC);
		int lrank = re.first;
		int nbrmax = re.second;
		Afterconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

		if (Afterconvergence_ind < Orgconvergence_ind)
		{
			OrgMS = AfterMakespan;
			OrgTEC = AfterTEC;
			Orgconvergence_ind = Afterconvergence_ind;
			ObjectMS = AfterMakespan;
			ObjectTEC = AfterTEC;
		}
		else
		{
			FacFamSeq[optfac].erase(FacFamSeq[optfac].begin() + Pos2);
			FacFamSeq[crifac].erase(FacFamSeq[crifac].begin() + Pos1);
			FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, Fam2);
			FacFamSeq[crifac].insert(FacFamSeq[crifac].begin() + Pos1, Fam1);
		}
	}

}

void QCMOEA_PM::BasedindLS_SetupInsert(vector<vector <int>>& FacFamSeq, vector <vector <int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int fac = 0;
	for (int i = 0; i < FacFamSeq.size(); i++)
	{
		if (ObjectMS == FacSpan[i])
		{
			fac = i;
		}
	}
	int CurFam = 0;
	int NowFam = 0;
	int Max_Setup = 0;
	int Pos = 0;
	vector<vector<int>> TempFacFamSeq = FacFamSeq;

	for (int fam = 0; fam < TempFacFamSeq[fac].size(); fam++)
	{
		int PreFam = FacFamSeq[fac][0];
		for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++)
		{

			CurFam = FacFamSeq[fac][Fam];
			if (this->m_SetupTime[0][PreFam][CurFam] > Max_Setup)
			{
				Max_Setup = this->m_SetupTime[0][PreFam][CurFam];
				Pos = Fam;
				NowFam = CurFam;
			}
			PreFam = CurFam;
		}

	}
	FacFamSeq[fac].erase(FacFamSeq[fac].begin() + Pos);
	int bestFac = -1, bestPos = -1;
	float ind = FindBestPosToInsertFamForAllFac_Ind(FacFamSeq, JobSeqInFam, SpeedMatrix, NowFam, bestFac, bestPos,
		nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacEC, ObjectMS, ObjectTEC, CMOEAPopulation);
	FacFamSeq[bestFac].insert(FacFamSeq[bestFac].begin() + bestPos, NowFam);

}

void QCMOEA_PM::BasedindLS_SetupSwap(vector<vector <int>>& FacFamSeq, vector <vector <int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int OrgMS = ObjectMS;
	float OrgTEC = ObjectTEC;

	if (OrgMS > nadirpointMS)
		nadirpointMS = OrgMS;
	if (OrgMS < idealpointMS)
		idealpointMS = OrgMS;

	if (OrgTEC > nadirpointTEC)
		nadirpointTEC = OrgTEC;
	if (OrgTEC < idealpointTEC)
		idealpointTEC = OrgTEC;

	float OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	float OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	float Orgconvergence_ind = 0;
	int prank = Individual_PM::GetPRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC, OrgMS, OrgTEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	Orgconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

	for (int i = 0; i < m_Families / m_Factories; i++)
	{
		int CriMS = 0;
		int OptMS = INT_MAX;
		int crifac = 0, optfac = 0;
		for (int i = 0; i < FacFamSeq.size(); i++)
		{
			if (CriMS < FacSpan[i])
			{
				crifac = i; 
				CriMS = FacSpan[i]; 
			}
			if (OptMS >= FacSpan[i]) 
			{
				optfac = i; 
				OptMS = FacSpan[i]; 
			}
		}

		int CurFam;
		int PreFam;
		int NowFam1;
		int NowFam2;
		int Max_Setup = 0;
		int Pos1 = 0;
		int Pos2 = 0;
		int Min_Setup = INT_MAX;


		PreFam = FacFamSeq[crifac][0];
		for (int Fam = 0; Fam < FacFamSeq[crifac].size(); Fam++)
		{
			CurFam = FacFamSeq[crifac][Fam];
			if (this->m_SetupTime[0][PreFam][CurFam] > Max_Setup)
			{
				Max_Setup = this->m_SetupTime[0][PreFam][CurFam];
				Pos1 = Fam;
				NowFam1 = CurFam;
			}
			PreFam = CurFam;
		}

		PreFam = FacFamSeq[optfac][0];
		for (int Fam = 0; Fam < FacFamSeq[optfac].size(); Fam++)
		{
			CurFam = FacFamSeq[optfac][Fam];
			if (this->m_SetupTime[0][PreFam][CurFam] < Min_Setup)
			{
				Min_Setup = this->m_SetupTime[0][PreFam][CurFam];
				Pos2 = Fam;
				NowFam2 = CurFam;
			}
			PreFam = CurFam;
		}

		FacFamSeq[crifac].erase(FacFamSeq[crifac].begin() + Pos1);
		FacFamSeq[optfac].erase(FacFamSeq[optfac].begin() + Pos2);

		FacFamSeq[crifac].insert(FacFamSeq[crifac].begin() + Pos1, NowFam2);
		FacSpan[crifac] = GetCTimeForPerFac(FacFamSeq[crifac], JobSeqInFam, SpeedMatrix);
		FacEC[crifac] = GetTECForPerFac(FacFamSeq[crifac], JobSeqInFam, SpeedMatrix);

		FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, NowFam1);
		FacSpan[optfac] = GetCTimeForPerFac(FacFamSeq[optfac], JobSeqInFam, SpeedMatrix);
		FacEC[optfac] = GetTECForPerFac(FacFamSeq[optfac], JobSeqInFam, SpeedMatrix);

		int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
		float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0);

		if (AfterMakespan > nadirpointMS)
			nadirpointMS = AfterMakespan;
		if (AfterMakespan < idealpointMS)
			idealpointMS = AfterMakespan;
		if (AfterTEC > nadirpointTEC)
			nadirpointTEC = AfterTEC;
		if (AfterTEC < idealpointTEC)
			idealpointTEC = AfterTEC;

		bool flag = true;

		if ((AfterMakespan < OrgMS && AfterTEC <= OrgTEC) || (AfterMakespan <= OrgMS && AfterTEC < OrgTEC))
		{
			for (int i = 0; i < CMOEAPopulation.size(); i++)
			{
				if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				Individual_PM tempIndi;
				tempIndi.m_FacFamSeqArray = FacFamSeq;
				tempIndi.m_JobSeqInFamArray = JobSeqInFam;
				tempIndi.MS = AfterMakespan;
				tempIndi.TEC = AfterTEC;
				tempIndi.m_SpeedVector = SpeedMatrix; 
				CMOEAPopulation.push_back(tempIndi);
			}
		}

		float AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
		float AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
		Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
		int prank = Individual_PM::GetPRank(CMOEAPopulation, AfternormalMS, AfternormalTEC, AfterMakespan, AfterTEC);
		pair<int, int> re1 = Individual_PM::GetLRank(CMOEAPopulation, AfternormalMS, AfternormalTEC);
		lrank = re1.first;
		nbrmax = re1.second;
		float Afterconvergence_ind = Afterconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

		if (Afterconvergence_ind < Orgconvergence_ind)
		{
			OrgMS = AfterMakespan;
			OrgTEC = AfterTEC;
			Orgconvergence_ind = Afterconvergence_ind;
			ObjectMS = AfterMakespan;
			ObjectTEC = AfterTEC;
		}
		else
		{
			FacFamSeq[crifac].erase(FacFamSeq[crifac].begin() + Pos1);
			FacFamSeq[optfac].erase(FacFamSeq[optfac].begin() + Pos2);
			FacFamSeq[crifac].insert(FacFamSeq[crifac].begin() + Pos1, NowFam1);
			FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, NowFam2);
		}

	}
}

void QCMOEA_PM::Basedind_DCFams(vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int Len = 2;
	int CriFac;
	float max = FLT_MIN;

	for (int i = 0; i < FacFamSeq.size(); i++)
	{
		if (max < FacEC[i])
		{
			max = FacEC[i];
			CriFac = i;
		}
	}
	vector<int> ExtractFamSeq;
	bool groupSelected = false;
	do {
		int Fac;
		if (ExtractFamSeq.size() < Len / 2) 
			Fac = CriFac;
		else  
			Fac = rand() % m_Factories;
		if (FacFamSeq[Fac].size() > 1)
		{
			int Pos = rand() % FacFamSeq[Fac].size();
			ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]); 
			FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);
			groupSelected = true;
		}
	} while (ExtractFamSeq.size() < Len && groupSelected);
	for (int Fam = 0; Fam < ExtractFamSeq.size(); Fam++)
	{

		int CurFam = ExtractFamSeq[Fam]; 
		int bestFac = -1, bestPos = -1; 

		if (Fam == ExtractFamSeq.size() - 1)
		{
			float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR(CMOEAPopulation, FacFamSeq, JobSeqInFam, SpeedMatrix, CurFam, bestFac, bestPos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
		}
		else
		{
			float randValue = wyt_rand(0, 1);

			if (randValue < 0.5)
			{
				FindBestPosToInsertFamForAllFac_Makespan(FacFamSeq, JobSeqInFam, SpeedMatrix, CurFam, bestFac, bestPos);
			}
			else
			{
				FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, SpeedMatrix, CurFam, bestFac, bestPos);
			}
		}
		FacFamSeq[bestFac].insert(FacFamSeq[bestFac].begin() + bestPos, CurFam);
	}

}


/*************************************工件进化*********************************************/

void QCMOEA_PM::BasedindCirJobInFamTobestPos(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC,
	int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int fac = -1;
	for (int i = 0; i < FacFamSeq.size(); i++)
	{
		if (ObjectMS == FacSpan[i])
		{
			fac = i;
		}
	}

	for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++)
	{
		int CurFam = FacFamSeq[fac][Fam];
		if (JobSeqInFam[CurFam].size() > 1)
		{
			Basedind_JobsInFam(fac, CurFam, FacFamSeq, JobSeqInFam, SpeedMatrix, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		}
	}
}

void QCMOEA_PM::BasedindSwapJob(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC,
	int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	int OrgMS = ObjectMS;
	float OrgTEC = ObjectTEC;

	if (OrgMS > nadirpointMS)
		nadirpointMS = OrgMS;
	if (OrgMS < idealpointMS)
		idealpointMS = OrgMS;

	if (OrgTEC > nadirpointTEC)
		nadirpointTEC = OrgTEC;
	if (OrgTEC < idealpointTEC)
		idealpointTEC = OrgTEC;

	float OrgnormalMS = -1;
	float OrgnormalTEC = -1;

	OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	float Orgconvergence_ind = 0.0;
	int prank = Individual_PM::GetPRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC, OrgMS, OrgTEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	Orgconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

	int fac = -1;
	int min = -1;

	for (int i = 0; i < FacFamSeq.size(); i++)
	{
		if (min < FacSpan[i])
		{
			min = FacSpan[i];
			fac = i;
		}
	}
	for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++)
	{
		int CurFam = FacFamSeq[fac][Fam]; 
		if (JobSeqInFam[CurFam].size() > 1)
		{
			int pt1 = 0;
			int pt2 = 0;

			if (JobSeqInFam[CurFam].size() > 2)
			{
				pt1 = rand() % JobSeqInFam[CurFam].size();
				do
				{
					pt2 = rand() % JobSeqInFam[CurFam].size();
				} while (pt1 == pt2);

				if (pt1 > pt2)
				{
					
					int t = pt1;
					pt1 = pt2;
					pt2 = t;
				}
			}
			else if (JobSeqInFam[CurFam].size() == 2)
			{
				pt1 = 0;
				pt2 = 1;
			}
			vector<vector<int>> TempJobSeqInFam = JobSeqInFam;

			int Job = JobSeqInFam[CurFam][pt1];
			JobSeqInFam[CurFam][pt1] = JobSeqInFam[CurFam][pt2];
			JobSeqInFam[CurFam][pt2] = Job;

			FacSpan[fac] = GetCTimeForPerFac(FacFamSeq[fac], JobSeqInFam, SpeedMatrix);

			int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
			if (AfterMakespan > nadirpointMS)
				nadirpointMS = AfterMakespan;
			if (AfterMakespan < idealpointMS)
				idealpointMS = AfterMakespan;


			float AfterTEC = 0;
			FacEC[fac] = GetTECForPerFac(FacFamSeq[fac], JobSeqInFam, SpeedMatrix);
			AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0);

			if (AfterTEC > nadirpointTEC)
				nadirpointTEC = AfterTEC;
			if (AfterTEC < idealpointTEC)
				idealpointTEC = AfterTEC;
			bool flag = true;
			if ((AfterMakespan <= OrgMS && AfterTEC < OrgTEC) || (AfterMakespan <= OrgMS && AfterTEC < OrgTEC))
			{
				for (int i = 0; i < CMOEAPopulation.size(); i++)
				{
					if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					Individual_PM tempIndi;
					tempIndi.m_FacFamSeqArray = FacFamSeq;
					tempIndi.m_JobSeqInFamArray = JobSeqInFam;
					tempIndi.MS = AfterMakespan;
					tempIndi.TEC = AfterTEC;
					tempIndi.m_SpeedVector = SpeedMatrix; 
					CMOEAPopulation.push_back(tempIndi);
				}
			}
		
			float AfternormalMS = -1;
			float AfternormalTEC = -1;
			AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
			AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

			Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
			float Afterconvergence_ind = 0;
			int prank = Individual_PM::GetPRank(CMOEAPopulation, AfternormalMS, AfternormalTEC, AfterMakespan, AfterTEC);
			pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, AfternormalMS, AfternormalTEC);
			int lrank = re.first;
			int nbrmax = re.second;
			Afterconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

			if (Afterconvergence_ind < Orgconvergence_ind)
			{
				OrgMS = AfterMakespan;
				OrgTEC = AfterTEC;
				Orgconvergence_ind = Afterconvergence_ind;
				ObjectMS = AfterMakespan;
				ObjectTEC = AfterTEC;
			}
			else
			{
				JobSeqInFam = TempJobSeqInFam;
			}
		}

	}
}


/*************************************档案集进化*************************************************************/
void QCMOEA_PM::BasedindRefInsert(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC,
	int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	calc_distribution_ind(CMOEAPopulation);
	int OrgMS = ObjectMS;
	float OrgTEC = ObjectTEC;
	if (OrgMS > nadirpointMS)
		nadirpointMS = OrgMS;
	if (OrgMS < idealpointMS)
		idealpointMS = OrgMS;

	if (OrgTEC > nadirpointTEC)
		nadirpointTEC = OrgTEC;
	if (OrgTEC < idealpointTEC)
		idealpointTEC = OrgTEC;

	float OrgnormalMS = -1;
	float OrgnormalTEC = -1;

	OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	int prank = Individual_PM::GetPRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC, OrgMS, OrgTEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, OrgnormalMS, OrgnormalTEC);;
	int lrank, nbrmax;
	lrank = re.first;
	nbrmax = re.second;
	float Orgconvergence_ind = 0;
	Orgconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);

	vector<int> FamsExtracted;
	unordered_map<int, vector<int>> JobsExtracted;
	vector<vector<int>> TempFacFamSeq = FacFamSeq;
	vector<vector<int>> TempJobSeqInFam = JobSeqInFam;
	vector<int> m_FamDArray = { 2, 3, 4, 5, 6 };
	this->m_FamD = m_FamDArray[rand() % m_FamDArray.size()];  // 随机选择
	BasedindDR_FamsAndJobs(CMOEAPopulation, FacFamSeq, JobSeqInFam, SpeedMatrix, FamsExtracted, JobsExtracted, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, m_FamD);
	vector<int> FacSpan(FacFamSeq.size(), 0);
	vector<float> FacEC(FacFamSeq.size(), 0);
	float AfterTEC = 0.0;
	int AfterMakespan = 0;

	GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, SpeedMatrix, FacSpan, FacEC, AfterMakespan, AfterTEC);


	if (AfterMakespan > nadirpointMS)
		nadirpointMS = AfterMakespan;
	if (AfterMakespan < idealpointMS)
		idealpointMS = AfterMakespan;

	if (AfterTEC > nadirpointTEC)
		nadirpointTEC = AfterTEC;
	if (AfterTEC < idealpointTEC)
		idealpointTEC = AfterTEC;

	bool flag = true;
	if ((AfterMakespan <= OrgMS && AfterTEC < OrgTEC) || (AfterMakespan < OrgMS && AfterTEC <= OrgTEC))
	{
		for (int i = 0; i < CMOEAPopulation.size(); i++)
		{
			if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			Individual_PM tempIndi;
			tempIndi.m_FacFamSeqArray = FacFamSeq;
			tempIndi.m_JobSeqInFamArray = JobSeqInFam;
			tempIndi.MS = AfterMakespan;
			tempIndi.TEC = AfterTEC;
			tempIndi.m_SpeedVector = SpeedMatrix; 
			CMOEAPopulation.push_back(tempIndi);
			m_FamDArray.push_back(m_FamD);
		}
	}

	float AfternormalMS = -1;
	float AfternormalTEC = -1;

	AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	float Afterconvergence_ind = 0;
	int Afterlrank, Afternbrmax;
	int Afterprank = Individual_PM::GetPRank(CMOEAPopulation, AfternormalMS, AfternormalTEC, AfterMakespan, AfterTEC);
	pair<int, int> re1 = Individual_PM::GetLRank(CMOEAPopulation, AfternormalMS, AfternormalTEC);
	Afterlrank = re1.first;
	Afternbrmax = re1.second;
	Afterconvergence_ind = 2 * (Afterprank - 1) + (1.0 / (Afternbrmax + 1)) * (Afterlrank - 1);

	if (Afterconvergence_ind < Orgconvergence_ind)
	{
		OrgMS = AfterMakespan;
		OrgTEC = AfterTEC;
		Orgconvergence_ind = Afterconvergence_ind;
		ObjectMS = AfterMakespan;
		ObjectTEC = AfterTEC;
	}
	else
	{
		FacFamSeq = TempFacFamSeq;
		JobSeqInFam = TempJobSeqInFam;
	}

}

void QCMOEA_PM::BasedindDR_FamsAndJobs(vector<Individual_PM>& CMOEAPopulation, vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FamsExtracted,
	unordered_map<int, vector<int>>& JobsExtracted, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int Fam_D)
{
	unordered_map<int, pair<int, int>> FamPosErasedFromFac;
	auto UpdateFamErasedFromFac = [&FamPosErasedFromFac](int Pos, int fac)
	{
		if (FamPosErasedFromFac.find(fac) == end(FamPosErasedFromFac))
		{
			FamPosErasedFromFac[fac] = { Pos, Pos - 1 };
		}
		else
		{
			if (Pos < FamPosErasedFromFac[fac].first)
			{
				FamPosErasedFromFac[fac].first = Pos;
			}
			if (Pos >= FamPosErasedFromFac[fac].second)
			{
				FamPosErasedFromFac[fac].second = Pos - 1;
			}
			else
			{
				FamPosErasedFromFac[fac].second = FamPosErasedFromFac[fac].second - 1;
			}
		}
	};

	FamsExtracted.clear();
	JobsExtracted.clear();

	while (FamsExtracted.size() < this->m_FamD)
	{
		int fac;
		do
		{
			fac = rand() % this->m_Factories;
		} while (Sol[fac].size() <= 1);

		int Pos = rand() % Sol[fac].size();
		int FamExt = Sol[fac][Pos];
		Sol[fac].erase(Sol[fac].begin() + Pos);
		FamsExtracted.push_back(FamExt);
		UpdateFamErasedFromFac(Pos, fac);
	}


	for (int i = 0; i < FamsExtracted.size(); i++)
	{
		int Fam = FamsExtracted[i];
		if (JobSeqInFam[Fam].size() >= 3)
		{
			for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++)
			{
				int JobPos = rand() % JobSeqInFam[Fam].size();
				int Job = JobSeqInFam[Fam][JobPos];
				JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
				JobsExtracted[Fam].push_back(Job);
			}
		}
	}


	while (FamsExtracted.size() > 0)
	{
		int BestFac = -1;
		int BestPos = -1;
		int Pos = rand() % FamsExtracted.size();
		int CurFam = FamsExtracted[Pos];
		FamsExtracted.erase(FamsExtracted.begin() + Pos);
		float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR(CMOEAPopulation, Sol, JobSeqInFam, SpeedMatrix, CurFam, BestFac, BestPos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
		Sol[BestFac].insert(Sol[BestFac].begin() + BestPos, CurFam);
		int FamPos = find(begin(Sol[BestFac]), end(Sol[BestFac]), CurFam) - begin(Sol[BestFac]);

		for (int i = 0; i < JobsExtracted[CurFam].size(); i++)
		{
			int BestJobPos = -1;
			float minFacInd = INT_MAX;
			int CurJob = JobsExtracted[CurFam][i];
			vector<vector<int>> TempJobSeqInFam = JobSeqInFam;
			for (int Pos = 0; Pos <= JobSeqInFam[CurFam].size(); Pos++)
			{
				TempJobSeqInFam[CurFam].insert(TempJobSeqInFam[CurFam].begin() + Pos, CurJob);
				int SpanAfterInsert = this->GetCTimeForPerFac(Sol[BestFac], TempJobSeqInFam, SpeedMatrix);
				float FacInd = GetindForPerFacAfterInsertJob_DR(CMOEAPopulation, BestFac, Sol, Sol[BestFac], JobSeqInFam, SpeedMatrix, SpanAfterInsert, CurFam, FamPos, JobsExtracted[CurFam][i], Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
				TempJobSeqInFam[CurFam].erase(TempJobSeqInFam[CurFam].begin() + Pos);
				if (FacInd < minFacInd)
				{
					minFacInd = FacInd;
					BestPos = Pos;
				}
			}
			JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);
		}
	}
}

/*************************************初始化*************************************************************/
void QCMOEA_PM::InitialPop()
{
	vector<vector<int>> speed_matrix_0 = vector<vector<int>>(m_Jobs, vector<int>(m_Machines, 0));

	//初始化档案集 Initialize Reference set and Set flags
	for (int PS = 0; PS < m_RefSize; PS++)//
	{
		vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//组，工厂完工时间
		vector<float> FacEC(m_Factories);
		vector<vector<int>> JobSeqInFam, FacFamSeq; //工件在组的序列，组在工厂的序列，新的组在工厂的序列

		//第一个参考解用LPT，其它随机生成 Generate job sequence in each family and family sequence using LPT
		if (PS == 0 || PS == 1)
		{
			this->SortJobsInFam(0, JobSeqInFam); //LPT
			this->SortFam(2, FamPrmu); //组平均准备时间

		}
		else if (PS == 2 || PS == 3)
		{
			this->SortJobsInFam(1, JobSeqInFam); //SPT
			this->SortFam(2, FamPrmu); //组平均准备时间
		}
		else // 随机生成工件序列和组序列 Generate job sequence in each family and family sequence randomly
		{
			JobSeqInFam = this->m_JobsInEachFamily;
			for (int fam = 0; fam < JobSeqInFam.size(); fam++)
			{
				shuffle(JobSeqInFam[fam].begin(), JobSeqInFam[fam].end(), rand_generator());//打乱组内工件顺序
			}
			for (int fam = 0; fam < FamPrmu.size(); fam++)
			{
				FamPrmu[fam] = fam;
			}
			shuffle(FamPrmu.begin(), FamPrmu.end(), rand_generator());//打乱组顺序
		}

		FacFamSeq.clear();
		FacFamSeq.resize(this->m_Factories);

		//int Makespan = 0;
		int CurFam = -1;
		int BestFac = -1, BestPos = -1;

		int Fam = 0;
		for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
		{
			CurFam = FamPrmu[Fam];
			FacFamSeq[Fac].insert(FacFamSeq[Fac].begin() + 0, CurFam);
			Fam++;
		}
		if (PS == 0 || PS == 2 || PS == 4)
		{
			for (; Fam < FamPrmu.size(); Fam++)
			{
				CurFam = FamPrmu[Fam];
				BestFac = -1;
				BestPos = -1;
				this->FindBestPosToInsertFamForAllFac_Makespan(FacFamSeq, JobSeqInFam, speed_matrix_0, CurFam, BestFac, BestPos);
				FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);

			}
		}
		else if (PS == 1 || PS == 3 || PS == 5)
		{
			//以ind为目标
			for (; Fam < FamPrmu.size(); Fam++)
			{
				CurFam = FamPrmu[Fam];
				BestFac = -1;
				BestPos = -1;
				this->FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, speed_matrix_0, CurFam, BestFac, BestPos);
				FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
			}
		}

		else
		{
			bool flag = wyt_rand(0.5);
			if (flag)
			{
				for (; Fam < FamPrmu.size(); Fam++)
				{
					CurFam = FamPrmu[Fam];
					BestFac = -1;
					BestPos = -1;
					this->FindBestPosToInsertFamForAllFac_Makespan(FacFamSeq, JobSeqInFam, speed_matrix_0, CurFam, BestFac, BestPos);
					FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
				}
			}
			else
			{
				for (; Fam < FamPrmu.size(); Fam++)
				{
					CurFam = FamPrmu[Fam];
					BestFac = -1;
					BestPos = -1;
					this->FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, speed_matrix_0, CurFam, BestFac, BestPos);
					FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
				}
			}
		}


		int MS = 0; float TEC = 0;

		GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, speed_matrix_0, FacSpan, FacEC, MS, TEC);

		//cout  << PS << ":" << MS << "," << TEC << endl;
		// 初始化档案集（赋值）
		m_RefJobSeqInFamArray[PS] = JobSeqInFam;  //工件组序列
		m_RefFacFamSeqArray[PS] = FacFamSeq;   //工厂组序列
		m_RefFacSpanArray[PS] = FacSpan;  //所有工厂的完工时间
		m_RefFacECArray[PS] = FacEC;  //所有工厂的能耗
		m_RefSpanArray[PS] = MS;  // 最大完工时间
		m_RefTECArray[PS] = TEC;  //总能耗
		m_RefSpeedVector[PS] = speed_matrix_0;


		//更新参考集
		m_CMOEAPopulation[PS].m_FacFamSeqArray = FacFamSeq;
		m_CMOEAPopulation[PS].m_JobSeqInFamArray = JobSeqInFam;
		m_CMOEAPopulation[PS].MS = MS;
		m_CMOEAPopulation[PS].TEC = TEC;
		m_CMOEAPopulation[PS].m_SpeedVector = speed_matrix_0;  //速度矩阵
	}

	// 初始化组种群 Initilize Familiy-sequence population, i.e., PS1
	for (int PS = 0; PS < m_PS1; PS++)
	{
		//将档案集的解赋给组序列
		m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[PS];
		m_SpanArray1[PS] = m_RefSpanArray[PS];
		m_TECArray1[PS] = m_RefTECArray[PS];
		m_Map1[PS] = PS;  //标记档案集中对应的jobpop序号
	}

	// 初始化工件种群 Initialize Job-Sequence population, i.e., PS2
	for (int PS = 0; PS < m_PS2; PS++)
	{
		//前AS个解
		m_JobSeqInFamArray[PS] = m_RefJobSeqInFamArray[PS];
		m_SpanArray2[PS] = m_RefSpanArray[PS];
		m_TECArray2[PS] = m_RefTECArray[PS];
		m_Map2[PS] = PS;  //标记对应的档案集中fampop序号
	}
	//归一化
	Normalize(m_CMOEAPopulation, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC);//种群，最低点MS，最低点TEC，理想点MS，理想点TEC
	//cout << "初始化结束！" << endl;

}

/*************************************解的更新*************************************************************/
void QCMOEA_PM::UpdateArchiveGroupJobSet(int Muu)
{
	//cout << "解更新开始！" << endl;

	//占优关系
	Individual_PM::Pareto_relation(m_CMOEAPopulation);

	//弱化支配解 de-emphasize dominated solutions
	for (int j = 0; j < m_CMOEAPopulation.size(); j++)
	{
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				if (m_CMOEAPopulation[i].pareto_rel[j] == 1)//j占优i
				{
					//cout << j << "占优" << i << endl;
					m_CMOEAPopulation[i].flag = 999;//flag+ 弱化支配解
				}
			}
		}
	}

	//cout << "弱化支配解计算完毕！" << endl;
	vector<Individual_PM> Temp;
	Temp.clear();
	//对参考集(NR)进行预处理，根据占优规则选择NR中的非支配解。
	for (int i = 0; i < m_CMOEAPopulation.size(); i++)
	{

		//根据占优规则选择NR中的非支配解
		if (m_CMOEAPopulation[i].flag == 0)
		{
			Temp.push_back(m_CMOEAPopulation[i]);
		}
	}
	//cout << "非支配解计算完毕！" << endl;
	//cout << Temp.size() << endl;

	//如果选择的个体数量小于mu×AN ，则继续选择，直到达到mu×AN个个体。这些mu×AN解形成一个新的NR来替换原来的NR。
	//参数：mu(Archive种群倍数)
	int mu = Muu;
	if (Temp.size() < mu * m_RefSize)
	{
		while (true)
		{
			for (int i = 0; i < m_CMOEAPopulation.size(); i++)
			{
				if (m_CMOEAPopulation[i].flag == 999)
				{
					Temp.push_back(m_CMOEAPopulation[i]);
					if (Temp.size() == 10 * m_RefSize)
						break;
				}
			}
			break;
		}
	}
	//cout << "mu×AN解决方案构成完毕！" << endl;

	//这些mu×AN解决方案形成一个新的NR，以取代原来的NR，并重新计算其收敛和分散指标。
	m_CMOEAPopulation.clear();
	for (int i = 0; i < Temp.size(); i++)
		m_CMOEAPopulation.push_back(Temp[i]);

	//归一化
	Normalize(m_CMOEAPopulation, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC);
	//cout << "nadirpointMS：" << m_nadirpointMS << "\tnadirpointTEC：" << m_nadirpointTEC << endl;

	//计算收敛指标
	Individual_PM::calc_dpconvergence_ind(m_CMOEAPopulation);
	//cout << "计算收敛指标后11111111111111" << endl;

	//计算分散指标
	calc_distribution_ind(m_CMOEAPopulation);
	//cout << "计算分布指标后11111111111111" << endl;
	int n, n1, n2, n3, nrank;
	int the_one;
	n = 0;   //NRs.size  预选择的个体集合
	n1 = 0;  //Q.size   参考集的个体数量
	n2 = 0;  //NRe.size 由distribution threshold弱化的集合
	n3 = 0;  //NRd.size  被支配个体的集合

	/*n1参考集 mu*AN个 从n1里面选取一个收敛指标最小的放进NRs，并在n1里面移除它，然后根据分布性指标选取离刚才收敛性指标最小的那个解距离比较近的（距离<thr_zeta）放进NRe，
	然后根据占优关系选取占有于一些解，那些解放进NRd，一直到NRs中数量=AN*/
	vector<Individual_PM> TempCMOEAPopulation;
	TempCMOEAPopulation.clear();

	//占优关系
	Pareto_relation(m_CMOEAPopulation);
	//cout << "占优关系后11111111111111" << endl;

	n1 = m_CMOEAPopulation.size();

	while (n1 > 0)
	{
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				the_one = i;
				break;
			}
		}

		m_CMOEAPopulation[the_one].flag = 1;
		TempCMOEAPopulation.push_back(m_CMOEAPopulation[the_one]);
		n++;


		//弱化附近解 de-emphasize neighbors
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				if (m_CMOEAPopulation[i].distribution_ind[the_one] < thr_zeta)
				{
					m_CMOEAPopulation[i].flag = 999;
					n2++;
				}
			}
		}

		//弱化支配解 de-emphasize dominated solutions
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				if (m_CMOEAPopulation[i].pareto_rel[the_one] == 1)//the_one占优于i
				{
					//不好的flag设为999
					m_CMOEAPopulation[i].flag = 999;
					n3++;
				}
			}
		}

		//剩余数量 number of the rest
		n1 = 0;
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				n1++;
			}
		}
	}
	//cout << "参数前11111111111111" << endl;

	//参数 threshold
	float oldthr_zeta = thr_zeta;
	float ratio = n * 1.0 / (m_RefSize * 1.0);
	if (n3 < (mu - 1) * m_RefSize)
		thr_zeta = thr_zeta * exp((ratio - 1.0) / (2 * 1.0));
	else
		thr_zeta = oldthr_zeta;


	while (n < m_RefSize)
	{
		//剩余数量 number of the rest
		n1 = 0;
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				n1++;
			}
		}

		if (n1 == 0)
		{
			for (int i = 0; i < m_CMOEAPopulation.size(); i++)
			{
				if (m_CMOEAPopulation[i].flag == 999)
				{
					m_CMOEAPopulation[i].flag = 0;
				}
			}
		}

		//选择一个
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				the_one = i;
				break;
			}
		}
		m_CMOEAPopulation[the_one].flag = 1;
		TempCMOEAPopulation.push_back(m_CMOEAPopulation[the_one]);
		n++;

		//弱化附近解 de-emphasize neighbors
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				if (m_CMOEAPopulation[i].distribution_ind[the_one] < oldthr_zeta)
				{
					m_CMOEAPopulation[i].flag = 999;
					n2++;
				}
			}
		}
		//cout << "11111111111111" << endl;
		//弱化支配解 de-emphasize dominated solutions
		for (int i = 0; i < m_CMOEAPopulation.size(); i++)
		{
			if (m_CMOEAPopulation[i].flag == 0)
			{
				if (m_CMOEAPopulation[i].pareto_rel[the_one] == 1)
				{
					m_CMOEAPopulation[i].flag = 999;
					n3++;
				}
			}
		}
	}

	vector<Individual_PM> Temp1;
	Temp1.clear();
	//cout << "TempCMOEAPopulation:" << TempCMOEAPopulation.size() << endl;

	for (int i = 0; i < TempCMOEAPopulation.size(); i++)
	{
		Temp1.push_back(TempCMOEAPopulation[i]);
	}
	//cout << "CMOEAPopulation:" << CMOEAPopulation.size() << endl;

	//变速
	Speed_mutation(TempCMOEAPopulation, Temp1, m_CMOEAPopulation);
	//cout <<"m_RefSize:"<< m_RefSize << endl;

	for (int PS = 0; PS < m_RefSize; PS++)
	{
		// 更新档案集（赋值） 
		m_RefJobSeqInFamArray[PS] = TempCMOEAPopulation[PS].m_JobSeqInFamArray;  //工件组序列
		m_RefFacFamSeqArray[PS] = TempCMOEAPopulation[PS].m_FacFamSeqArray;   //工厂组序列		
		m_RefSpanArray[PS] = TempCMOEAPopulation[PS].MS;  // 最大完工时间
		m_RefTECArray[PS] = TempCMOEAPopulation[PS].TEC;  //总能耗
		m_RefSpeedVector[PS] = TempCMOEAPopulation[PS].m_SpeedVector;
		m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[PS];
		m_JobSeqInFamArray[PS] = m_RefJobSeqInFamArray[PS];
		m_SpeedVector = m_RefSpeedVector[PS];
	}
}

void QCMOEA_PM::SaveTECandDeMS(vector<Individual_PM>& CMOEAPopulation)
{
	vector<Individual_PM> tempCMOEAPopulation;
	tempCMOEAPopulation.clear();
	vector<Individual_PM> tempCMOEAPopulation2;
	tempCMOEAPopulation2.clear();
	vector<Individual_PM> tempCMOEAPopulation3;
	tempCMOEAPopulation3.clear();
	vector<Individual_PM> tempCMOEAPopulation4;
	tempCMOEAPopulation4.clear();


	int orgSize = CMOEAPopulation.size();
	cout << orgSize << endl;
	for (int PS = 0; PS < orgSize; PS++)
	{
		tempCMOEAPopulation.push_back(CMOEAPopulation[PS]);
		bool flag2 = false;
		vector<int> FacSpan(m_Factories);//组，工厂完工时间
		//节能策略1 减能耗
		vector<vector<int>> DelayTime;
		DelayTime.clear();
		DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
		GetDelayTimeAllFac(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, tempCMOEAPopulation[PS].m_SpeedVector, FacSpan, DelayTime);
		vector<vector<int>> JCTime(this->m_Jobs, vector<int>(this->m_Machines, 0));
		for (int Fac = 0; Fac < m_Factories; Fac++)
		{
			GetCTimeForPerFac(tempCMOEAPopulation[PS].m_FacFamSeqArray[Fac], tempCMOEAPopulation[PS].m_JobSeqInFamArray, tempCMOEAPopulation[PS].m_SpeedVector);
		}

		bool sign = false;

		for (int j = 0; j < m_Jobs; j++)
		{
			for (int i = 1; i < m_Machines; i++)
			{
				if ((DelayTime[j][i] > 0) && (i < m_Machines - 1))//中间机器
				{
					if ((JCTime[j][i] < (JCTime[j][i + 1] - m_TureJobOpertime[j][i + 1][tempCMOEAPopulation[PS].m_SpeedVector[j][i + 1]])))
					{
						int Speedlevel = tempCMOEAPopulation[PS].m_SpeedVector[j][i] - 1;

						for (int level = Speedlevel; level > 0; level--)
						{
							if ((m_TureJobOpertime[j][i][level] - m_TureJobOpertime[j][i][tempCMOEAPopulation[PS].m_SpeedVector[j][i]]) < DelayTime[j][i])
							{
								//如果速度慢一档后的处理时间能够弥补DelayTime， 就减速（减能耗）
								if ((JCTime[j][i] + (m_TureJobOpertime[j][i][level] - m_TureJobOpertime[j][i][tempCMOEAPopulation[PS].m_SpeedVector[j][i]])) < (JCTime[j][i + 1] - m_TureJobOpertime[j][i + 1][tempCMOEAPopulation[PS].m_SpeedVector[j][i + 1]]))
								{
									tempCMOEAPopulation[PS].m_SpeedVector[j][i] = level;
									sign = true;
								}
								else
									break;
							}
							else
								break;
						}
					}

				}

				else if ((DelayTime[j][i] > 0) && (i == m_Machines - 1))
				{
					int Speedlevel = tempCMOEAPopulation[PS].m_SpeedVector[j][i] - 1;
					for (int level = Speedlevel; level > 0; level--)
					{
						if ((m_TureJobOpertime[j][i][level] - m_TureJobOpertime[j][i][tempCMOEAPopulation[PS].m_SpeedVector[j][i]]) < DelayTime[j][i])
						{
							tempCMOEAPopulation[PS].m_SpeedVector[j][i] = level;
							sign = true;
						}
						else
							break;
					}
				}
			}//end machine
		}//end job

		if (sign == true)
		{
			//如果第一个策略执行 就重新计算makespan和tec
			int Makespan1 = 0;
			float TotalEC1 = 0.0;

			Makespan1 = GetCTimeForAllFac(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, tempCMOEAPopulation[PS].m_SpeedVector);
			TotalEC1 = GetTECForAllFac(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, tempCMOEAPopulation[PS].m_SpeedVector);
			tempCMOEAPopulation[PS].MS = Makespan1;
			tempCMOEAPopulation[PS].TEC = TotalEC1;

			//判断
			if ((CMOEAPopulation[PS].MS >= tempCMOEAPopulation[PS].MS) && (CMOEAPopulation[PS].TEC > tempCMOEAPopulation[PS].TEC))
			{
				CMOEAPopulation.push_back(tempCMOEAPopulation[PS]);
				flag2 = true;
			}
		}


		//节能策略2：在1的基础上，在不超过Makespan的情况下，降低非关键工厂的v，来降低TEC
		vector<int> FacSpan3(m_Factories);
		if (flag2)//如果执行了1，并且得到的解非支配
		{
			tempCMOEAPopulation3.push_back(tempCMOEAPopulation[PS]);
			vector<int> notcirfac;
			notcirfac.clear();
			for (int fac = 0; fac < m_Factories; fac++)
			{
				if (tempCMOEAPopulation3.back().MS > FacSpan[fac])
				{
					notcirfac.push_back(fac);
				}
				FacSpan3[fac] = FacSpan[fac];
			}
			int teration = 0;//迭代
			float tempTEC = 0.0;

			do {
				teration++;
				if (teration == 1)
					tempTEC = tempCMOEAPopulation[PS].TEC;

				if (teration == 2)
					tempTEC = tempCMOEAPopulation3.back().TEC;
				if (teration == 3)
					break;

				bool judgeV = true;
				vector<vector<int>> tempSpeedMatrix;
				tempSpeedMatrix.clear();
				tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
				tempSpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;


				//对非关键工厂进行变速
				for (int f = 0; f < notcirfac.size(); f++)
				{
					judgeV = true;
					//int speedLevel;
					for (int fam = 0; fam < tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]].size(); fam++)
					{
						int CurFam = tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]][fam];
						for (int job = 0; job < tempCMOEAPopulation3.back().m_JobSeqInFamArray[CurFam].size(); job++)
						{
							int CurJob = tempCMOEAPopulation3.back().m_JobSeqInFamArray[CurFam][job];
							for (int m = 0; m < m_Machines; m++)
							{
								if (tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] > 0)
								{
									tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] = tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] - 1;
									judgeV = false;

									//真实的处理时间
								}

							}
						}
					}
					//若有速度有变化
					if (judgeV == false)
					{
						//得到工厂f变速后的Span
						int tempnotcirSpan = GetCTimeForPerFac(tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]], tempCMOEAPopulation3.back().m_JobSeqInFamArray, tempCMOEAPopulation3.back().m_SpeedVector);//得到每个工厂的完工时间
						//如果改变速度后未超过Makspan，才真正变速
						if (tempnotcirSpan > tempCMOEAPopulation3.back().MS)
						{
							tempCMOEAPopulation3.back().m_SpeedVector = tempSpeedMatrix;
						}
						else
						{
							tempSpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;
						}

					}

				}


				int Makespan3 = 0;
				float TotalEC3 = 0.0;

				Makespan3 = GetCTimeForAllFac(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, tempCMOEAPopulation3.back().m_SpeedVector);
				TotalEC3 = GetTECForAllFac(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, tempCMOEAPopulation3.back().m_SpeedVector);

				tempCMOEAPopulation3.back().MS = Makespan3;
				tempCMOEAPopulation3.back().TEC = TotalEC3;

				//判断
				if ((CMOEAPopulation[PS].MS >= tempCMOEAPopulation3.back().MS) && (tempCMOEAPopulation[PS].TEC > tempCMOEAPopulation3.back().TEC))
				{
					CMOEAPopulation.push_back(tempCMOEAPopulation3.back());
				}

			} while (tempTEC > tempCMOEAPopulation3.back().TEC);

		}

		//节能策略3：在2的基础上，在不超过原TEC的情况下，提高关键工厂的v，来降低MS
		vector<int> FacSpan4(m_Factories);
		if (flag2)
		{
			tempCMOEAPopulation4.push_back(tempCMOEAPopulation3.back());
			vector<int> cirfac;
			cirfac.clear();
			for (int fac = 0; fac < m_Factories; fac++)
			{
				if (tempCMOEAPopulation4.back().MS == FacSpan[fac])
				{
					cirfac.push_back(fac);
				}
				FacSpan4[fac] = FacSpan[fac];
			}
			int teration = 0;
			int tempMS = 0;

			do {
				teration++;
				if (teration == 1)
					tempMS = CMOEAPopulation[PS].MS;

				if (teration == 2)
				{
					tempMS = tempCMOEAPopulation4.back().MS;
					//cirfac.clear();
					vector<int> tempcirfac;
					tempcirfac.clear();
					//判断关键工厂是否改变
					for (int fac = 0; fac < m_Factories; fac++)
					{
						if (tempCMOEAPopulation4.back().MS == FacSpan4[fac])
						{
							tempcirfac.push_back(fac);
						}
					}
					bool judgecirfacifchange = true;
					for (int i = 0; i < tempcirfac.size(); i++)
					{
						if (tempcirfac[i] != cirfac[i])
						{
							judgecirfacifchange = false;
							break;
						}
					}
					if (!judgecirfacifchange)
					{
						teration = 1;
						cirfac.clear();
						cirfac = tempcirfac;
					}

				}

				if (teration == 3)
					break;

				bool judgeV = true;
				vector<vector<int>> tempSpeedMatrix;
				tempSpeedMatrix.clear();
				tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
				tempSpeedMatrix = tempCMOEAPopulation4.back().m_SpeedVector;


				//对非关键工厂进行变速
				for (int f = 0; f < cirfac.size(); f++)
				{
					judgeV = true;
					for (int fam = 0; fam < tempCMOEAPopulation4.back().m_FacFamSeqArray[cirfac[f]].size(); fam++)
					{
						int CurFam = tempCMOEAPopulation4.back().m_FacFamSeqArray[cirfac[f]][fam];
						for (int job = 0; job < tempCMOEAPopulation4.back().m_JobSeqInFamArray[CurFam].size(); job++)
						{
							int CurJob = tempCMOEAPopulation4.back().m_JobSeqInFamArray[CurFam][job];
							for (int m = 0; m < m_Machines; m++)
							{
								if (tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] < 2)
								{
									tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] = tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] + 1;
									judgeV = false;
								}

							}
						}
					}
					//若有速度有变化
					if (judgeV == false)
					{
						//得到工厂f变速后的TEC
						GetCTimeForAllFac(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, tempCMOEAPopulation4.back().m_SpeedVector);
						float tempTotalEC = GetTECForAllFac(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, tempCMOEAPopulation4.back().m_SpeedVector);


						//如果改变速度后未超过原TEC，才真正变速
						if (tempTotalEC > CMOEAPopulation[PS].TEC)
						{
							tempCMOEAPopulation4.back().m_SpeedVector = tempSpeedMatrix;
							break;
						}

					}

				}

				int Makespan4 = 0;
				float TotalEC4 = 0.0;

				vector<vector<int>> TempMC;
				Makespan4 = GetCTimeForAllFac(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, tempCMOEAPopulation4.back().m_SpeedVector);
				TotalEC4 = GetTECForAllFac(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, tempCMOEAPopulation4.back().m_SpeedVector);
				tempCMOEAPopulation4.back().MS = Makespan4;
				tempCMOEAPopulation4.back().TEC = TotalEC4;

				//判断
				if ((CMOEAPopulation[PS].MS > tempCMOEAPopulation4.back().MS) && (CMOEAPopulation[PS].TEC >= tempCMOEAPopulation4.back().TEC))
				{
					CMOEAPopulation.push_back(tempCMOEAPopulation4.back());
				}

			} while (tempMS > tempCMOEAPopulation4.back().MS);

		}
	}

}

void QCMOEA_PM::EvolutionProcess(int mu)
{
	thr_zeta = 1.0;
	InitialPop(); //初始化种群
	long InitTime = Base::GetElapsedProcessTime();
	m_InitTime = InitTime;
	int count = 0;	
	vector<vector<vector<double>>> Qfam(m_PS1, vector<vector<double>>(3, vector < double>(5)));//q表
	while (::GetElapsedProcessTime() - InitTime < m_TimeLimit)
	{
		//协同进化--组
		for (int PS = 0; PS < m_PS1; PS++)
		{
			vector<vector<int>> FacFamSeq = m_FacFamSeqArray[PS];
			int Map1 = wyt_rand(m_RefSize);// 从档案集中随机挑选一个合作者（工件序列） form a solution by randomly selecting a RefJobSeq;
			vector<int> FacSpan(m_Factories, 0);
			vector<float> FacTEC(m_Factories, 0);
			int ObjectMS;
			float ObjectTEC;
			GetMSandTECForPerandToalFac(FacFamSeq, m_RefJobSeqInFamArray[Map1], m_RefSpeedVector[PS], FacSpan, FacTEC, ObjectMS, ObjectTEC);
			int newstate = 0;
			int oldstate = 0;
			int OrgMS = ObjectMS;
			float OrgTEC = ObjectTEC;
			//使用q-learning进行组序列选择策略
			int act = QIG_DetermineAction(oldstate, Qfam[PS]);
			QIG_PerformAction_Fams(act, FacFamSeq, m_RefJobSeqInFamArray[Map1], m_RefSpeedVector[PS], FacSpan, FacTEC, ObjectMS, ObjectTEC);
			double Reward = QIG_CaculateReward(OrgMS, ObjectMS, OrgTEC, ObjectTEC);
			newstate = QIG_DetermineState(OrgMS, ObjectMS, OrgTEC, ObjectTEC);
			QIG_UpdateQ(oldstate, newstate, act, Reward, Qfam[PS]);
			oldstate = newstate;
		}
		//协同进化--工件
		for (int PS = 0; PS < m_PS2; PS++)
		{
			vector<vector<int>> JobSeqInFam = m_JobSeqInFamArray[PS];
			int Map2 = wyt_rand(m_RefSize);// 从档案集中随机挑选一个合作者（组序列） form a solution by randomly selecting a ReFacFamseq
			vector<int> FacSpan(m_Factories, 0);
			vector<float> FacEC(m_Factories, 0);
			int ObjectMS;
			float ObjectTEC;
			GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[Map2], JobSeqInFam, m_RefSpeedVector[PS], FacSpan, FacEC, ObjectMS, ObjectTEC);
			BasedindCirJobInFamTobestPos(m_RefFacFamSeqArray[Map2], JobSeqInFam, m_RefSpeedVector[PS], m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
			BasedindSwapJob(m_RefFacFamSeqArray[Map2], JobSeqInFam, m_RefSpeedVector[PS], m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		}
		//协同进化--档案集
		for (int PS = 0; PS < m_RefSize; PS++)
		{
			BasedindRefInsert(m_RefFacFamSeqArray[PS], m_RefJobSeqInFamArray[PS], m_RefSpeedVector[PS], m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
				m_CMOEAPopulation, m_RefSpanArray[PS], m_RefTECArray[PS]);
		}
		UpdateArchiveGroupJobSet(mu);
		count++;

	}//END While

	SaveTECandDeMS(m_CMOEAPopulation);
}

int QCMOEA_PM::RunEvolution(string Dir, int CPUTime, int Reps, vector<vector<Individual_PM>>& CMOEAFinalAfterRepParetoSet, int AN, int mu, int Instances)// 10, 4
{

	int FamD = 6;
	this->ReadInstanceFileNameList(Dir);

	ofstream ofile;
	ofile.open("..\\Result\\Model_Result_QCMOEA_PM.csv");
	if (!ofile.is_open())
	{
		cout << "..\\Result\\Model_Result_QCMOEA_PM.csv" << "\t is not open" << endl;
		exit(0);
	}
	for (int ins = 0; ins < Instances; ins++)
	{
		this->ReadInstance(ins);
		this->GetJobTotalPTime();
		this->GetFamTotalPTime();
		this->GetFamAvgSetupTime();

		vector<Individual_PM> FinalCMOEAParetoSet;//最终的paretoset
		vector<Individual_PM> AfterRepParetoSet;
		vector<Individual_PM> TempAfterRepParetoSet;
		TempAfterRepParetoSet.clear();

		for (int r = 0; r < Reps; r++)
		{
			long TimeLimit = CPUTime * m_Machines * m_Families; //original: 20 * m_Jobs * m_Machines

			this->SetParameters(AN, FamD, TimeLimit, m_AllJobTotalPTime);
			long StartTime_IG = ::GetElapsedProcessTime();
			//进化过程
			EvolutionProcess(mu);

			//非支配解
			m_CMOEAParetoSet.clear();
			Individual_PM::Pareto_relation(m_CMOEAPopulation);

			//弱化支配解 de-emphasize dominated solutions
			for (int j = 0; j < m_CMOEAPopulation.size(); j++)
			{
				for (int i = 0; i < m_CMOEAPopulation.size(); i++)
				{
					if (m_CMOEAPopulation[i].flag == 0)
					{
						if (m_CMOEAPopulation[i].pareto_rel[j] == 1)
						{
							m_CMOEAPopulation[i].flag = 999;
						}
					}
				}
			}

			for (int i = 0; i < m_CMOEAPopulation.size(); i++)
			{
				if (m_CMOEAPopulation[i].flag == 0)
					m_CMOEAParetoSet.push_back(m_CMOEAPopulation[i]);
			}

			//去除重复
			FinalCMOEAParetoSet.clear();
			for (int i = 0; i < m_CMOEAParetoSet.size(); i++)
			{
				bool fg = true;
				for (int j = 0; j < FinalCMOEAParetoSet.size(); j++)
				{
					if ((FinalCMOEAParetoSet[j].MS == m_CMOEAParetoSet[i].MS) && (FinalCMOEAParetoSet[j].TEC == m_CMOEAParetoSet[i].TEC))
					{
						fg = false;
						break;
					}
				}
				if (fg)
				{
					FinalCMOEAParetoSet.push_back(m_CMOEAParetoSet[i]);
				}

			}

			long EndTime_IG = ::GetElapsedProcessTime();

			string Str = this->m_InstanceNameList[ins];

			cout << Str << "\t" << this->m_Factories << "\t" << this->m_Machines << "\t" << this->m_Families << "\t"
				<< this->m_SetupType << "\t" << this->m_Jobs << "\t";


			for (int PS = 0; PS < FinalCMOEAParetoSet.size(); PS++)
			{
				cout << FinalCMOEAParetoSet[PS].MS << "\t" << FinalCMOEAParetoSet[PS].TEC << "\t";

				TempAfterRepParetoSet.push_back(FinalCMOEAParetoSet[PS]);
			}
			cout << r + 1 << "\t" << EndTime_IG - StartTime_IG << endl;


		}//end rep

		//非支配解
		AfterRepParetoSet.clear();

		Individual_PM::Pareto_relation(TempAfterRepParetoSet);

		//弱化支配解 de-emphasize dominated solutions
		for (int j = 0; j < TempAfterRepParetoSet.size(); j++)
		{
			for (int i = 0; i < TempAfterRepParetoSet.size(); i++)
			{
				if (TempAfterRepParetoSet[i].flag == 0)
				{
					if (TempAfterRepParetoSet[i].pareto_rel[j] == 1)
					{
						TempAfterRepParetoSet[i].flag = 999;
					}
				}
			}
		}

		for (int i = 0; i < TempAfterRepParetoSet.size(); i++)
		{
			if (TempAfterRepParetoSet[i].flag == 0)
				AfterRepParetoSet.push_back(TempAfterRepParetoSet[i]);
		}

		//去除重复
		vector<Individual_PM> FinalAfterRepParetoSet;
		FinalAfterRepParetoSet.clear();
		for (int i = 0; i < AfterRepParetoSet.size(); i++)
		{
			bool fg = true;
			for (int j = 0; j < FinalAfterRepParetoSet.size(); j++)
			{
				if ((FinalAfterRepParetoSet[j].MS == AfterRepParetoSet[i].MS) && (FinalAfterRepParetoSet[j].TEC == AfterRepParetoSet[i].TEC))
				{
					fg = false;
					break;
				}
			}
			if (fg)
			{
				FinalAfterRepParetoSet.push_back(AfterRepParetoSet[i]);
			}
		}

		string Str = this->m_InstanceNameList[ins];

		cout << "AfterRep:" << "\t" << Str << "\t" << this->m_Factories << "\t" << this->m_Machines << "\t" << this->m_Families << "\t"
			<< this->m_SetupType << "\t" << this->m_Jobs << "\t";
		cout << endl;

		CMOEAFinalAfterRepParetoSet[ins].clear();
		for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++)
		{
			cout << FinalAfterRepParetoSet[PS].MS << "\t" << FinalAfterRepParetoSet[PS].TEC << "\t";

			ofile << FinalAfterRepParetoSet[PS].MS << "," << FinalAfterRepParetoSet[PS].TEC << ",";

			CMOEAFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);

		}
		ofile << endl;
		cout << endl;

	}//end ins
	ofile.close();
	return 0;
}


void QCMOEA_PM::QIG_PerformAction_Fams(int act, vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	switch (act) {
	case 0:
		//cout << "0" << endl;
		BasedindRandFamInFacTobestPos(FacFamSeq, JobSeqInFam, SpeedMatrix, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
			m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		break;
	case 1:
		//cout << "1" << endl;

		BasedindSwapFam(FacFamSeq, JobSeqInFam, SpeedMatrix, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
			m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		break;
	case 2:
		//cout << "2" << endl;

		BasedindLS_SetupInsert(FacFamSeq, JobSeqInFam, SpeedMatrix, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
			m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		break;
	case 3:
		//cout << "3" << endl;

		BasedindLS_SetupSwap(FacFamSeq, JobSeqInFam, SpeedMatrix, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
			m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		break;
	case 4:
		//cout << "4" << endl;

		Basedind_DCFams(FacFamSeq, JobSeqInFam, SpeedMatrix, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC,
			m_CMOEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
		break;

	default:
		std::cerr << "Error: Invalid action " << act << " in QIG_PerformAction_Fams." << std::endl;
	}

	// 重新计算总的 MS 和 TEC 确保最新解
	GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, SpeedMatrix, FacSpan, FacEC, ObjectMS, ObjectTEC);

}

int QCMOEA_PM::QIG_DetermineAction(int state, const vector<vector<double>>& Q) {
	double rande = wyt_rand_include_right(0.0, 1.0);
	if (rande < e)
	{
		vector<int> QI = findMaxIndices(Q[state]);
		return int(QI[wyt_rand(QI.size())]);
	}
	else
	{
		return wyt_rand(Q[state].size());
	}
}

int QCMOEA_PM::QIG_DetermineState(int PreMS, int actMS, float PreTEC, float actTEC)
{
	if (actMS < PreMS && actTEC < PreTEC)
	{
		return 0;//好
	}
	else if (actTEC < PreTEC || actMS < PreMS)
	{
		return 1;//改进一个解 较好
	}
	else
	{
		return 2;//没改进
	}
}

void QCMOEA_PM::QIG_UpdateQ(int oldstate, int newstate, int act, double reward, vector<vector<double>>& Q)
{
	//double e = 0.75;贪婪策略  discountFactor = 0.2;折扣因子 rewardFactor = 0.75;学习率
	if (newstate == 0 || newstate == 1)
	{
		Q[oldstate][act] = (1.0 - rewardFactor) * Q[oldstate][act] + rewardFactor * (reward + discountFactor * *max_element(Q[newstate].begin(), Q[newstate].end()));
	}
}

double QCMOEA_PM::QIG_CaculateReward(int preMS, int newMS, float preTEC, float newTEC)
{
	double reward = (preMS - newMS) * 10 + (preTEC - newTEC) * 5;
	return reward;
}


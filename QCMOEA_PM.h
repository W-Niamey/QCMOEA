#ifndef QCMOEA_PM_H
#define QCMOEA_PM_H
#include <unordered_map>
#include "NOperator_PM.h"
#include "Individual_PM.h"

class QCMOEA_PM : public NOperator_PM, public Individual_PM
{
public:
	QCMOEA_PM();
	virtual ~QCMOEA_PM();
	int RunEvolution(string Dir, int CPUTime, int Reps, vector<vector<Individual_PM>>& CMOEAFinalAfterRepParetoSet, int AN, int mu, int Instances);

private:
	long m_InitTime;
	long m_TimeLimit;

	int m_Temperature;
	int d;//破坏
	double e = 0.75; //贪婪策略
	double discountFactor = 0.2; //折扣因子
	double rewardFactor = 0.75;//学习率

	//asdf
	int m_RefSize;
	int m_Popsize;

	float thr_zeta;

	int m_nadirpointMS;    //最低点
	float m_nadirpointTEC;

	int m_idealpointMS;	  //理想点
	float m_idealpointTEC;

	int m_FamD;  // 随机选择

	vector<int> m_OperType1;
	vector<int> m_OperType2;
	//档案集
	vector<vector<vector<int>>> m_RefSpeedVector;
	vector<vector<vector<int>>> m_RefFacFamSeqArray;
	vector<vector<vector<int>>> m_RefJobSeqInFamArray;
	vector<vector<int>> m_RefFacSpanArray;
	vector<vector<float>> m_RefFacECArray;
	vector<int> m_RefSpanArray;
	vector<float> m_RefTECArray;
	vector <int> m_nRefCriFacArray;

	vector<Individual_PM> m_CMOEAPopulation;
	vector<Individual_PM> m_CMOEAParetoSet;  //ParetoSet解集

	vector<vector<vector<int>>> m_FacFamSeqArray;
	vector<vector<vector<int>>> m_JobSeqInFamArray;
	vector<int> m_SpanArray1;
	vector<int> m_SpanArray2;
	vector<float> m_TECArray1;
	vector<float> m_TECArray2;

	int m_DesLen;
	int m_PS1;
	int m_PS2;
	vector<int> m_Map1;
	vector<int> m_Map2;

	vector<int> m_Age1;
	vector<int> m_Age2;
	int m_AgeLimit;
	double m_Rate_FamOp;
	vector<int> m_RecordSpan;

	void SetParameters(int AN,int NumberOfExtractedFams, long TimeLimit, int AllJobTotalPTime);

	void InitialPop();

	void EvolutionProcess(int mu);

	double GetTemperature();

	/*********************组进化************************/
	void BasedindRandFamInFacTobestPos(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacTEC, int& ObjectMS, float& ObjectTEC);

	void BasedindSwapFam(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	void BasedindLS_SetupInsert(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	void BasedindLS_SetupSwap(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	void Basedind_DCFams(vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	/**********************工件进化********************************/
	void BasedindCirJobInFamTobestPos(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC,
		int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	void BasedindSwapJob(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);


	/*******************************档案集进化************************************/
	void BasedindRefInsert(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

	void BasedindDR_FamsAndJobs(vector<Individual_PM>& CMOEAPopulation, vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FamsExtracted, unordered_map<int, vector<int>>& JobsExtracted, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int m_FamD);

	void UpdateArchiveGroupJobSet(int mu);

	void SaveTECandDeMS(vector<Individual_PM>& CMOEAPopulation);

	int QIG_DetermineAction(int state, const vector<vector<double>>& Q);

	int QIG_DetermineState(int MS, int actMS, float TEC, float actTEC);

	void QIG_UpdateQ(int oldstate, int newstate, int act, double reward, vector<vector<double>>& Q);

	double QIG_CaculateReward(int MS, int actMS, float TEC, float actTEC);

	void QIG_PerformAction_Fams(int act, vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC);

};

#endif

#ifndef NOPERATOR_PM_H
#define NOPERATOR_PM_H

#include "Individual_PM.h"
#include "Problem_PM.h"

class NOperator_PM : public Problem_PM
{
public:
	NOperator_PM();

	virtual ~NOperator_PM();

	void SortJobsInFam(int SortMethod, vector<vector<int>>& JobSeqInFam);
	void SortFam(int SortMethod, vector<int>& FamPermu);
	vector<int> findMaxIndices(const vector<double>& vec);
	void Speed_mutation(vector<Individual_PM>& CMOEAPopulation, vector<Individual_PM>& Temp, vector<Individual_PM>& tureCMOEAPopulation);

	/**********************************************************************************************************************/
/********************************  ��������깤ʱ��FacCTime  **************************************/
	int GetCTimeForAllFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix);

	int GetCTimeForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix);

	/***********************************************************************************************************************/
/********************************  �����ܺ�TEC  **************************************/
	float GetTECForAllFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix);

	float GetTECForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix);

	/***********************************************************************************************************************/
/********************************  �������ʱ�䣨����ٶȲ���ʱ����ʹ�ã�  **************************************/
	int GetDelayTimeAllFac(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FacSpan, vector<vector<int>>& DelayTime);//Forward pass calculation

	int GetDelayTimePerFac(vector <int> FamSeqInFac, vector <vector <int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<vector<int>>& DelayTime);

	/***********************************************************************************************************************/
/********************************  ����ؼ����������ܺ�  **************************************/
	void GetMSandTECForPerandToalFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
		vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC);

	pair<int, float> GetCTimeAndTECForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix);

	/************************************************************************************************************************/
/****************************  ����makespan����Ѱ��Ѱ�����λ�ò���  **************************************/

	int FindBestPosToInsertFamForPerFac_Makespan(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos);

	int FindBestPosToInsertFamForAllFac_Makespan(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos);

	/************************************************************************************************************************/
/****************************  ����TEC����Ѱ��Ѱ�����λ�ò���  **************************************/

	float FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos);

	float FindBestPosToInsertFamForPerFac_TEC(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos);

	/************************************************************************************************************************/
/****************************  ��������ָ��Թ���InsJob�ڹ�����JobSeqInFam[Fam]��Ѱ�����λ��  **************************************/

	float FindBestPosToInsertFamForAllFac_Ind(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacTEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation);

	float FindBestPosToInsertFamForPerFac_Ind(int Fac, const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FacSpan, vector<float>& FacTEC, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation);

	float GetindForPerFacAfterInsertFam(int Fac, const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation);

	float FindBestPosToInsertFamForAllFac_Ind_DR(vector<Individual_PM>& CMOEAPopulation, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC);

	float FindBestPosToInsertFamForPerFac_Ind_DR(vector<Individual_PM>& CMOEAPopulation, int Fac, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC);

	float GetindForPerFacAfterInsertFam_DR(vector<Individual_PM>& CMOEAPopulation, int Fac, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
		const vector<vector<int>>& SpeedMatrix, int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC);

	/************************************************************************************************************************/
/****************************  ��������ָ��Թ���InsJob�ڹ�����JobSeqInFam[Fam]��Ѱ�����λ��  **************************************/

	void Basedind_JobsInFam(int fac, int Fam, const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC);

	float GetindForPerFacAfterInsertJob(int fac, const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int InsJob, int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacEC, int& OrgMS, float& OrgTEC, vector<Individual_PM>& CMOEAPopulation);

	float GetindForPerFacAfterInsertJob_DR(vector<Individual_PM>& CMOEAPopulation, int fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
		const vector<vector<int>>& SpeedMatrix, int SpanAfterInsert, int InsFam, int FamPos, int InsJob, int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC);

	int FindBestPosToInsertJobInCurFam_MS(const vector<int>& FamSeqInFac, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int Fam, int InsJob, int& BestPos);

	float FindBestPosToInsertJobInCurFam_TEC(const vector<int>& FamSeqInFac, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int Fam, int InsJob, int& BestPos, vector<float>& FacEC);

	//��һ��
	void Normalize(vector<Individual_PM>& Pop, int& nadirpointMS, float& nadirpointTEC, int& idealpointMS, float& idealpointTEC);
};
#endif

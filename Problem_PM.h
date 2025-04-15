#ifndef PROBLEM_PM_H
#define PROBLEM_PM_H
#include <string>
#include <vector>
using namespace std;
class Problem_PM
{
public:
	Problem_PM();
	virtual ~Problem_PM();

	void GenerateInstance();
	void ReadInstanceFileNameList(string Dir);
	void ReadInstance(int InsNo);
	void ReadSingleInstance();

protected:
	int m_Factories;
	int m_Machines;
	int m_Families;
	int m_Jobs;
	int m_SetupType;
	vector<string> m_InstanceFileNameList;
	vector<string> m_InstanceNameList;
	vector<vector<int>> m_JobsInEachFamily; //jobs in each family
	vector<vector<int>> m_JobOperPTime; //operation time on each machine for each job
	vector<int> m_MaxOfML;				//Maximum of maintenance level of machines
	vector<int> m_MT;						//The maintenance time on machines
	vector<vector<vector<int>>> m_SetupTime; //setup time between families
	vector<vector<vector<int>>> m_TureJobOpertime;  //变速之后工件的处理时间
	vector<float> m_Speed;  //速度
	float m_UnitIdleEC;//空闲能耗系数
	float m_UnitSetupEC;//准备能耗系数
	float m_UnitPMEC;//WE能耗系数
	vector<float> m_UnitPEC;
	vector<float> m_TEC;

	vector<vector<int>> m_JobSeqInFam; //Job sequence in each Family
	vector<int> m_JobTotalPTime; //all the operation time on all machines for a job;
	int m_AllJobTotalPTime;
	vector<int> m_FamTotalPTime; //all the operation time on all machines for the jobs in a family
	vector<double> m_FamAvgSetupTime; //average all the setup times for families



	void GetJobTotalPTime();
	void GetFamTotalPTime();
	void GetFamAvgSetupTime();
};
#endif

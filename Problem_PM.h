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
	vector<vector<vector<int>>> m_TureJobOpertime;  //����֮�󹤼��Ĵ���ʱ��
	vector<float> m_Speed;  //�ٶ�
	float m_UnitIdleEC;//�����ܺ�ϵ��
	float m_UnitSetupEC;//׼���ܺ�ϵ��
	float m_UnitPMEC;//WE�ܺ�ϵ��
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

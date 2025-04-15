#include <iostream>
#include <numeric>
#include <algorithm>

#include "NOperator_PM.h"
#include "Base.h"

using namespace std;
using namespace Base;

NOperator_PM::NOperator_PM()
{
}
NOperator_PM::~NOperator_PM()
{
}
vector<int> NOperator_PM::findMaxIndices(const vector<double>& vec)
{
	vector<int> maxIndices;
	if (vec.empty())
		return maxIndices; // 返回空的索引向量

	double maxVal = vec[0];
	maxIndices.push_back(0);

	for (int i = 1; i < vec.size(); ++i) {
		if (vec[i] > maxVal) {
			maxVal = vec[i];
			maxIndices.clear();
			maxIndices.push_back(i);
		}
		else if (vec[i] == maxVal) {
			maxIndices.push_back(i);
		}
	}

	return maxIndices;
}

/***********************************************************************************************************************/
/******************************************  对工件组里面的工件进行排序  ***********************************************/
void NOperator_PM::SortJobsInFam(int SortMethod, vector<vector<int>>& JobSeqInFam) //0:LPT,   1:SPT,  2:JobWeightTotalPTime, 
{
	JobSeqInFam = this->m_JobsInEachFamily; //initialize memeory
	//0:LPT,   1:SPT,  2:JobWeightTotalPTime, 3:first, 4:last, 5:first-last, 6:黄颖颖师姐
	if (SortMethod == 0) //LPT
	{
		for (auto& JobSeq : JobSeqInFam)
		{
			sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
				{
					return this->m_JobTotalPTime[j1] > this->m_JobTotalPTime[j2];
				});
		}
	}
	if (SortMethod == 1) //SPT
	{
		for (auto& JobSeq : JobSeqInFam)
		{
			sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
				{
					return this->m_JobTotalPTime[j1] < this->m_JobTotalPTime[j2];
				});
		}
	}

	
}

/***********************************************************************************************************************/
/***********************************************  对初始工件组进行排序  ************************************************/
void NOperator_PM::SortFam(int SortMethod, vector<int>& FamPermu)
{
	FamPermu.clear();
	FamPermu.resize(this->m_Families);
	for (int Fam = 0; Fam < this->m_Families; Fam++)
		FamPermu[Fam] = Fam;
	Pair<double>* ch = new Pair<double>[FamPermu.size()];

	if (SortMethod == 0) //LPT
	{
		for (int Fam = 0; Fam < FamPermu.size(); Fam++)
		{
			ch[Fam].dim = FamPermu[Fam];
			ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
		}
		sort(ch, ch + this->m_Families, PairGreater<double>());
	}
	if (SortMethod == 1) //SPT
	{
		for (int Fam = 0; Fam < FamPermu.size(); Fam++)
		{
			ch[Fam].dim = FamPermu[Fam];
			ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
		}
		sort(ch, ch + this->m_Families, PairLess<double>());
	}
	if (SortMethod == 2) //组平均准备时间，in increasing order
	{
		for (int Fam = 0; Fam < FamPermu.size(); Fam++)
		{
			ch[Fam].dim = FamPermu[Fam];
			ch[Fam].value = this->m_FamAvgSetupTime[FamPermu[Fam]];
		}
		sort(ch, ch + this->m_Families, PairLess<double>());
	}
	for (int Fam = 0; Fam < this->m_Families; Fam++)
		FamPermu[Fam] = ch[Fam].dim;
	delete[]ch;
}

/****************************************************************无快评**********************************************************************************************/
void NOperator_PM::Speed_mutation(vector<Individual_PM>& CMOEAPopulation, vector<Individual_PM>& Temp, vector<Individual_PM>& tureCMOEAPopulation)
{
	//对于NR_s中的个体进行速度突变，(γ,η,v)to(γ,η,v′) ，若新解占优原解，则用v'更新v 

	for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
	{
		//变速 及单位加工时间能耗
		for (int j = 0; j < m_Jobs; j++)
		{
			for (int m = 0; m < m_Machines; m++)
			{
				Temp[PS].m_SpeedVector[j][m] = rand() % 3;
			}
		}

		vector<int> FacSpan(m_Factories, 0);
		vector<float> FacEC(m_Factories, 0);
		int Makespan = 0;
		float TotalEC = 0;
		GetMSandTECForPerandToalFac(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, Temp[PS].m_SpeedVector, FacSpan, FacEC, Makespan, TotalEC);

		Temp[PS].MS = Makespan;
		Temp[PS].TEC = TotalEC;
		if (((Temp[PS].MS < CMOEAPopulation[PS].MS) && (Temp[PS].TEC < CMOEAPopulation[PS].TEC)) || ((Temp[PS].MS < CMOEAPopulation[PS].MS) && (Temp[PS].TEC == CMOEAPopulation[PS].TEC)) || ((Temp[PS].MS == CMOEAPopulation[PS].MS) && (Temp[PS].TEC < CMOEAPopulation[PS].TEC)))
		{
			CMOEAPopulation[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
			CMOEAPopulation[PS].MS = Temp[PS].MS;
			CMOEAPopulation[PS].TEC = Temp[PS].TEC;
			tureCMOEAPopulation.push_back(Temp[PS]);  //加入参考集
		}
	}

}

/**********************************************************************************************************************/
/********************************  计算最大完工时间FacCTime  **************************************/
int NOperator_PM::GetCTimeForAllFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix)
{
	int Makespan = 0;

	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		int TempMakespan = GetCTimeForPerFac(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix);
		if (Makespan < TempMakespan)
		{
			Makespan = TempMakespan;
		}
	}
	return Makespan;
}

int NOperator_PM::GetCTimeForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix)
{
	if (FamSeqInFac.empty() || JobSeqInFam.empty() || JobSeqInFam[0].empty()) {
		cerr << "Error: Invalid input sequence(s)." << endl;
		return -1; // 返回一个错误代码，或采取适当的错误处理措施
	}
	vector<int>MC(m_Machines, 0);
	vector<int> TempML = m_MaxOfML;

	vector<int> MachReadyTime(this->m_Machines, 0);
	for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
	{
		int CurFam = FamSeqInFac[Fam];

		if (Fam == 0) //the first group
		{	//第一组的准备时间
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
			}
		}
		else //from the second group to the end
		{
			int PreFam = FamSeqInFac[Fam - 1];
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
			}
		}
		for (int j = 0; j < JobSeqInFam[CurFam].size(); j++) //scheduling jobs in CurFam
		{
			int CurJob = JobSeqInFam[CurFam][j];

			for (int m = 0; m < this->m_Machines; m++) //on the rest machines
			{
				int StandOperTime = m_JobOperPTime[CurJob][m];//标准加工时间 没有除以速度之前的加工时间
				int CurOperTime = m_TureJobOpertime[CurJob][m][SpeedMatrix[CurJob][m]];//实际加工时间 标准加工时间除以速度之后的时间

				if (TempML[m] >= StandOperTime)
				{
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //on the first machine
					TempML[m] -= StandOperTime;
				}
				else
				{
					MachReadyTime[m] += m_MT[m];
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //加入维修时间
					TempML[m] = m_MaxOfML[m];
					TempML[m] -= StandOperTime;
					MC[m] += m_MT[m];
				}
			}
		}
	}

	int FacCTime = MachReadyTime[this->m_Machines - 1];
	return FacCTime;
}

/***********************************************************************************************************************/
/********************************  计算能耗TEC  **************************************/
float NOperator_PM::GetTECForAllFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix)
{
	float TEC = 0;
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		float TempTEC = GetTECForPerFac(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix);
		TEC += TempTEC;
		//cout << "工厂" << Fac << "的能耗：" << TempTEC << endl;
	}
	//cout << "总能耗：" << TEC << endl;
	return TEC;
}

float NOperator_PM::GetTECForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix)
{
	vector<int> MC(m_Machines, 0);
	vector<int> TempML = m_MaxOfML;
	vector<int> MachReadyTime(this->m_Machines, 0);
	for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
	{
		int CurFam = FamSeqInFac[Fam];

		if (Fam == 0) //the first group
		{
			//第一组的准备时间
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
				//cout << m << "机器：" << MachReadyTime[m] << endl;
			}
		}
		else //from the second group to the end
		{
			int PreFam = FamSeqInFac[Fam - 1];
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
			}
		}
		for (int j = 0; j < JobSeqInFam[CurFam].size(); j++) //scheduling jobs in CurFam
		{
			int CurJob = JobSeqInFam[CurFam][j];

			for (int m = 0; m < this->m_Machines; m++) //on the rest machines
			{
				int StandOperTime = m_JobOperPTime[CurJob][m];//标准加工时间 没有除以速度之前的加工时间
				int CurOperTime = m_TureJobOpertime[CurJob][m][SpeedMatrix[CurJob][m]];//实际加工时间 标准加工时间除以速度之后的时间
				//cout << j << "工件在" << m << "机器上的CurOperTime：" << CurOperTime << endl;
				if (TempML[m] >= StandOperTime)
				{
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //on the first machine
					TempML[m] -= StandOperTime;
					//cout << m << "机器TempML：" << TempML[m] << endl;
				}
				else
				{
					//cout << "维修了" << MT[m] << endl;
					MachReadyTime[m] += m_MT[m];
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //加入维修时间
					TempML[m] = m_MaxOfML[m];
					TempML[m] -= StandOperTime;
					MC[m] += m_MT[m];
					//cout << m << "机器TempML：" << TempML[m] << endl;
				}
			}
		}
	}

	//计算能耗
	float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

	for (int m = 0; m < m_Machines; m++)
	{
		int preFam = -1;
		int Tptime = 0, Tsetuptime = 0, Tidletime = 0, Tpmtime = 0;
		int Fam = FamSeqInFac[0];
		for (int g = 0; g < FamSeqInFac.size(); g++)
		{
			Fam = FamSeqInFac[g];

			if (preFam == -1)
			{
				TSetupTimeEC += m_SetupTime[m][Fam][Fam] * m_UnitSetupEC;
				Tsetuptime += m_SetupTime[m][Fam][Fam];
			}
			else
			{
				TSetupTimeEC += m_SetupTime[m][preFam][Fam] * m_UnitSetupEC;
				Tsetuptime += m_SetupTime[m][preFam][Fam];
			}
			preFam = Fam;
			for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
			{
				int Job = JobSeqInFam[Fam][j];
				TpTimeEC += m_TureJobOpertime[Job][m][SpeedMatrix[Job][m]] * m_UnitPEC[SpeedMatrix[Job][m]];
				Tptime += m_TureJobOpertime[Job][m][SpeedMatrix[Job][m]];
			}

		}
		Tpmtime += MC[m];
		if (m > 0)
		{
			Tidletime += (MachReadyTime[m] - Tptime - Tsetuptime - Tpmtime);
			TIdleTimeEC += Tidletime * m_UnitIdleEC;
		}

	}
	float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;

	return FacEC;
}

/***********************************************************************************************************************/
/********************************  计算空闲时间（最后速度策略时计算使用）  **************************************/
int NOperator_PM::GetDelayTimeAllFac(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<int>& FacSpan, vector<vector<int>>& DelayTime)//Forward pass calculation
{
	int min = -1;
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		FacSpan[Fac] = GetDelayTimePerFac(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix, DelayTime);//得到每个工厂的完工时间
		if (FacSpan[Fac] > min)
		{
			min = FacSpan[Fac];
		}
	}
	return min; //得到关键工厂
}

int NOperator_PM::GetDelayTimePerFac(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, vector<vector<int>>& DelayTime)
{
	vector<vector<int>> JCTime(this->m_Jobs);
	for (int j = 0; j < JCTime.size(); j++)
		JCTime[j].resize(this->m_Machines, 0);
	vector<int> MachReadyTime(this->m_Machines, 0);   //机器准备时间初始化为0
	vector<int> TempML = m_MaxOfML;   // 用来记录机器的剩余操作时间
	vector<int> MC(this->m_Machines, 0); // 用于统计维修次数

	for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
	{
		int CurFam = FamSeqInFac[Fam];  //当前组
		if (Fam == 0) //the first group of jobs 第一个组
		{
			for (int m = 0; m < this->m_Machines; m++)
				MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
		}
		else //第二到最后组 from the second group of jobs to the end;
		{
			int PreFam = FamSeqInFac[Fam - 1];
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
			}

		}

		for (int j = 0; j < JobSeqInFam[CurFam].size(); j++) //当前组的工件调度 Scheduling Jobs in CurFam
		{
			int CurJob = JobSeqInFam[CurFam][j];  //当前工件
			int StandOperTime = m_JobOperPTime[CurJob][0];//标准加工时间 没有除以速度之前的加工时间
			int CurOperTime = m_TureJobOpertime[CurJob][0][SpeedMatrix[CurJob][0]];//实际加工时间 标准加工时间除以速度之后的时间
			if (TempML[0] >= StandOperTime)
			{
				JCTime[CurJob][0] = MachReadyTime[0] + CurOperTime; //在第一个机器上的完工时间 on the first machine
				TempML[0] -= StandOperTime;
			}
			else
			{
				bool speedAdjusted = false;
				for (int s = SpeedMatrix[CurJob][0] + 1; s < m_Speed.size(); s++) // Try increasing speed
				{
					int adjustedCurOperTime = m_TureJobOpertime[CurJob][0][s];
					if (TempML[0] >= StandOperTime)
					{
						// Speed adjustment successful
						m_TureJobOpertime[CurJob][0][SpeedMatrix[CurJob][0]] = adjustedCurOperTime;
						JCTime[CurJob][0] = MachReadyTime[0] + m_TureJobOpertime[CurJob][0][SpeedMatrix[CurJob][0]];
						TempML[0] -= StandOperTime;
						speedAdjusted = true;
						break;
					}
				}
				if (!speedAdjusted)
				{
					// No feasible speed adjustment, proceed with maintenance
					JCTime[CurJob][0] = MachReadyTime[0] + m_MT[0] + m_TureJobOpertime[CurJob][0][SpeedMatrix[CurJob][0]];
					TempML[0] = m_MaxOfML[0];
					TempML[0] -= StandOperTime;
				}
			}
			for (int m = 1; m < this->m_Machines; m++) //on the rest machines
			{
				StandOperTime = m_JobOperPTime[CurJob][m];//标准加工时间 没有除以速度之前的加工时间
				CurOperTime = m_TureJobOpertime[CurJob][m][SpeedMatrix[CurJob][m]];//实际加工时间 标准加工时间除以速度之后的时间
				if (TempML[m] >= StandOperTime)
				{
					JCTime[CurJob][m] = (m == 0 ? MachReadyTime[m] : max(JCTime[CurJob][m - 1], MachReadyTime[m])) + CurOperTime; //on the first machine
					if ((j > 0) && (JCTime[CurJob][m - 1] > MachReadyTime[m]))//准备时间没办法改变速度,所以j>0
					{
						int PreJob = JobSeqInFam[CurFam][j - 1];
						DelayTime[PreJob][m] = JCTime[CurJob][m - 1] - MachReadyTime[m];
					}
					TempML[m] -= StandOperTime;
				}
				else
				{
					bool speedAdjusted = false;
					for (int s = SpeedMatrix[CurJob][m] + 1; s < m_Speed.size(); s++) // Try increasing speed
					{
						int adjustedCurOperTime = m_TureJobOpertime[CurJob][m][s];
						if (TempML[m] >= StandOperTime)
						{
							// Speed adjustment successful
							m_TureJobOpertime[CurJob][m][SpeedMatrix[CurJob][m]] = adjustedCurOperTime;
							JCTime[CurJob][m] = (m == 0 ? MachReadyTime[m] : max(JCTime[CurJob][m - 1], MachReadyTime[m])) + adjustedCurOperTime;
							TempML[m] -= StandOperTime;
							speedAdjusted = true;
							break;
						}
					}
					if (!speedAdjusted)
					{
						// No feasible speed adjustment, proceed with maintenance
						JCTime[CurJob][m] = (m == 0 ? MachReadyTime[m] : max(JCTime[CurJob][m - 1], MachReadyTime[m])) + m_MT[m] + CurOperTime;
						if ((j > 0) && (JCTime[CurJob][m - 1] > MachReadyTime[m]))
						{
							int PreJob = JobSeqInFam[CurFam][j - 1];
							DelayTime[PreJob][m] = JCTime[CurJob][m - 1] - MachReadyTime[m] + m_MT[m]; // 加入维修时间的延迟
						}
						TempML[m] = m_MaxOfML[m];
						TempML[m] -= StandOperTime;
					}
				}
			}
			MachReadyTime = JCTime[CurJob];
		}
	}

	// 得到工件完工时间
	return MachReadyTime[this->m_Machines - 1];
}

/***********************************************************************************************************************/
/********************************  计算关键工厂和总能耗  **************************************/
void NOperator_PM::GetMSandTECForPerandToalFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC)
{
	//得到关键工厂和总的TEC
	int max = -1;
	for (int fac = 0; fac < FacFamSeq.size(); fac++)
	{
		auto re = GetCTimeAndTECForPerFac(FacFamSeq[fac], JobSeqInFam, SpeedMatrix);

		FacSpan[fac] = re.first;
		FacEC[fac] = re.second;

		if (FacSpan[fac] > max)
		{
			max = FacSpan[fac];
		}
	}
	Makespan = max; //得到关键工厂
	TotalEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);
}

pair<int, float> NOperator_PM::GetCTimeAndTECForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix)
{
	vector<int> MC = vector<int>(m_Machines, 0);
	vector<int> TempML = m_MaxOfML;
	vector<int> MachReadyTime(this->m_Machines, 0);
	for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
	{
		int CurFam = FamSeqInFac[Fam];
		if (CurFam < 0 || CurFam >= JobSeqInFam.size()) {
			std::cerr << "Error: CurFam index out of bounds. CurFam = " << CurFam << std::endl;
			continue; // 跳过这个组，防止越界错误
		}

		if (Fam == 0) //the first group
		{	//第一组的准备时间
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
			}
		}
		else //from the second group to the end
		{
			int PreFam = FamSeqInFac[Fam - 1];
			for (int m = 0; m < this->m_Machines; m++)
			{
				MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
			}
		}
		for (int j = 0; j < JobSeqInFam[CurFam].size(); j++) //scheduling jobs in CurFam
		{
			int CurJob = JobSeqInFam[CurFam][j];

			for (int m = 0; m < this->m_Machines; m++) //on the rest machines
			{
				int StandOperTime = m_JobOperPTime[CurJob][m];//标准加工时间 没有除以速度之前的加工时间
				int CurOperTime = m_TureJobOpertime[CurJob][m][SpeedMatrix[CurJob][m]];//实际加工时间 标准加工时间除以速度之后的时间

				if (TempML[m] >= StandOperTime)
				{
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //on the first machine
					TempML[m] -= StandOperTime;
				}
				else
				{
					MachReadyTime[m] += m_MT[m];
					MachReadyTime[m] = (m == 0 ? MachReadyTime[m] : max(MachReadyTime[m - 1], MachReadyTime[m])) + CurOperTime; //加入维修时间
					TempML[m] = m_MaxOfML[m];
					TempML[m] -= StandOperTime;
					MC[m] += m_MT[m];
				}
			}
		}
	}

	//计算能耗
	float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

	for (int m = 0; m < m_Machines; m++)
	{
		int preFam = -1;
		int Tptime = 0, Tsetuptime = 0, Tidletime = 0, Tpmtime = 0;
		int Fam = FamSeqInFac[0];
		for (int g = 0; g < FamSeqInFac.size(); g++)
		{
			Fam = FamSeqInFac[g];

			if (preFam == -1)
			{
				TSetupTimeEC += m_SetupTime[m][Fam][Fam] * m_UnitSetupEC;
				Tsetuptime += m_SetupTime[m][Fam][Fam];
			}
			else
			{
				TSetupTimeEC += m_SetupTime[m][preFam][Fam] * m_UnitSetupEC;
				Tsetuptime += m_SetupTime[m][preFam][Fam];
			}
			preFam = Fam;
			for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
			{
				int Job = JobSeqInFam[Fam][j];
				TpTimeEC += m_TureJobOpertime[Job][m][SpeedMatrix[Job][m]] * m_UnitPEC[SpeedMatrix[Job][m]];
				Tptime += m_TureJobOpertime[Job][m][SpeedMatrix[Job][m]];
			}

		}
		Tpmtime += MC[m];
		if (m > 0)
		{
			Tidletime += (MachReadyTime[m] - Tptime - Tsetuptime - Tpmtime);
			TIdleTimeEC += Tidletime * m_UnitIdleEC;
		}

	}
	float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;

	int FacCTime = MachReadyTime[this->m_Machines - 1];
	return { FacCTime, FacEC };
}

/**********************************************************************************************************************************/
/*********************************  对待调度工件组InsFam在工厂FacFamSeq中寻找最好位置(Makespan)  ******************************/
int NOperator_PM::FindBestPosToInsertFamForAllFac_Makespan(const vector <vector <int>>& FacFamSeq, const vector <vector <int>>& JobSeqInFam,
	const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos)
{
	int minMakespan = INT_MAX;
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		int Pos = -1;
		int Makespan = this->FindBestPosToInsertFamForPerFac_Makespan(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix, InsFam, Pos);
		if (Makespan < minMakespan)
		{
			minMakespan = Makespan;
			BestFac = Fac;
			BestPos = Pos;
		}
	}
	return minMakespan;
}

int NOperator_PM::FindBestPosToInsertFamForPerFac_Makespan(const vector <int>& NewFamSeqInFac, const vector <vector <int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos)
{
	int minFacMakespan = INT_MAX;
	vector <int> TempFamSeqInFac = NewFamSeqInFac;
	for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
	{
		TempFamSeqInFac.insert(TempFamSeqInFac.begin() + Pos, InsFam);
		int FacMakespan = this->GetCTimeForPerFac(TempFamSeqInFac, JobSeqInFam, SpeedMatrix);
		if (FacMakespan < minFacMakespan)
		{
			minFacMakespan = FacMakespan;
			BestPos = Pos;
		}
		TempFamSeqInFac.erase(TempFamSeqInFac.begin() + Pos);
	}
	return minFacMakespan;
}

/*************************************************************************************************************************/
/*******************************  对待调度工件组InsFam在工厂FamSeqInFac中寻找最好位置(TEC)  ******************/
float NOperator_PM::FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos)
{
	float minTEC = INT_MAX;
	//计算未插入之前各工厂的EC
	vector <float> FacEC(FacFamSeq.size(), 0.0);
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		FacEC[Fac] = GetTECForPerFac(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix);
		//cout << "未插入时工厂" << Fac << "的EC：" << FacEC[Fac] << endl;
	}

	vector <float> TempFacEC(FacFamSeq.size(), 0.0);
	TempFacEC = FacEC;

	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		int Pos = -1;
		TempFacEC[Fac] = this->FindBestPosToInsertFamForPerFac_TEC(FacFamSeq[Fac], JobSeqInFam, SpeedMatrix, InsFam, Pos);

		float TEC = accumulate(TempFacEC.begin(), TempFacEC.end(), 0);

		if (TEC < minTEC)
		{
			minTEC = TEC;
			BestFac = Fac;
			BestPos = Pos;
		}
		TempFacEC[Fac] = FacEC[Fac];
	}
	return minTEC;
}

float NOperator_PM::FindBestPosToInsertFamForPerFac_TEC(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos)
{
	float minFacEC = INT_MAX;
	vector <int> TempFamSeqInFac = NewFamSeqInFac;
	for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
	{
		TempFamSeqInFac.insert(TempFamSeqInFac.begin() + Pos, InsFam);
		float FacEC = GetTECForPerFac(TempFamSeqInFac, JobSeqInFam, SpeedMatrix);
		if (FacEC < minFacEC)
		{
			minFacEC = FacEC;
			BestPos = Pos;
		}
		TempFamSeqInFac.erase(TempFamSeqInFac.begin() + Pos);
	}
	return minFacEC;
}

/***********************************************************************************************************************/
/********************************  基于收敛指标寻找最好位置插入（组）  **************************************/
float NOperator_PM::FindBestPosToInsertFamForAllFac_Ind(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacTEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation)
{
	float minInd = FLT_MAX;
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		int Pos = -1;
		float Ind = this->FindBestPosToInsertFamForPerFac_Ind(Fac, FacFamSeq, JobSeqInFam, SpeedMatrix, FacSpan, FacTEC, InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);
		if (Ind < minInd)
		{
			minInd = Ind;
			BestFac = Fac;
			BestPos = Pos;
		}
	}
	return minInd;
}

float NOperator_PM::FindBestPosToInsertFamForPerFac_Ind(int Fac, const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	vector<int>& FacSpan, vector<float>& FacTEC, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation)
{
	float minFacInd = FLT_MAX;
	vector<vector<int>> TempFacFamSeq = FacFamSeq;
	for (int Pos = 0; Pos <= FacFamSeq[Fac].size(); Pos++)
	{
		TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);
		float FacInd = GetindForPerFacAfterInsertFam(Fac, TempFacFamSeq, JobSeqInFam, SpeedMatrix, InsFam, Pos,
			nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacTEC, OrgMS, OrgTEC, CMOEAPopulation);
		TempFacFamSeq[Fac].erase(TempFacFamSeq[Fac].begin() + Pos);
		if (FacInd < minFacInd)
		{
			minFacInd = FacInd;
			BestPos = Pos;
		}
	}
	return minFacInd;
}

float NOperator_PM::GetindForPerFacAfterInsertFam(int Fac, const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int Pos,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacEC, int OrgMS, float OrgTEC, vector<Individual_PM>& CMOEAPopulation)
{

	int MSAfterInsert = 0;
	float TotalTECAfterInsert = 0.0;
	GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, SpeedMatrix, FacSpan, FacEC, MSAfterInsert, TotalTECAfterInsert);

	//cout << "总能耗：" << TEC << endl;
	//判断是否替换最低点和理想点
	if (MSAfterInsert > nadirpointMS)
		nadirpointMS = MSAfterInsert;
	if (MSAfterInsert < idealpointMS)
		idealpointMS = MSAfterInsert;
	if (TotalTECAfterInsert > nadirpointTEC)
		nadirpointTEC = TotalTECAfterInsert;
	if (TotalTECAfterInsert < idealpointTEC)
		idealpointTEC = TotalTECAfterInsert;

	//判断是否比操作前改进，若改进则加入参考集
	bool flag = true;
	if ((MSAfterInsert < OrgMS && TotalTECAfterInsert <= OrgTEC) || (MSAfterInsert <= OrgMS && TotalTECAfterInsert < OrgTEC))
	{
		for (int i = 0; i < CMOEAPopulation.size(); i++)
		{
			if ((MSAfterInsert == CMOEAPopulation[i].MS) && (TotalTECAfterInsert == CMOEAPopulation[i].TEC))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			//cout << "*********************" << endl;
			Individual_PM tempIndi;
			tempIndi.m_FacFamSeqArray = FacFamSeq;
			tempIndi.m_JobSeqInFamArray = JobSeqInFam;
			tempIndi.MS = MSAfterInsert;
			tempIndi.TEC = TotalTECAfterInsert;
			tempIndi.m_SpeedVector = SpeedMatrix;  //速度矩阵
			CMOEAPopulation.push_back(tempIndi);
		}
	}

	//归一化
	float normalMS = -1;
	float normalTEC = -1;
	if (nadirpointMS != idealpointMS)
		normalMS = (static_cast<float>(MSAfterInsert - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	else
		normalMS = 0.0f; // 或者设置为另一个合适的值

	if (nadirpointTEC != idealpointTEC)
		normalTEC = (TotalTECAfterInsert - idealpointTEC) / (nadirpointTEC - idealpointTEC);
	else
		normalTEC = 0.0f; // 或者设置为另一个合适的值

	//cout << "个体归一化后的MS：" << normalMS << "\tTEC：" << normalTEC << endl;
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	int prank = Individual_PM::GetPRank(CMOEAPopulation, normalMS, normalTEC, MSAfterInsert, TotalTECAfterInsert);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, normalMS, normalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	float Orgconvergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);
	//cout << "Orgconvergence_ind:" << Orgconvergence_ind << endl;
	return Orgconvergence_ind;
}

float NOperator_PM::FindBestPosToInsertFamForAllFac_Ind_DR(vector<Individual_PM>& CMOEAPopulation, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
	float minInd = INT_MAX;
	for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
	{
		int Pos = -1;
		float Ind = this->FindBestPosToInsertFamForPerFac_Ind_DR(CMOEAPopulation, Fac, FacFamSeq, JobSeqInFam, SpeedMatrix, InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
		if (Ind < minInd)
		{
			minInd = Ind;
			BestFac = Fac;
			BestPos = Pos;
		}
	}
	return minInd;
}

float NOperator_PM::FindBestPosToInsertFamForPerFac_Ind_DR(vector<Individual_PM>& CMOEAPopulation, int Fac, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
	float minFacInd = INT_MAX;
	for (int Pos = 0; Pos <= FacFamSeq[Fac].size(); Pos++)
	{
		float FacInd = GetindForPerFacAfterInsertFam_DR(CMOEAPopulation, Fac, FacFamSeq, JobSeqInFam, SpeedMatrix, InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
		FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);
		if (FacInd < minFacInd)
		{
			minFacInd = FacInd;
			BestPos = Pos;
		}
	}
	return minFacInd;
}

float NOperator_PM::GetindForPerFacAfterInsertFam_DR(vector<Individual_PM>& CMOEAPopulation, int Fac, vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
	const vector<vector<int>>& SpeedMatrix, int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
	FacFamSeq[Fac].insert(FacFamSeq[Fac].begin() + Pos, InsFam);//InsFam就是经过LPT排序之后的那个组序列中按顺序的某个组
	vector <int> FacSpan(FacFamSeq.size(), 0);
	vector<float> FacEC(FacFamSeq.size(), 0);
	int Makespan;
	float TEC;
	GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, SpeedMatrix, FacSpan, FacEC, Makespan, TEC);

	//判断是否替换最低点和理想点
	if (Makespan > nadirpointMS)
		nadirpointMS = Makespan;
	if (Makespan < idealpointMS)
		idealpointMS = Makespan;

	//判断是否替换最低点和理想点
	if (TEC > nadirpointTEC)
		nadirpointTEC = TEC;
	if (TEC < idealpointTEC)
		idealpointTEC = TEC;
	//归一化
	float normalMS = -1;
	float normalTEC = -1;
	normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	//计算指标Ind
	float convergence_ind = 0;
	int prank = Individual_PM::GetPRank(CMOEAPopulation, normalMS, normalTEC, Makespan, TEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, normalMS, normalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	convergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);
	return convergence_ind;
}

/************************************************************************************************************************/
/****************************  对工件InsJob在工件组JobSeqInFam[Fam]中寻找最佳位置  **************************************/
void NOperator_PM::Basedind_JobsInFam(int fac, int Fam, const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual_PM>& CMOEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int& ObjectMS, float& ObjectTEC)
{
	vector<int> JobSeqsInFamForExtracted = JobSeqInFam[Fam];
	random_shuffle(JobSeqsInFamForExtracted.begin(), JobSeqsInFamForExtracted.end());
	int FamPos = find(begin(FacFamSeq[fac]), end(FacFamSeq[fac]), Fam) - begin(FacFamSeq[fac]);
	int nCnt = 0, CurPos = 0, BestPos = -1;
	int Makespan = -1;

	while (nCnt < JobSeqInFam[Fam].size())
	{
		CurPos = CurPos % JobSeqsInFamForExtracted.size();
		int CurJob = JobSeqsInFamForExtracted[CurPos]; 
		vector<int>::iterator it = find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), CurJob);
		JobSeqInFam[Fam].erase(it);

		float minFacInd = INT_MAX;
		for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
		{
			float FacInd = GetindForPerFacAfterInsertJob(fac, FacFamSeq, JobSeqInFam, SpeedMatrix, Fam, CurJob, Pos,
				nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacEC, ObjectMS, ObjectTEC, CMOEAPopulation);
			JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + Pos);
			if (FacInd < minFacInd)
			{
				minFacInd = FacInd;
				BestPos = Pos;
			}

		}
		JobSeqInFam[Fam].insert(JobSeqInFam[Fam].begin() + BestPos, CurJob);


		vector<int> AfterInsertSpanFac(FacFamSeq.size(), 0);
		vector<float> AfterInsertECFac(FacFamSeq.size(), 0);
		int AfterInsertMS = -1;
		float AfterInsertTEC = -1;
		GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, SpeedMatrix, AfterInsertSpanFac, AfterInsertECFac, AfterInsertMS, AfterInsertTEC);

		ObjectMS = AfterInsertMS;
		ObjectTEC = AfterInsertTEC;

		CurPos++;
		nCnt++;

	}

}

float NOperator_PM::GetindForPerFacAfterInsertJob(int fac, const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& SpeedMatrix,
	int InsFam, int InsJob, int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<int>& FacSpan, vector<float>& FacEC, int& OrgMS, float& OrgTEC, vector<Individual_PM>& CMOEAPopulation)
{
	JobSeqInFam[InsFam].insert(JobSeqInFam[InsFam].begin() + JobPos, InsJob);
	pair<int, float> MsandTEC = GetCTimeAndTECForPerFac(FacFamSeq[fac], JobSeqInFam, SpeedMatrix);
	FacSpan[fac] = MsandTEC.first;
	FacEC[fac] = MsandTEC.second;
	float TEC = accumulate(FacEC.begin(), FacEC.end(), 0);
	int Makespan = *max_element(FacSpan.begin(), FacSpan.end());

	//判断是否替换最低点和理想点
	if (Makespan > nadirpointMS)
		nadirpointMS = Makespan;
	if (Makespan < idealpointMS)
		idealpointMS = Makespan;

	//判断是否替换最低点和理想点
	if (TEC > nadirpointTEC)
		nadirpointTEC = TEC;
	if (TEC < idealpointTEC)
		idealpointTEC = TEC;

	//判断是否比操作前改进，若改进则加入参考集
	bool flag = true;
	if ((Makespan < OrgMS && TEC <= OrgTEC) || (Makespan <= OrgMS && TEC < OrgTEC))
	{
		for (int i = 0; i < CMOEAPopulation.size(); i++)
		{
			if ((Makespan == CMOEAPopulation[i].MS) && (TEC == CMOEAPopulation[i].TEC))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			//cout << "*********************" << endl;
			Individual_PM tempIndi;
			tempIndi.m_FacFamSeqArray = FacFamSeq;
			tempIndi.m_JobSeqInFamArray = JobSeqInFam;
			tempIndi.m_SpeedVector = SpeedMatrix;  //速度矩阵
			tempIndi.MS = Makespan;
			tempIndi.TEC = TEC;
			CMOEAPopulation.push_back(tempIndi);
			OrgMS = Makespan;
			OrgTEC = TEC;
		}
	}

	//归一化
	float normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	float normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	//cout << "个体归一化后的MS：" << normalMS << "\tTEC：" << normalTEC << endl;
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	int prank = Individual_PM::GetPRank(CMOEAPopulation, normalMS, normalTEC, Makespan, TEC);
	pair<int, int> re = Individual_PM::GetLRank(CMOEAPopulation, normalMS, normalTEC);
	int lrank = re.first;
	int nbrmax = re.second;
	float convergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);
	return convergence_ind;
}

float NOperator_PM::GetindForPerFacAfterInsertJob_DR(vector<Individual_PM>& CMOEAPopulation, int fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
	const vector<vector<int>>& SpeedMatrix, int SpanAfterInsert, int InsFam, int FamPos, int InsJob, int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
	//计算未插入之前各工厂的span
	vector <int> FacSpan(FacFamSeq.size(), 0);
	for (int f = 0; f < FacFamSeq.size(); f++)
	{
		FacSpan[f] = GetCTimeForPerFac(FamSeqInFac, JobSeqInFam, SpeedMatrix);
		//cout << "未插入时工厂" << f << "的Span：" << FacSpan[f] << endl;
	}
	FacSpan[fac] = SpanAfterInsert;

	//cout << "插入工厂" << fac << "后的Span：" << FacSpan[fac] << endl;
	int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
	//cout << "插入的MakeSpan：" << Makespan << endl;

	//判断是否替换最低点和理想点
	if (Makespan > nadirpointMS)
		nadirpointMS = Makespan;
	if (Makespan < idealpointMS)
		idealpointMS = Makespan;

	vector<float> FacEC(FacFamSeq.size(), 0);
	float TEC;
	//计算能耗
	vector<vector<int>> TempFacFamSeq;
	TempFacFamSeq.clear();
	TempFacFamSeq.resize(FacFamSeq.size());
	TempFacFamSeq = FacFamSeq;

	vector<vector<int>> TempJobSeqInFam;
	TempJobSeqInFam.clear();
	TempJobSeqInFam.resize(JobSeqInFam.size());
	TempJobSeqInFam = JobSeqInFam;
	TempJobSeqInFam[InsFam].insert(TempJobSeqInFam[InsFam].begin() + JobPos, InsJob);

	for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
	{
		FacEC[fac] = GetTECForPerFac(TempFacFamSeq[fac], TempJobSeqInFam, SpeedMatrix);
	}

	TEC = accumulate(FacEC.begin(), FacEC.end(), 0);
	//cout << "总能耗：" << TEC << endl;

	//判断是否替换最低点和理想点
	if (TEC > nadirpointTEC)
		nadirpointTEC = TEC;
	if (TEC < idealpointTEC)
		idealpointTEC = TEC;


	//归一化
	float normalMS = -1;
	float normalTEC = -1;
	//cout << endl << "normalize" << endl;

	normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
	normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
	//cout << "个体归一化后的MS：" << normalMS << "\tTEC：" << normalTEC << endl;
	Normalize(CMOEAPopulation, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
	//计算指标Ind
	float convergence_ind = 0;
	int prank = Individual_PM::GetPRank(CMOEAPopulation, normalMS, normalTEC, Makespan, TEC);
	pair<int, int> re1 = Individual_PM::GetLRank(CMOEAPopulation, normalMS, normalTEC);
	int lrank = re1.first;
	int nbrmax = re1.second;
	convergence_ind = 2 * (prank - 1) + (1.0 / (nbrmax + 1)) * (lrank - 1);
	return convergence_ind;
}

int NOperator_PM::FindBestPosToInsertJobInCurFam_MS(const vector<int>& FamSeqInFac, vector<vector<int>> JobSeqInFam, const vector<vector<int>>& SpeedMatrix, int Fam, int InsJob, int& BestPos)
{
	int min = INT_MAX;
	for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
	{
		//插入
		JobSeqInFam[Fam].insert(JobSeqInFam[Fam].begin() + Pos, InsJob);
		int span = this->GetCTimeForPerFac(FamSeqInFac, JobSeqInFam, SpeedMatrix);
		if (span < min)
		{
			min = span;
			BestPos = Pos;
		}
		JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + Pos);
	}
	return min;
}

float NOperator_PM::FindBestPosToInsertJobInCurFam_TEC(const vector<int>& FamSeqInFac, vector<vector<int>> JobSeqInFam,
	const vector<vector<int>>& SpeedMatrix, int Fam, int InsJob, int& BestPos, vector<float>& FacEC)
{
	float min = FLT_MAX;
	for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
	{
		//插入
		JobSeqInFam[Fam].insert(JobSeqInFam[Fam].begin() + Pos, InsJob);
		float Factec = this->GetTECForPerFac(FamSeqInFac, JobSeqInFam, SpeedMatrix);
		float TotalTEC = accumulate(FacEC.begin(), FacEC.end(), 0);
		if (TotalTEC < min)
		{
			min = TotalTEC;
			BestPos = Pos;
		}
		JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + Pos);
	}
	return min;
}

//归一化
void NOperator_PM::Normalize(vector<Individual_PM>& Pop, int& nadirpointMS, float& nadirpointTEC, int& idealpointMS, float& idealpointTEC)
{

	int minMS = INT_MAX;
	float minTEC = FLT_MAX;

	int maxMS = INT_MIN;
	float maxTEC = FLT_MIN;

	for (int i = 0; i < Pop.size(); i++)
	{
		if (Pop[i].MS < minMS)
			minMS = Pop[i].MS;

		if (Pop[i].MS > maxMS)
			maxMS = Pop[i].MS;

		if (Pop[i].TEC < minTEC)
			minTEC = Pop[i].TEC;

		if (Pop[i].TEC > maxTEC)
			maxTEC = Pop[i].TEC;

		//标志设为0
		Pop[i].flag = 0;

	}
	nadirpointMS = maxMS;
	nadirpointTEC = maxTEC;
	idealpointMS = minMS;
	idealpointTEC = minTEC;


	//归一化：将值放在0-1范围内
	for (int i = 0; i < Pop.size(); i++)
	{
		if (maxMS == minMS)
			Pop[i].normalMS = 0.0f;  // 所有 MS 相等，归一化为0
		else
			Pop[i].normalMS = static_cast<float>(Pop[i].MS - minMS) / static_cast<float>(maxMS - minMS);

		if (maxTEC == minTEC)
			Pop[i].normalTEC = 0.0f;  // 所有 TEC 相等，归一化为0
		else
			Pop[i].normalTEC = (Pop[i].TEC - minTEC) / (maxTEC - minTEC);
		//cout << Pop[i].normalMS << "  " << Pop[i].normalTEC << endl;
	}
}
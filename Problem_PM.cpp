#include "Problem_PM.h"
#include <fstream>
#include <iostream>
#include <sstream>

Problem_PM::Problem_PM()
{
}
Problem_PM::~Problem_PM()
{
}
void Problem_PM::GenerateInstance()
{
	//vector<int> FactoryArray = { 2,3,4,5,6,7 }; //共产生6*3*3*3*5=810个算例
	vector<int> FactoryArray = { 2,3,4 }; //共产生3*3*3*3*10=810个算例
	vector<int> MachineArray = { 2,4,6 };
	vector<int> FamilyArray = { 20,40,60 };
	vector<int> SetupTypeArray = { 0,1,2 }; //Class 0: [1,20]; Class 1: [1,50]; Class 2: [1,100];
	int minSetupTime = 1;
	vector<int> maxSetupTimeArray = { 20,50,100 };
	int minJobsInFamily = 1;
	int maxJobsInFamily = 10;
	int minProcessingTime = 1;
	int maxProcessingTime = 10;
	int InsRep = 5;

	ofstream ofile1, ofile2;
	ofile1.open("..\\Benchmark\\InstanceFileNameList.txt");

	int Counter_Files = 1;
	for (int i = 0; i < FactoryArray.size(); i++)
	{
		for (int j = 0; j < MachineArray.size(); j++)
		{
			for (int k = 0; k < FamilyArray.size(); k++)
			{
				for (int s = 0; s < SetupTypeArray.size(); s++)
				{
					for (int rep = 0; rep < InsRep; rep++)
					{
						int Factories = FactoryArray[i];
						int Machines = MachineArray[j];
						int Families = FamilyArray[k];
						int SetupType = SetupTypeArray[s];
						int maxSetupTime = maxSetupTimeArray[SetupType];

						ostringstream FileName;
						FileName << "FG_" << Factories << "_" << Machines << "_" << Families << "_" << SetupType + 1 << "_" << rep + 1 << ".txt";
						cout << Counter_Files << "\t" << FileName.str() << endl;
						ofile1 << Counter_Files << "\t" << FileName.str() << endl;
						ofile2.open("..\\Benchmark\\" + FileName.str());
						ofile2 << "Factories" << "\t" << Factories << endl;
						ofile2 << "Machines" << "\t" << Machines << endl;
						ofile2 << "Families" << "\t" << Families << endl;
						ofile2 << "SetupType" << "\t" << SetupType + 1 << endl;

						//Number of jobs in each family
						ofile2 << endl << "Number of Jobs in each Family" << endl;
						vector<vector<int>> JobsInEachFamily(Families);
						int nJobs = 0;
						for (int Fam = 0; Fam < Families; Fam++)
						{
							int nJobsInFamily = rand() % (maxJobsInFamily - minJobsInFamily + 1) + minJobsInFamily;
							for (int xx = 0; xx < nJobsInFamily; xx++)
								JobsInEachFamily[Fam].push_back(nJobs + xx);
							nJobs += nJobsInFamily;
							ofile2 << nJobsInFamily << "\t";
						}
						ofile2 << endl;

						int maxjobsingroup = 0;
						int minjobsingroup = 0;
						//Jobs in each Family
						ofile2 << endl << "Jobs in each Family" << endl;
						for (int Fam = 0; Fam < Families; Fam++)
						{
							for (int Job = 0; Job < JobsInEachFamily[Fam].size(); Job++)
								ofile2 << JobsInEachFamily[Fam][Job] << "\t";
							ofile2 << endl;
						}

						for (int Fam = 0; Fam < Families; Fam++)
						{
							if (maxjobsingroup < JobsInEachFamily[Fam].size()) {
								maxjobsingroup = JobsInEachFamily[Fam].size();
							}
							if (minjobsingroup > JobsInEachFamily[Fam].size()) {
								minjobsingroup = JobsInEachFamily[Fam].size();
							}
						}
						//Total number of jobs
						ofile2 << endl << "Total number of jobs" << endl << nJobs << endl;

						vector<vector<int>> Processingtimes(nJobs, vector<int>(Machines, 0));
						int times = 0;
						//Processing time for jobs
						ofile2 << endl << "Processing times of jobs" << endl;
						for (int Job = 0; Job < nJobs; Job++)
						{
							for (int machine = 0; machine < Machines; machine++) {
								times = rand() % (maxProcessingTime - minProcessingTime + 1) + minProcessingTime;
								ofile2 << times << "\t";
								Processingtimes[Job][machine] = times;
							}
							ofile2 << endl;
						}

						vector<vector<int>> timeofmachine(Machines, vector<int>(2, 0));//每台机器上最大的加工时间和最小的加工时间 0大1小
						for (int machine = 0; machine < Machines; machine++) {
							for (int Job = 0; Job < nJobs; Job++)
							{
								if (Processingtimes[Job][machine] > timeofmachine[machine][0]) {
									timeofmachine[machine][0] = Processingtimes[Job][machine];
								}
								if (Processingtimes[Job][machine] < timeofmachine[machine][1]) {
									timeofmachine[machine][0] = Processingtimes[Job][machine];
								}
							}
						}


						vector<vector<int>> setuptimes(Families, vector<int>(Families, 0));
						int setup = 0;
						//Setup Times between Families
						ofile2 << endl << "Stepup times between Families, where [f][f] represents initial setup time" << endl;
						for (int mac = 0; mac < Machines; mac++)
						{
							ofile2 << "On machine\t" << mac + 1 << endl;
							for (int Fam1 = 0; Fam1 < Families; Fam1++)
							{
								for (int Fam2 = 0; Fam2 < Families; Fam2++) {
									setup = rand() % (maxSetupTime - minSetupTime + 1) + minSetupTime;
									ofile2 << setup << "\t";
									setuptimes[Fam1][Fam2] = setup;
								}
								ofile2 << endl;
							}
						}
						int setupmin = 0;
						int setupmax = 0;

						for (int Fam1 = 0; Fam1 < Families; Fam1++)
						{
							for (int Fam2 = 0; Fam2 < Families; Fam2++) {
								if (setuptimes[Fam1][Fam2] > setupmax) {
									setupmax = setuptimes[Fam1][Fam2];
								}
								if (setuptimes[Fam1][Fam2] < setupmin) {
									setupmin = setuptimes[Fam1][Fam2];
								}
							}

						}

						//Maximum of maintenance level of machines
						ofile2 << endl << "Maximum of maintenance level of machines" << endl;
						//加工时间50%-70%

						for (int mac = 0; mac < Machines; mac++)
						{
							int startML = (Families / Factories) * ((maxjobsingroup + minjobsingroup) / 2) * ((timeofmachine[mac][0] + timeofmachine[mac][1]) / 2) * 0.5;
							int endML = (Families / Factories) * ((maxjobsingroup + minjobsingroup) / 2) * ((timeofmachine[mac][0] + timeofmachine[mac][1]) / 2) * 0.75;
							ofile2 << rand() % (endML - startML + 1) + startML << "\t";
						}

						//The maintenance time on machines
						ofile2 << endl << "The maintenance time on machines" << endl;

						for (int mac = 0; mac < Machines; mac++)
						{
							int startMT = setupmax / 2;
							int endMT = timeofmachine[mac][0] + setupmax;
							ofile2 << rand() % (endMT - startMT + 1) + startMT << "\t";
						}

						ofile2.close();
						Counter_Files++;
					}
				}
			}
		}
	}
	ofile1.close();
}

void Problem_PM::ReadInstanceFileNameList(string Dir)
{
	m_InstanceFileNameList.clear();
	m_InstanceNameList.clear();
	ifstream ifile;
	ifile.open(Dir + "\\" + "InstanceFileNameList.txt");

	while (true)
	{
		int x;
		string FName;
		ifile >> x >> FName;
		if (ifile.peek() != EOF)
		{
			m_InstanceFileNameList.push_back(Dir + "\\" + FName);
			m_InstanceNameList.push_back(FName);
		}
		else
			break;
	}
	ifile.close();
}

void Problem_PM::ReadInstance(int InsNo)
{
	ifstream ifile;
	ifile.open(m_InstanceFileNameList[InsNo]);

	string Str;
	//int SetupTimeType;

	//Read Configuration information
	ifile >> Str >> m_Factories;
	ifile >> Str >> m_Machines;
	ifile >> Str >> m_Families;
	ifile >> Str >> m_SetupType;

	//Read Number of Jobs in each Family;
	vector<int> NumbOfJobsInEachFamily(m_Families);
	for (int i = 0; i < 6; i++)
	{
		ifile >> Str;
	}
	for (int Fam = 0; Fam < m_Families; Fam++)
	{
		ifile >> NumbOfJobsInEachFamily[Fam];

	}
	cout << endl;
	//Read Jobs in each Family;
	for (int i = 0; i < 4; i++)
	{
		ifile >> Str;
	}
	m_JobsInEachFamily.clear();
	m_JobsInEachFamily.resize(m_Families);
	for (int Fam = 0; Fam < m_Families; Fam++)
	{
		m_JobsInEachFamily[Fam].resize(NumbOfJobsInEachFamily[Fam]);
	}
	for (int Fam = 0; Fam < m_Families; Fam++)
	{
		for (int i = 0; i < m_JobsInEachFamily[Fam].size(); i++)
		{
			ifile >> m_JobsInEachFamily[Fam][i];
		}
	}

	//Read Total Number of Jobs
	for (int i = 0; i < 4; i++)
	{
		ifile >> Str;
	}
	ifile >> m_Jobs;
	//Read Processing times of Jobs
	for (int i = 0; i < 4; i++)
	{
		ifile >> Str;
	}
	m_JobOperPTime.clear();
	m_JobOperPTime.resize(m_Jobs);
	for (int j = 0; j < m_Jobs; j++)
	{
		m_JobOperPTime[j].resize(m_Machines);
	}
	for (int j = 0; j < m_Jobs; j++)
	{
		for (int m = 0; m < m_Machines; m++)
		{
			ifile >> m_JobOperPTime[j][m];
		}
	}

	//Read Setup times between Families
	for (int i = 0; i < 10; i++)
	{
		ifile >> Str;
	}
	m_SetupTime.clear();
	m_SetupTime.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++)
	{
		m_SetupTime[mac].resize(m_Families);
		for (int Fam = 0; Fam < m_Families; Fam++)
		{
			m_SetupTime[mac][Fam].resize(m_Families);
		}
	}
	for (int mac = 0; mac < m_Machines; mac++)
	{
		for (int i = 0; i < 3; i++)
		{
			ifile >> Str;
		}
		for (int Fam1 = 0; Fam1 < m_Families; Fam1++)
		{
			for (int Fam2 = 0; Fam2 < m_Families; Fam2++)
			{
				ifile >> m_SetupTime[mac][Fam1][Fam2];
			}
		}
	}


	for (int i = 0; i < 6; i++)
	{
		ifile >> Str;
	}
	m_MaxOfML.clear();
	m_MaxOfML.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++)
	{
		ifile >> m_MaxOfML[mac];
	}
	//Read the maintenance time on machines
	for (int i = 0; i < 5; i++)
	{
		ifile >> Str;
	}
	m_MT.clear();
	m_MT.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++)
	{
		ifile >> m_MT[mac];
		//cout << "MT[" << mac << "] = " << MT[mac] << endl;
	}

	m_UnitIdleEC = 1.0;
	m_UnitSetupEC = 0.50;
	m_UnitPMEC = 0.50;

	m_Speed.clear();
	m_Speed.resize(3);
	m_Speed[0] = 1.0;
	m_Speed[1] = 1.5;
	m_Speed[2] = 2.0;

	m_UnitPEC.clear();
	m_UnitPEC.resize(3);

	m_UnitPEC[0] = 4 * m_Speed[0] * m_Speed[0];
	m_UnitPEC[1] = 4 * m_Speed[1] * m_Speed[1];
	m_UnitPEC[2] = 4 * m_Speed[2] * m_Speed[2];

	//变速后的处理时间（真正的处理时间）
	m_TureJobOpertime.clear();
	m_TureJobOpertime.resize(this->m_Jobs, vector<vector<int>>(this->m_Machines, vector<int>(this->m_Speed.size())));

	for (int j = 0; j < m_Jobs; j++)
	{
		for (int i = 0; i < m_Machines; i++)
		{
			for (int s = 0; s < this->m_Speed.size(); ++s)
			{
				m_TureJobOpertime[j][i][s] = static_cast<int> (m_JobOperPTime[j][i] / this->m_Speed[s]);//规定工件的实际加工时间=标准加工时间/速度
			}
		}
	}
}

void Problem_PM::ReadSingleInstance()
{
	ifstream ifile;
	ifile.open("..\\FG_2_3_5_1_1.txt");
	string Str;

	//Read Configuration information
	ifile >> Str >> m_Factories;
	ifile >> Str >> m_Machines;
	ifile >> Str >> m_Families;
	ifile >> Str >> m_SetupType;

	//Read Number of Jobs in each Family;
	vector<int> NumbOfJobsInEachFamily(m_Families);
	for (int i = 0; i < 6; i++)
		ifile >> Str;
	for (int Fam = 0; Fam < m_Families; Fam++) {
		ifile >> NumbOfJobsInEachFamily[Fam];
	}

	//Read Jobs in each Family;
	for (int i = 0; i < 4; i++)
		ifile >> Str;
	m_JobsInEachFamily.clear();
	m_JobsInEachFamily.resize(m_Families);
	for (int Fam = 0; Fam < m_Families; Fam++)
		m_JobsInEachFamily[Fam].resize(NumbOfJobsInEachFamily[Fam]);
	for (int Fam = 0; Fam < m_Families; Fam++) {
		for (int i = 0; i < m_JobsInEachFamily[Fam].size(); i++) {
			ifile >> m_JobsInEachFamily[Fam][i];
		}
	}

	//Read Total Number of Jobs
	for (int i = 0; i < 4; i++)
		ifile >> Str;
	ifile >> m_Jobs;
	//Read Processing times of Jobs
	for (int i = 0; i < 4; i++)
		ifile >> Str;
	m_JobOperPTime.clear();
	m_JobOperPTime.resize(m_Jobs);
	for (int j = 0; j < m_Jobs; j++)
		m_JobOperPTime[j].resize(m_Machines);
	for (int j = 0; j < m_Jobs; j++) {
		for (int m = 0; m < m_Machines; m++) {
			ifile >> m_JobOperPTime[j][m];
		}
	}


	//Read Setup times between Families
	for (int i = 0; i < 10; i++)
		ifile >> Str;
	m_SetupTime.clear();
	m_SetupTime.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++)
	{
		m_SetupTime[mac].resize(m_Families);
		for (int Fam = 0; Fam < m_Families; Fam++)
			m_SetupTime[mac][Fam].resize(m_Families);
	}
	for (int mac = 0; mac < m_Machines; mac++)
	{
		for (int i = 0; i < 3; i++)
			ifile >> Str;
		for (int Fam1 = 0; Fam1 < m_Families; Fam1++) {
			for (int Fam2 = 0; Fam2 < m_Families; Fam2++) {
				ifile >> m_SetupTime[mac][Fam1][Fam2];
			}
		}
	}
	//Read Maximum of maintenance level of machines
	for (int i = 0; i < 6; i++)
	{
		ifile >> Str;
	}
	m_MaxOfML.clear();
	m_MaxOfML.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++) {
		ifile >> m_MaxOfML[mac];
		//cout << "MaxofML[" << mac << "] = " << MaxofML[mac] << endl;
	}
	//Read the maintenance time on machines
	for (int i = 0; i < 5; i++)
		ifile >> Str;
	m_MT.clear();
	m_MT.resize(m_Machines);
	for (int mac = 0; mac < m_Machines; mac++) {
		ifile >> m_MT[mac];
		//cout << "MT[" << mac << "] = " << MT[mac] << endl;
	}

	m_UnitIdleEC = 1.0;
	m_UnitSetupEC = 0.50;
	m_UnitPMEC = 0.50;

	m_Speed.clear();
	m_Speed.resize(3);
	m_Speed[0] = 1.0;
	m_Speed[1] = 1.5;
	m_Speed[2] = 2.0;

	m_UnitPEC.clear();
	m_UnitPEC.resize(3);

	m_UnitPEC[0] = 4 * m_Speed[0] * m_Speed[0];
	m_UnitPEC[1] = 4 * m_Speed[1] * m_Speed[1];
	m_UnitPEC[2] = 4 * m_Speed[2] * m_Speed[2];

	//变速后的处理时间（真正的处理时间）
	m_TureJobOpertime.clear();
	m_TureJobOpertime.resize(this->m_Jobs, vector<vector<int>>(this->m_Machines, vector<int>(this->m_Speed.size())));

	for (int j = 0; j < m_Jobs; j++)
	{
		for (int i = 0; i < m_Machines; i++)
		{
			for (int s = 0; s < this->m_Speed.size(); ++s)
			{
				m_TureJobOpertime[j][i][s] = static_cast<int> (m_JobOperPTime[j][i] / this->m_Speed[s]);//规定工件的实际加工时间=标准加工时间/速度
			}
		}
	}

}

void Problem_PM::GetJobTotalPTime()
{
	this->m_JobTotalPTime.clear();
	this->m_JobTotalPTime.resize(this->m_Jobs, 0);
	m_AllJobTotalPTime = 0;
	for (int j = 0; j < this->m_JobTotalPTime.size(); j++)
	{
		for (int m = 0; m < this->m_Machines; m++)
		{
			this->m_JobTotalPTime[j] += this->m_JobOperPTime[j][m];
			m_AllJobTotalPTime += this->m_JobOperPTime[j][m];
		}
	}
}

//得到组内工件的加工时间之和
void Problem_PM::GetFamTotalPTime()
{
	this->m_FamTotalPTime.clear();
	this->m_FamTotalPTime.resize(this->m_JobsInEachFamily.size(), 0);
	for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
	{
		for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
		{
			this->m_FamTotalPTime[Fam] += this->m_JobTotalPTime[this->m_JobsInEachFamily[Fam][j]];
		}
	}
}

//得到组平均准备时间
void Problem_PM::GetFamAvgSetupTime()
{
	this->m_FamAvgSetupTime.clear();
	this->m_FamAvgSetupTime.resize(this->m_Families, 0);
	//当前组与所有组在所有机器上准备时间的和 / 组数
	for (int CurFam = 0; CurFam < this->m_FamAvgSetupTime.size(); CurFam++)
	{
		for (int PreFam = 0; PreFam < this->m_FamAvgSetupTime.size(); PreFam++)
		{
			for (int m = 0; m < this->m_Machines; m++)
			{
				this->m_FamAvgSetupTime[CurFam] += this->m_SetupTime[m][PreFam][CurFam];
			}
		}
		this->m_FamAvgSetupTime[CurFam] /= this->m_FamAvgSetupTime.size();
	}
}


// implementation of class DetectorSimData
#include "DetectorSimData.h"

void
DetectorSimData::MakeOptDeviceSimData(unsigned int id)
{
	fOptDeviceSimMap.emplace(id, OptDeviceSimData(id));
}

void 
DetectorSimData::SetEnergyDeposit(const double eDep)
{
	fEnergyDeposit.push_back(eDep);
}

void
DetectorSimData::SetEnergyDeposit(Particle::Component comp, const double eDep)
{
	auto it = fEDepComponentMap.emplace(std::make_pair(comp, double()));
	it.first->second.push_back(eDep);
}

void
DetectorSimData::AddBarEnergyDeposit(const int barId, const double eDep)
{
	fCurrentBarEnergyDeposit[barId] += eDep;
}

void
DetectorSimData::AddBarFirstTime(const int barId, const double time)
{
	auto it = fCurrentBarFirstTime.find(barId);
	if (it == fCurrentBarFirstTime.end() || time < it->second) {
		fCurrentBarFirstTime[barId] = time;
	}
}

void
DetectorSimData::EndBarEvent()
{
	fBarEnergyDeposit.push_back(fCurrentBarEnergyDeposit);
	fBarFirstTime.push_back(fCurrentBarFirstTime);
	fCurrentBarEnergyDeposit.clear();
	fCurrentBarFirstTime.clear();
}

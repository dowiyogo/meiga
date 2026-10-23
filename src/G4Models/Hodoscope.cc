#include "Hodoscope.h"
#include "Geometry.h"
#include "Materials.h"
#include "G4MScintillatorBarAction.h"
#include "G4MOptDeviceAction.h"

#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace std;

Hodoscope::Hodoscope(const int id, const Detector::DetectorType type) :
	Scintillator(id, type)
{
	fName = "Hodoscope";
}

void
Hodoscope::BuildDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps)
{
	if (IsNA64Mini()) {
		BuildNA64MiniDetector(logMother, aEvent, fCheckOverlaps);
		return;
	}

	BuildLegacyDetector(logMother, aEvent, fCheckOverlaps);
}

void
Hodoscope::SetDetectorProperties(const boost::property_tree::ptree &aTree)
{
	Scintillator::SetDetectorProperties(aTree);

	for (const auto& property : aTree) {
		string xmlLabel = property.first;
		string xmlValue = property.second.data();

		if (xmlLabel == "hodoscopeModel") {
			boost::algorithm::trim(xmlValue);
			fHodoscopeModel = xmlValue;
			if (IsNA64Mini()) {
				fName = "HodoscopeNA64Mini";
				SetNBars(fBarsPerAxis);
			}
		}
		else if (xmlLabel == "hodoscope_material_detail") {
			boost::algorithm::trim(xmlValue);
			fHodoscopeMaterialDetail = xmlValue;
		}
		else if (xmlLabel == "lower_module_reflector") {
			boost::algorithm::trim(xmlValue);
			fLowerModuleReflector = xmlValue;
		}
		else if (xmlLabel == "upper_module_reflector") {
			boost::algorithm::trim(xmlValue);
			fUpperModuleReflector = xmlValue;
		}
	}
}

double
Hodoscope::GetNA64BarCenterCoordinate(int index, double /*barWidth*/, double subplaneOffset) const
{
	if (index < 0 || index >= fBarsPerAxis) {
		throw std::runtime_error("NA64-mini hodoscope channel index out of range; expected 0..15.");
	}

	const int subplane = index / fBarsPerSubplane;
	const int localIndex = index % fBarsPerSubplane;
	const double centeredOrigin = 0.5 * (fBarsPerSubplane - 1) * fBarPitchNominal;
	return localIndex * fBarPitchNominal + subplane * subplaneOffset - centeredOrigin;
}

void
Hodoscope::AddNA64BarMetadata(const std::string& module, const std::string& axis, int channel, double centerX, double centerY, double centerZ)
{
	const int globalBase = (module == "LOWER") ? 0 : 32;
	const int axisBase = (axis == "X") ? 0 : 16;
	const int globalChannel = globalBase + axisBase + channel;

	BarMetadata meta;
	meta.module = module;
	meta.axis = axis;
	meta.subplane = (channel < fBarsPerSubplane) ? "A" : "B";
	meta.localIndex = channel % fBarsPerSubplane;
	meta.channelIndexAxis = channel;
	meta.globalChannelId = globalChannel;
	meta.sensitiveBarId = globalChannel + 1;
	std::ostringstream barId;
	barId << module << "_" << axis << "_" << std::setw(2) << std::setfill('0') << channel;
	meta.barId = barId.str();
	meta.sipmId = meta.barId + "_SIPM";
	meta.centerX = centerX;
	meta.centerY = centerY;
	meta.centerZ = centerZ;
	meta.length = GetBarLength();
	meta.width = GetBarWidth();
	meta.thickness = GetBarThickness();
	fNA64BarMetadata.push_back(meta);
}

void
Hodoscope::ValidateNA64MiniGeometry() const
{
	if (static_cast<int>(fNA64BarMetadata.size()) != GetNA64TotalActiveBars()) {
		throw std::runtime_error("NA64-mini hodoscope geometry error: total active bars is not 64.");
	}
	if (fBarsPerAxis != 16) {
		throw std::runtime_error("NA64-mini hodoscope geometry error: every axis must have exactly 16 channels.");
	}
	if (fBarsPerAxis * fBarsPerAxis != GetNA64PixelsPerModule()) {
		throw std::runtime_error("NA64-mini hodoscope geometry error: each module must have 256 pixels.");
	}
	if (std::abs((fLowerModuleZMax - fLowerModuleZMin) - fModuleTotalThickness) > 1e-6 * CLHEP::mm ||
	    std::abs((fUpperModuleZMax - fUpperModuleZMin) - fModuleTotalThickness) > 1e-6 * CLHEP::mm) {
		throw std::runtime_error("NA64-mini hodoscope geometry error: a module does not respect z_min/z_max thickness.");
	}
	if (fUpperModuleReflector == "VikuitiESR") {
		throw std::runtime_error("NA64-mini hodoscope geometry error: upper module must not use VikuitiESR by default.");
	}
}

void
Hodoscope::BuildNA64MiniDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps)
{
	if (fHodoscopeMaterialDetail != "active_only") {
		Logger::Print("NA64-mini dead material modes are TODO. Building active_only geometry.", WARNING, "BuildDetector");
	}

	fNA64BarMetadata.clear();
	SetBarLength(75.0 * CLHEP::mm);
	SetBarWidth(3.0 * CLHEP::mm);
	SetBarThickness(1.0 * CLHEP::mm);

	auto sdMan = G4SDManager::GetSDMpointer();
	const int detectorId = GetId();
	G4ThreeVector detectorPos = Geometry::ToG4Vector(GetDetectorPosition(), 1.);

	ostringstream namedetector;
	namedetector << "HodoscopeNA64Mini_" << detectorId;
	Logger::Print("Building NA64-mini hodoscope with 64 active BC-408 bars", INFO, "BuildDetector");

	G4VisAttributes green(G4Colour::Green());
	G4VisAttributes blue(G4Colour::Blue());

	struct ModuleSpec {
		const char* name;
		double zMin;
		double zMax;
	};
	const ModuleSpec modules[] = {
		{"LOWER", fLowerModuleZMin, fLowerModuleZMax},
		{"UPPER", fUpperModuleZMin, fUpperModuleZMax}
	};

	for (const auto& module : modules) {
		for (int axisIndex = 0; axisIndex < 2; ++axisIndex) {
			const string axis = (axisIndex == 0) ? "X" : "Y";
			for (int channel = 0; channel < fBarsPerAxis; ++channel) {
				const int subplane = channel / fBarsPerSubplane;
				const string subplaneName = (subplane == 0) ? "A" : "B";
				const double transverse = GetNA64BarCenterCoordinate(channel, GetBarWidth(), fSubplaneOffset);
				const double localZ =
					(axis == "X")
						? ((subplane == 0) ? 12.5 * CLHEP::mm : 9.5 * CLHEP::mm)
						: ((subplane == 0) ? 3.5 * CLHEP::mm : 0.5 * CLHEP::mm);
				const double centerZ = module.zMin + localZ;
				const double centerX = (axis == "X") ? transverse : 0.0;
				const double centerY = (axis == "Y") ? transverse : 0.0;

				AddNA64BarMetadata(module.name, axis, channel, centerX, centerY, centerZ);
				const BarMetadata& meta = fNA64BarMetadata.back();

				G4Box* solidScinBar = nullptr;
				if (axis == "X") {
					solidScinBar = new G4Box(meta.barId, 0.5 * GetBarWidth(), 0.5 * GetBarLength(), 0.5 * GetBarThickness());
				}
				else {
					solidScinBar = new G4Box(meta.barId, 0.5 * GetBarLength(), 0.5 * GetBarWidth(), 0.5 * GetBarThickness());
				}

				G4LogicalVolume* logScinBar = new G4LogicalVolume(solidScinBar, Materials().ScinPlastic, meta.barId, 0, 0, 0);
				logScinBar->SetVisAttributes(axis == "X" ? green : blue);

				const G4ThreeVector barPosition(
					detectorPos.getX() + centerX,
					detectorPos.getY() + centerY,
					detectorPos.getZ() + centerZ);
				new G4PVPlacement(nullptr, barPosition, logScinBar, meta.barId, logMother, false, meta.sensitiveBarId, fCheckOverlaps);

				G4MScintillatorBarAction* const scinBarSD =
					new G4MScintillatorBarAction(namedetector.str() + "/" + meta.barId, detectorId, meta.sensitiveBarId, aEvent);
				sdMan->AddNewDetector(scinBarSD);
				logScinBar->SetSensitiveDetector(scinBarSD);

				MakeOptDevice(meta.sensitiveBarId, OptDevice::eSiPM);
			}
		}
	}

	ValidateNA64MiniGeometry();
}

void
Hodoscope::BuildLegacyDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps)
{
	// solids
	G4Box* solidEnclosure = nullptr;
	G4Box* solidCasing = nullptr;
	G4Box* solidCoating = nullptr;
	G4Box* solidScinBar = nullptr;

	G4Tubs* solidClad2 = nullptr;
	G4Tubs* solidClad1 = nullptr;
	G4Tubs* solidFiber = nullptr;

	G4Box* solidPixel = nullptr;

	// logical volumes
	G4LogicalVolume* logEnclosure = nullptr;
	G4LogicalVolume* logCasing = nullptr;
	G4LogicalVolume* logCoating = nullptr;
	G4LogicalVolume* logScinBar = nullptr;

	G4LogicalVolume* logClad2 = nullptr;
	G4LogicalVolume* logClad1 = nullptr;
	G4LogicalVolume* logFiber = nullptr;

	G4LogicalVolume* logPixel = nullptr;

	auto sdMan = G4SDManager::GetSDMpointer();

	// number of scintillator bars
	int nBars = GetNBars();
	
	// Scintillator bars
	double barWidth = GetBarWidth();
	double barLength = GetBarLength();
	double barThickness = GetBarThickness();
	double barCoatingThickness = GetBarCoatingThickness();
	double detectorWidth = nBars * (barWidth + barCoatingThickness);
	double halfWidth = 0.5 * nBars * (barWidth + barCoatingThickness);

	// optical fibers
	double claddingThickness = GetCladdingThickness();
	double fiberRadius = GetFiberRadius();

	// optical device
	OptDevice pixel = GetOptDevice(OptDevice::eMChPMT);
	double pixelSizeX = pixel.GetLength();
	double pixelSizeY = pixel.GetWidth();
	double pixelSizeZ = pixel.GetThickness();

	// detector
	G4ThreeVector detectorPos = Geometry::ToG4Vector(GetDetectorPosition(), 1.);
	int detectorId = GetId();
	ostringstream namedetector;
	namedetector.str("");
	namedetector << "Hodoscope_"+to_string(detectorId);
	ostringstream msg; 
	msg << "Building detector: " << namedetector.str() << ", ID: " << detectorId;
	Logger::Print(msg, INFO, "BuildDetector");

	// determine casing size
	double maxLength;
	maxLength = (barLength > detectorWidth ? barLength : detectorWidth);

	const double fCasingThickness = 2 * CLHEP::mm;
	const double fCasingSizeX = maxLength + fCasingThickness;
	const double fCasingSizeY = maxLength + fCasingThickness;
	const double fCasingSizeZ = barThickness * 4;

	const double fEnclosureSizeX = fCasingSizeX;
	const double fEnclosureSizeY = fCasingSizeY;
	const double fEnclosureSizeZ = fCasingSizeZ;

	// Detector enclosure
	solidEnclosure = new G4Box("Enclosure", 0.5*fEnclosureSizeX, 0.5*fEnclosureSizeY, 0.5*fEnclosureSizeZ);
	// Panel casing
	solidCasing = new G4Box("Casing", 0.5*fCasingSizeX, 0.5*fCasingSizeY, 0.5*fCasingSizeZ);
	// Scintillator bars
	solidCoating = new G4Box("BarCoating", 0.5*barLength, 0.5*barWidth, 0.5*barThickness);
	solidScinBar = new G4Box("BarScin", 0.5*(barLength-barCoatingThickness), 0.5*(barWidth-barCoatingThickness), 0.5*(barThickness-barCoatingThickness));
	// Optical Fibers
	solidClad2  = new G4Tubs("Clad2", 0, fiberRadius, barLength/2-2*barCoatingThickness-pixelSizeZ, 0, 360*deg);
	solidClad1  = new G4Tubs("Clad1", 0, fiberRadius - claddingThickness, barLength/2-2*barCoatingThickness-pixelSizeZ, 0, 360*deg);
	solidFiber  = new G4Tubs("Fiber", 0, fiberRadius - 2*claddingThickness, barLength/2-2*barCoatingThickness-pixelSizeZ, 0, 360*deg);
	// Pixel
	solidPixel   = new G4Box("Pixel", pixelSizeX, pixelSizeY, pixelSizeZ);
	
	// Defne BarCoating as LogicalSkingSurface
	new G4LogicalSkinSurface("BarCoating", logCoating, Materials().ScinOptSurf);

	// useful rotation matrix for bars orientation
	G4RotationMatrix* rotationEnclosure = new G4RotationMatrix();
	G4RotationMatrix* rotationFiber = new G4RotationMatrix();
	G4RotationMatrix* rotationBot = new G4RotationMatrix();
	double alpha = GetRotationAngle();
	Logger::Print("Rotation angle: " + std::to_string(alpha), INFO, "BuildDetector");
	rotationEnclosure->rotateX(alpha);
	rotationFiber->rotateY(M_PI/2);
	rotationBot->rotateZ(M_PI/2);

	logEnclosure = new G4LogicalVolume(solidEnclosure, Materials().Air, "Enclosure", 0, 0, 0);
	new G4PVPlacement(rotationEnclosure, G4ThreeVector(detectorPos.getX(), detectorPos.getY(), detectorPos.getZ()), logEnclosure, "Enclosure", logMother, false, detectorId, fCheckOverlaps);

	// color definitions for visualization
	G4VisAttributes green(G4Colour::Green());
	G4VisAttributes blue(G4Colour::Blue());
	G4VisAttributes black(G4Colour::Black());

	int panelId = detectorId+1;

	string nameCasing = "Casing_"+to_string(panelId);

	logCasing = new G4LogicalVolume(solidCasing, Materials().Air, nameCasing, 0, 0, 0);
	new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logCasing, nameCasing, logEnclosure, false, panelId, fCheckOverlaps);

	double fFiberPosX = 0.5*barLength - 2*barCoatingThickness;
	double fFiberPosZ = fiberRadius + claddingThickness;

	// bars of the top panel
	for (int bIt=0; bIt<nBars; ++bIt) {

		string gridName = "X";
		int barId = bIt+1;
		int barId2 = 100*panelId + barId;

		double fBarPosX = 0;
		double fBarPosY = (bIt - nBars/2) * barWidth + barWidth/2;
		// center bar at top of the casing
		double fBarPosZ = 0.5*(0.5*barThickness*2 + fCasingThickness);
		
		G4ThreeVector barPosition = G4ThreeVector(fBarPosX, fBarPosY, fBarPosZ);
		
		string nameCoating = "BarCoating_"+gridName;
		string nameScinBar = "BarScin_"+gridName;
		string nameClad2 = "FiberClad2_"+gridName;
		string nameClad1 = "FiberClad1_"+gridName;
		string nameFiber = "Fiber_"+gridName;
		string namePixel = "Pixel_"+gridName+to_string(barId2);
		// register pixel in the detector class
		MakeOptDevice(barId, OptDevice::eMChPMT);

		// logical volumes
		logCoating = new G4LogicalVolume(solidCoating, Materials().ScinCoating, nameCoating, 0, 0, 0);
		logScinBar = new G4LogicalVolume(solidScinBar, Materials().ScinPlastic, nameScinBar, 0, 0, 0);
		logClad2 = new G4LogicalVolume(solidClad2, Materials().FPethylene, nameClad2, 0, 0, 0);
		logClad2->SetVisAttributes(green);
		logClad1 = new G4LogicalVolume(solidClad1, Materials().Pethylene, nameClad1, 0, 0, 0);
		logClad1->SetVisAttributes(green);
		logFiber = new G4LogicalVolume(solidFiber, Materials().PMMA, nameFiber, 0, 0, 0);
		logFiber->SetVisAttributes(green);
		logPixel = new G4LogicalVolume(solidPixel, Materials().Pyrex, namePixel, 0, 0, 0);
		logPixel->SetVisAttributes(blue);

		// physical volumes
		new G4PVPlacement(nullptr, barPosition, logCoating, nameCoating, logCasing, false, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(), logScinBar, nameScinBar, logCoating, false, barId2, fCheckOverlaps);
		// register scintillator bar as Sensitive Detector to store hits
		G4MScintillatorBarAction* const scinBarSD = new G4MScintillatorBarAction(namedetector.str() + "/" + nameScinBar+"_"+to_string(barId), detectorId, barId, aEvent);
		sdMan->AddNewDetector(scinBarSD);
		logScinBar->SetSensitiveDetector(scinBarSD);

		new G4PVPlacement(rotationFiber, G4ThreeVector(0., 0., fFiberPosZ), logClad2, nameClad2, logScinBar, true, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logClad1, nameClad1, logClad2, false, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logFiber, nameFiber, logClad1, false, barId2, fCheckOverlaps);

		new G4PVPlacement(rotationFiber, G4ThreeVector(fFiberPosX, 0, fFiberPosZ), logPixel, namePixel, logScinBar, false, barId2, fCheckOverlaps);

		// registration of pixel as a Sensitive Volume
		G4MOptDeviceAction* const pixelTopSD = new G4MOptDeviceAction(namedetector.str() + "/" + namePixel, detectorId, barId, aEvent);
		sdMan->AddNewDetector(pixelTopSD);
		logPixel->SetSensitiveDetector(pixelTopSD);

	} // end loop of bars of the top panel

	// bars of the bottom panel
	for (int bIt=0; bIt<nBars; ++bIt) {

		string gridName = "Y";
		int barId = bIt+1+nBars;
		int barId2 = 100*panelId + barId;

		double fBarPosX = (bIt - nBars/2) * barWidth + barWidth/2;
		double fBarPosY = 0;
		// center bar at top of the casing
		double fBarPosZ = -0.5*(0.5*barThickness*2 + fCasingThickness);

		G4ThreeVector barPosition = G4ThreeVector(fBarPosX, fBarPosY, fBarPosZ);

		string nameCoating = "BarCoating_"+gridName;
		string nameScinBar = "BarScin_"+gridName;
		string nameClad2 = "FiberClad2_"+gridName;
		string nameClad1 = "FiberClad1_"+gridName;
		string nameFiber = "Fiber_"+gridName;
		string namePixel = "Pixel_"+gridName+to_string(barId2);
		// register pixel in the detector class
		MakeOptDevice(barId, OptDevice::eMChPMT);

		// logical volumes
		logCoating = new G4LogicalVolume(solidCoating, Materials().ScinCoating, nameCoating, 0, 0, 0);
		logScinBar = new G4LogicalVolume(solidScinBar, Materials().ScinPlastic, nameScinBar, 0, 0, 0);
		logClad2 = new G4LogicalVolume(solidClad2, Materials().FPethylene, nameClad2, 0, 0, 0);
		logClad2->SetVisAttributes(green);
		logClad1 = new G4LogicalVolume(solidClad1, Materials().Pethylene, nameClad1, 0, 0, 0);
		logClad1->SetVisAttributes(green);
		logFiber = new G4LogicalVolume(solidFiber, Materials().PMMA, nameFiber, 0, 0, 0);
		logFiber->SetVisAttributes(green);
		logPixel = new G4LogicalVolume(solidPixel, Materials().Pyrex, namePixel, 0, 0, 0);
		logPixel->SetVisAttributes(blue);

		// physical volumes
		new G4PVPlacement(rotationBot, barPosition, logCoating, nameCoating, logCasing, false, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(), logScinBar, nameScinBar, logCoating, false, barId2, fCheckOverlaps);
		
		G4MScintillatorBarAction* const scinBarSD = new G4MScintillatorBarAction(namedetector.str() + "/" + nameScinBar+"_"+to_string(barId), detectorId, barId, aEvent);
		sdMan->AddNewDetector(scinBarSD);
		logScinBar->SetSensitiveDetector(scinBarSD);
		// register scintillator bar as Sensitive Detector to store hits
		new G4PVPlacement(rotationFiber, G4ThreeVector(0, 0, fFiberPosZ), logClad2, nameClad2, logScinBar, true, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logClad1, nameClad1, logClad2, false, barId2, fCheckOverlaps);
		new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logFiber, nameFiber, logClad1, false, barId2, fCheckOverlaps);

		new G4PVPlacement(rotationFiber, G4ThreeVector(fFiberPosX, 0, fFiberPosZ), logPixel, namePixel, logScinBar, false, barId2, fCheckOverlaps);

		// registration of pixel as a Sensitive Volume
		G4MOptDeviceAction* const pixelBotSD = new G4MOptDeviceAction(namedetector.str() + "/" + namePixel, detectorId, barId, aEvent);
		sdMan->AddNewDetector(pixelBotSD);
		logPixel->SetSensitiveDetector(pixelBotSD);


	} // end loop of batrs of the bottom panel

}

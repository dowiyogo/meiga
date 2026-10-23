#ifndef Hodoscope_h
#define Hodoscope_h 1

#include "G4LogicalVolume.hh"

#include "Event.h"
#include "Detector.h"
#include "Scintillator.h"

#include <string>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4OpticalSkinSurface;

class Hodoscope : public Scintillator
{
	/*
		Model class for an hodoscope-like detector.
		This type of detectors are used in muon detection techniques.
		
		Consists of a given number of planes composed by scintillator bars.
		Scintillator bars are displayed along the x- and y-axis conforming
		grids of (2N_x - 1) * (2N_y - 1) number of pixels per plane.


	*/

public:
	struct BarMetadata {
		std::string module;
		std::string axis;
		std::string subplane;
		std::string barId;
		std::string sipmId;
		int localIndex = 0;
		int channelIndexAxis = 0;
		int globalChannelId = 0;
		int sensitiveBarId = 0;
		double centerX = 0;
		double centerY = 0;
		double centerZ = 0;
		double length = 0;
		double width = 0;
		double thickness = 0;
	};

	Hodoscope(const int id, const Detector::DetectorType type);
	void BuildDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps = false) override;
	void SetDetectorProperties(const boost::property_tree::ptree &aTree) override;

	bool IsNA64Mini() const { return fHodoscopeModel == "NA64Mini"; }
	int GetNA64TotalActiveBars() const { return 64; }
	int GetNA64PixelsPerModule() const { return 256; }
	const std::vector<BarMetadata>& GetNA64BarMetadata() const { return fNA64BarMetadata; }
	double GetNA64BarCenterCoordinate(int index, double barWidth, double subplaneOffset) const;

private:
	void BuildLegacyDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps);
	void BuildNA64MiniDetector(G4LogicalVolume* logMother, Event& aEvent, G4bool fCheckOverlaps);
	void ValidateNA64MiniGeometry() const;
	void AddNA64BarMetadata(const std::string& module, const std::string& axis, int channel, double centerX, double centerY, double centerZ);

	std::string fHodoscopeModel = "Legacy";
	std::string fHodoscopeMaterialDetail = "active_only";
	std::string fLowerModuleReflector = "VikuitiESR";
	std::string fUpperModuleReflector = "TiO2Paint";
	std::vector<BarMetadata> fNA64BarMetadata;

	double fLowerModuleZMin = 0.00 * CLHEP::mm;
	double fLowerModuleZMax = 13.00 * CLHEP::mm;
	double fUpperModuleZMin = 98.40 * CLHEP::mm;
	double fUpperModuleZMax = 111.40 * CLHEP::mm;
	double fModuleTotalThickness = 13.00 * CLHEP::mm;
	double fModuleCenterSeparation = 98.40 * CLHEP::mm;
	double fXYAxisSeparation = 6.0 * CLHEP::mm;
	double fSubplaneOffset = 2.0 * CLHEP::mm;
	double fBarPitchNominal = 4.0 * CLHEP::mm;
	int fBarsPerAxis = 16;
	int fBarsPerSubplane = 8;
};

#endif

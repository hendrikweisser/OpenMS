// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMINFILE_H
#define OPENMS_FORMAT_XTANDEMINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
	/**
		@brief XTandem input file adapter
		
		This class is able to create a X!Tandem configuration file for a search
  	  	
  	@ingroup FileIO
	*/
  class XTandemInfile
  {
    public:

			enum ERROR_UNIT
			{
				DALTONS = 0,
				PPM
			};
						
			enum MASS_TYPE
			{
				MONOISOTOPIC = 0,
				AVERAGE
			};
			
			/// constructor
			XTandemInfile();

			/// constructor
			virtual ~XTandemInfile();

			//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
			//<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
			//<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
			//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
			//<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
			//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
			//<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
			//<note>values are monoisotopic|average </note>
			void setPeakMassTolerance(double tolerance);

			double getPeakMassTolerance() const;

			void setPrecursorMassTolerancePlus(double tol);

			double getPrecursorMassTolerancePlus() const;

			void setPrecursorMassToleranceMinus(double tol);

			double getPrecursorMassToleranceMinus() const;

			void setPrecursorMonoisotopicError(MASS_TYPE mono_isotopic);

			MASS_TYPE getPrecursorMonoisotopicError() const;

			void setPeakMassErrorUnit(ERROR_UNIT unit);

			ERROR_UNIT getPeakMassErrorUnit() const;

			void setPrecursorMassErrorUnit(ERROR_UNIT unit);

			ERROR_UNIT getPrecursorMassErrorUnit() const;
	

			//<note>spectrum conditioning parameters</note>
			//<note type="input" label="spectrum, dynamic range">100.0</note>
			//<note>The peaks read in are normalized so that the most intense peak
			//is set to the dynamic range value. All peaks with values of less that
			//1, using this normalization, are not used. This normalization has the
			//overall effect of setting a threshold value for peak intensities.</note>
			//<note type="input" label="spectrum, total peaks">50</note> 
			//<note>If this value is 0, it is ignored. If it is greater than zero (lets say 50),
			//then the number of peaks in the spectrum with be limited to the 50 most intense
			//peaks in the spectrum. X! tandem does not do any peak finding: it only
			//limits the peaks used by this parameter, and the dynamic range parameter.</note>
			//<note type="input" label="spectrum, maximum parent charge">4</note>
			//<note type="input" label="spectrum, use noise suppression">yes</note>
			//<note type="input" label="spectrum, minimum parent m+h">500.0</note>
			//<note type="input" label="spectrum, minimum fragment mz">150.0</note>
			//<note type="input" label="spectrum, minimum peaks">15</note> 
			//<note type="input" label="spectrum, threads">1</note>
			//<note type="input" label="spectrum, sequence batch size">1000</note>

			void setDynamicRange(double range);

			double getDynamicRange() const;

			void setTotalNumberOfPeaks(UInt number);

			UInt getTotalNumberOfPeaks() const;

			void setMaximumPrecursorCharge(UInt charge);

			UInt getMaximumPrecursorCharge() const;

			void setUseNoiseSupression(bool use_supression = true);

			bool getUseNoiseSupression() const;
			
			void setPrecursorLowerMZ(double mz);

			double getPrecursorLowerMZ() const;

			void setPeakLowerMZ(double mz);

			double getPeakLowerMZ() const;

			void setMinimalNumberOfPeaks(UInt number);

			UInt getMinimalNumberOfPeaks() const;

			void setNumberOfThreads(UInt threads);

			UInt getNumberOfThreads() const;

			void setBatchSize(UInt size);

			UInt getBatchSize();
		
			//<note>residue modification parameters</note>
			//<note type="input" label="residue, modification mass">57.022@C</note>
			//<note>The format of this parameter is m@X, where m is the modfication
			//mass in Daltons and X is the appropriate residue to modify. Lists of
			//modifications are separated by commas. For example, to modify M and C
			//with the addition of 16.0 Daltons, the parameter line would be
			//+16.0@M,+16.0@C
			//Positive and negative values are allowed.
			//</note>
			//<note type="input" label="residue, potential modification mass"></note>
			//<note>The format of this parameter is the same as the format
			//for residue, modification mass (see above).</note>
			//<note type="input" label="residue, potential modification motif"></note>
			//<note>The format of this parameter is similar to residue, modification mass,
			//with the addition of a modified PROSITE notation sequence motif specification.
			//For example, a value of 80@[ST!]PX[KR] indicates a modification
			//of either S or T when followed by P, and residue and the a K or an R.
			//A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it
			//is NOT followed by a P, then either an S or a T, NOT followed by a P.
			//Positive and negative values are allowed.
			//</note>
			void setFixedModifications(const String& mods);

			const String& getFixedMofications() const;

			void setVariableModifications(const String& mods);

			const String& getVariableModifications() const;

			void setVariableModificationMotif(const String& motif);

			const String& getVariableModificationMotif() const;
			
			void setOutputFilename(const String& output);

			const String& getOutputFilename() const;
	
			void setInputFilename(const String& input_file);

			const String& getInputFilename() const;

			void write(const String& filename) throw (Exception::UnableToCreateFile);

			void load(const String& filename) throw (Exception::FileNotFound, Exception::ParseError);

    protected:

			void writeTo_(std::ostream& os);

			void writeNote_(std::ostream& os, const String& type, const String& label, const String& value);

			double peak_mass_tolerance_;

			double precursor_mass_tolerance_plus_;

			double precursor_mass_tolerance_minus_;

			MASS_TYPE precursor_monoisotopic_error_;

			ERROR_UNIT precursor_mass_error_unit_;

			ERROR_UNIT peak_mass_error_unit_;

			double dynamic_range_;

			UInt total_number_peaks_;

			UInt max_precursor_charge_;
			
			bool noise_supression_;

			double precursor_lower_mz_;

			double peak_lower_mz_;

			UInt min_number_peaks_;

			UInt number_of_threads_;

			UInt batch_size_;
			
			String fixed_modifications_;

			String variable_modifications_;

			String variable_modification_motif_;

			String input_filename_;

			String output_filename_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTINFILE_H

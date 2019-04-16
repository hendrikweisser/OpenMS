// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_StandalonePrecursorCorrector StandalonePrecursorCorrector

    @brief Perform MS2 precursor mass correction and charge determination

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_StandalonePrecursorCorrector.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_StandalonePrecursorCorrector.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPStandalonePrecursorCorrector :
  public TOPPBase
{
public:

  TOPPStandalonePrecursorCorrector() :
    TOPPBase("StandalonePrecursorCorrector", "Perform MS2 precursor mass correction and charge determination", false)
  {
  }

protected:
  stringstream debug_stream_;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file containing centroided LC-MS/MS data");
    setValidFormats_("in", vector<String>(1, "mzML"));
    registerOutputFile_("out", "<file>", "", "Output file containing LC-MS/MS data with corrected precursor masses/charges", false);
    setValidFormats_("out", vector<String>(1, "mzML"));

    registerOutputFile_("qc_out", "<file>", "", "Output file containing spectra for quality control", false);
    setValidFormats_("qc_out", vector<String>(1, "mzML"));
    registerOutputFile_("qc_id_out", "<file>", "", "Output file containing annotations for quality control", false);
    setValidFormats_("qc_id_out", vector<String>(1, "idXML"));

    registerStringOption_("molecule_type", "<choice>", "protein", "Type of molecule that was analyzed in the experiment", false);
    setValidStrings_("molecule_type", ListUtils::create<String>("protein,RNA"));

    registerIntOption_("assign_charge", "<choice>", 0, "Assign charge states of spectra? Options: -1 - never, 0 - only if charge is 0 (unknown), 1 - always", false);
    setMinInt_("assign_charge", -1);
    setMaxInt_("assign_charge", 1);
    registerIntOption_("min_charge", "<num>", 2, "Minimum charge to consider (if 'assign_charge' is 0 or 1)", false);
    setMinInt_("min_charge", 1);
    registerIntOption_("max_charge", "<num>", 6, "Maximum charge to consider (if 'assign_charge' is 0 or 1)", false);
    setMinInt_("max_charge", 1);

    registerIntOption_("n_isotopes", "<num>", 5, "Number of isotopologue peaks to consider", false);
    setMinInt_("n_isotopes", 2);

    registerDoubleOption_("mz_tolerance", "<tolerance>", 10.0, "m/z tolerance when matching peaks", false);
    setMinFloat_("mz_tolerance", 0.0);
    registerStringOption_("mz_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance", false);
    setValidStrings_("mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
  }


  vector<double> score_(const MSSpectrum& real_spec,
                        const MSSpectrum& template_spec, Int charge,
                        double prec_mz, double mz_tol, bool tol_unit_ppm,
                        PeptideIdentification* qc_id_ptr)
  {
    Size n_isotopes = template_spec.size();
    vector<double> result(n_isotopes, -1.0);
    if (tol_unit_ppm)
    {
      mz_tol = prec_mz * mz_tol * 1e-6;
    }
    // extract all relevant peak intensities from real spectrum:
    vector<double> intensities;
    vector<double> peak_mz; // for QC annotations
    intensities.reserve(n_isotopes * 2 - 1);
    peak_mz.reserve(n_isotopes * 2 - 1);
    // shift template spectrum to the left (precursor peak is at the right):
    double shift = template_spec.back().getMZ() - prec_mz;
    for (Size i = 0; i < n_isotopes - 1; ++i)
    {
      Int index = real_spec.findNearest(template_spec[i].getMZ() - shift,
                                        mz_tol);
      if (index >= 0)
      {
        intensities.push_back(real_spec[index].getIntensity());
        peak_mz.push_back(real_spec[index].getMZ());
      }
      else
      {
        intensities.push_back(0.0);
        peak_mz.push_back(0.0);
      }
    }
    // shift template spectrum such that precursor peak is at the left:
    shift = template_spec[0].getMZ() - prec_mz;
    for (Size i = 0; i < n_isotopes; ++i)
    {
      Int index = real_spec.findNearest(template_spec[i].getMZ() - shift,
                                        mz_tol);
      if (index >= 0)
      {
        intensities.push_back(real_spec[index].getIntensity());
        peak_mz.push_back(real_spec[index].getMZ());
      }
      else
      {
        intensities.push_back(0.0);
        peak_mz.push_back(0.0);
      }
    }
    Size peak_count = count_if(intensities.begin(), intensities.end(),
                               [](double x){return x > 0.0;});
    if (peak_count <= 1) // not enough peaks to do anything useful
    {
      return result;
    }
    vector<double> template_intensities;
    template_intensities.reserve(n_isotopes);
    for (const auto& peak : template_spec)
    {
      template_intensities.push_back(peak.getIntensity());
    }
    debug_stream_ << "Charge " << charge <<  ":\nReal intensities: "
                  << ListUtils::concatenate(intensities, ", ")
                  << "\nIso. intensities: "
                  << ListUtils::concatenate(template_intensities, ", ") << endl;

    Size start_pos = n_isotopes - 1; // monoisotopic pos. in "intensities"
    for (Int offset = 0; offset < n_isotopes; ++offset)
    {
      auto start_it = intensities.begin() + start_pos;
      double cor = Math::pearsonCorrelationCoefficient(
        template_intensities.begin(), template_intensities.end(), start_it,
        start_it + n_isotopes);
      if (cor != cor) cor = -1.0; // NaN
      result[offset] = cor;
      if (qc_id_ptr) // @TODO: check for "valid" correlation?
      {
        vector<PeptideHit::PeakAnnotation> annotations;
        for (Size i = 0; i < n_isotopes; ++i)
        {
          double mz = peak_mz[start_pos + i];
          if (mz == 0.0) continue; // no matching peak
          PeptideHit::PeakAnnotation ann;
          ann.annotation = "i" + String(i);
          ann.charge = charge;
          ann.mz = mz;
          ann.intensity = intensities[start_pos + i];
          annotations.push_back(ann);
        }
        if (!annotations.empty())
        {
          PeptideHit hit;
          hit.setCharge(charge);
          hit.setMetaValue("label", "z=" + String(charge) + "/i=" +
                           String(offset));
          hit.setScore(cor);
          hit.setPeakAnnotations(annotations);
          qc_id_ptr->insertHit(hit);
        }
      }
      --start_pos;
    }
    debug_stream_ << "Correlations: " << ListUtils::concatenate(result, ", ")
                  << endl;

    return result;
  }

  
  ExitCodes main_(int, const char**) override
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in = getStringOption_("in"), out = getStringOption_("out");

    MSExperiment data;
    MzMLFile().load(in, data);
    data.sortSpectra();

    // build map from fragment spectrum (index) to precursor spectrum (index):
    map<Size, Size> precursor_indexes;
    Size previous_precursor = Size(-1);
    for (Size i = 0; i < data.size(); ++i)
    {
      if (data[i].getMSLevel() < 2)
      {
        previous_precursor = i;
      }
      else if ((data[i].getMSLevel() == 2) && // @TODO: include MS3 etc.?
               (!data[i].getPrecursors().empty()) &&
               (previous_precursor != Size(-1)))
      {
        precursor_indexes[i] = previous_precursor;
      }
    }

    Int assign_charge = getIntOption_("assign_charge");
    Int min_charge = getIntOption_("min_charge");
    Int max_charge = getIntOption_("max_charge");
    if (min_charge > max_charge) swap(min_charge, max_charge);
    vector<Int> charges;
    charges.reserve(max_charge - min_charge + 1);
    for (Int z = min_charge; z <= max_charge; ++z)
    {
      charges.push_back(z);
    }
    
    Int n_isotopes = getIntOption_("n_isotopes");
    String molecule_type = getStringOption_("molecule_type");
    CoarseIsotopePatternGenerator iso_gen(n_isotopes, false);
    auto averagine_fun = (molecule_type == "protein") ?
      &CoarseIsotopePatternGenerator::estimateFromPeptideWeight :
      &CoarseIsotopePatternGenerator::estimateFromRNAWeight;

    double mz_tol = getDoubleOption_("mz_tolerance");
    bool tol_unit_ppm = getStringOption_("mz_tolerance_unit") == "ppm";

    String qc_out = getStringOption_("qc_out");
    String qc_id_out = getStringOption_("qc_id_out");
    MSExperiment qc_spectra;
    vector<PeptideIdentification> qc_ids;
    
    Size counter = 0;
    progresslogger.startProgress(0, precursor_indexes.size(),
                                 "Correcting precursor information...");
    for (auto index_pair : precursor_indexes)
    {
      debug_stream_.str(""); // clear content
      // @TODO: what if there's more than one precursor?
      const Precursor& prec = data[index_pair.first].getPrecursors()[0];
      double prec_mz = prec.getMZ(); // @TODO: check where the peak actually is?
      Int prec_charge = prec.getCharge();
      const MSSpectrum& real_spec = data[index_pair.second];
      vector<Int> charge_candidates;
      if ((assign_charge == 1) || ((assign_charge == 0) && (prec_charge == 0)))
      {
        charge_candidates = charges;
      }
      else
      {
        charge_candidates.push_back(prec_charge);
      }
      if (charge_candidates[0] == 0) continue; // keep unknown charge

      if (!qc_out.empty()) // store QC data
      {
        MSSpectrum qc_spec;
        double mz_min = prec_mz - n_isotopes / min_charge;
        double mz_max = prec_mz + n_isotopes / min_charge;
        double mz_range = mz_max - mz_min;
        for (auto it = real_spec.MZBegin(mz_min - 0.1 * mz_range);
             it != real_spec.MZEnd(mz_max + 0.1 * mz_range); ++it)
        {
          if (it->getIntensity() > 0) qc_spec.push_back(*it);
        }
        if (qc_spec.empty()) continue; // no data
        qc_spec.setMSLevel(2);
        qc_spec.setRT(data[index_pair.first].getRT());
        qc_spec.getPrecursors().push_back(prec);
        qc_spectra.addSpectrum(qc_spec);
      }
      
      vector<pair<Size, double>> best_results; // best offset/score per charge
      best_results.reserve(charge_candidates.size());
      for (Int charge : charge_candidates)
      {
        double mass = prec_mz * charge;
        IsotopeDistribution iso_dist = (iso_gen.*averagine_fun)(mass);
        MSSpectrum template_spec;
        for (auto peak : iso_dist)
        {
          peak.setMZ(peak.getMZ() / charge);
          template_spec.push_back(peak);
        }
        PeptideIdentification* qc_id_ptr = nullptr;
        if (!qc_id_out.empty()) // add QC annotations
        {
          template_spec.getIntegerDataArrays().resize(1);
          template_spec.getIntegerDataArrays()[0].resize(iso_dist.size(),
                                                         charge);
          template_spec.getStringDataArrays().resize(1);
          for (Size i = 0; i < iso_dist.size(); ++i)
          {
            template_spec.getStringDataArrays()[0].push_back("i" + String(i));
          }
          PeptideIdentification qc_id;
          qc_id.setIdentifier("id");
          qc_id.setRT(data[index_pair.first].getRT());
          qc_id.setMZ(prec_mz);
          qc_id.setScoreType("correlation");
          qc_ids.push_back(qc_id);
          qc_id_ptr = &qc_ids.back();
        }
        vector<double> cors = score_(real_spec, template_spec, charge, prec_mz,
                                     mz_tol, tol_unit_ppm, qc_id_ptr);
        best_results.push_back(make_pair(0, cors[0]));
        for (Size offset = 1; offset < cors.size(); ++offset)
        {
          if (cors[offset] > best_results.back().second)
          {
            best_results.back().first = offset;
            best_results.back().second = cors[offset];
          }
        }
      }
      // find overall best match:
      Int best_charge = charge_candidates[0];
      Size best_offset = best_results[0].first;
      double best_score = best_results[0].second;
      for (Size i = 1; i < charge_candidates.size(); ++i)
      {
        if (best_results[i].second > best_score)
        {
          best_charge = charge_candidates[i];
          best_offset = best_results[i].first;
          best_score = best_results[i].second;
        }
      }
      if ((best_score > -1.0) && ((best_charge != prec_charge) ||
                                  (best_offset != 0)))
      {
        ++counter;
        LOG_INFO << "Spectrum " << index_pair.first << " (RT: "
                 << data[index_pair.first].getRT() << "): assigned charge "
                 << prec_charge << ", estimated charge " << best_charge
                 << " (isotopologue offset " << best_offset << ")" << endl;
        LOG_DEBUG << debug_stream_.str() << endl;
      }
      progresslogger.nextProgress();
    }
    progresslogger.endProgress();

    LOG_INFO << "Found " << counter << " suggestions." << endl;

    if (!qc_out.empty())
    {
      MzMLFile().store(qc_out, qc_spectra);
    }
    if (!qc_id_out.empty())
    {
      vector<ProteinIdentification> proteins(1);
      proteins[0].setIdentifier("id");
      IdXMLFile().store(qc_id_out, proteins, qc_ids);
    }
    
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPStandalonePrecursorCorrector tool;
  return tool.main(argc, argv);
}

/// @endcond

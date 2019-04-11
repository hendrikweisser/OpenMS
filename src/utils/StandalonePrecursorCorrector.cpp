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
    registerIntOption_("max_charge", "<num>", 5, "Maximum charge to consider (if 'assign_charge' is 0 or 1)", false);
    setMinInt_("max_charge", 1);

    registerIntOption_("n_isotopes", "<num>", 6, "Number of isotopologue peaks to consider", false);
    setMinInt_("n_isotopes", 2);

    registerDoubleOption_("mz_tolerance", "<tolerance>", 10.0, "m/z tolerance when matching peaks", false);
    setMinFloat_("mz_tolerance", 0.0);
    registerStringOption_("mz_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance", false);
    setValidStrings_("mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
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
      
      vector<pair<Size, double>> best_scores;
      best_scores.reserve(charge_candidates.size());
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
        }
        vector<PeptideHit::PeakAnnotation> annotations;
        Size best_offset = 0;
        double best_score = 0.0;
        if (!qc_id_out.empty())
        {
          PeptideIdentification qc_id;
          qc_id.setIdentifier("id");
          qc_id.setRT(data[index_pair.first].getRT());
          qc_id.setMZ(prec_mz);
          qc_id.setScoreType("hyperscore");
          qc_ids.push_back(qc_id);
        }
        for (Size offset = 0; offset < n_isotopes; ++offset)
        {
          // move template spectrum such that precursor m/z coincides with an
          // isotopologue peak (given by offset):
          double diff = template_spec[offset].getMZ() - prec_mz;
          for (auto& peak : template_spec)
          {
            peak.setMZ(peak.getMZ() - diff);
          }

          // @TODO: implement a better scoring function that includes m/z
          // deviations and correlation with averagine intensities
          double score = MetaboliteSpectralMatching::computeHyperScore(
            mz_tol, tol_unit_ppm, real_spec, template_spec, annotations);
          if (!qc_id_out.empty() && (score > 0))
          {
            PeptideIdentification& qc_id = qc_ids.back();
            PeptideHit hit;
            hit.setCharge(charge);
            hit.setMetaValue("label", "z=" + String(charge) + "/i=" +
                             String(offset));
            hit.setScore(score);
            hit.setPeakAnnotations(annotations);
            qc_id.insertHit(hit);
          }
          annotations.clear();
          if (score > best_score)
          {
            best_score = score;
            best_offset = offset;
          }          
        }
        best_scores.push_back(make_pair(best_offset, best_score));
      }
      // find overall best match:
      Int best_charge = charge_candidates[0];
      Size best_offset = best_scores[0].first;
      double best_score = best_scores[0].second;
      for (Size i = 1; i < charge_candidates.size(); ++i)
      {
        if (best_scores[i].second > best_score)
        {
          best_charge = charge_candidates[i];
          best_offset = best_scores[i].first;
          best_score = best_scores[i].second;
        }
      }
      if ((best_score > 0) && ((best_charge != prec_charge) ||
                               (best_offset != 0)))
      {
        ++counter;
        LOG_INFO << "Spectrum " << index_pair.first << " (RT: "
                 << data[index_pair.first].getRT() << "): assigned charge "
                 << prec_charge << ", estimated charge " << best_charge
                 << " (isotopologue offset " << best_offset << ")" << endl;
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

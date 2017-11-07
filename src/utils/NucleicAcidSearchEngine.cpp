// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Samuel Wein, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/IdentificationData.h>

// file types
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>

// digestion enzymes
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>

// ribonucleotides
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

// spectra comparison
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <boost/regex.hpp>

#include <QtCore/QProcess>

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>


#ifdef _OPENMP
#include <omp.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//#define DEBUG_NASEARCH

using namespace OpenMS;
using namespace std;

class NucleicAcidSearchEngine :
  public TOPPBase
{
  using ConstRibonucleotidePtr = const Ribonucleotide*;

public:
  NucleicAcidSearchEngine() :
    TOPPBase("NucleicAcidSearchEngine", "Annotate nucleic acid identifications to MS/MS spectra.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file: spectra");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("database", "<file>", "", "Input file: sequence database");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "Output file: mzTab");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_("id_out", "<file>", "", "Output file: idXML (for visualization in TOPPView)");
    setValidFormats_("id_out", ListUtils::create<String>("idXML"));

    registerOutputFile_("theo_ms2_out", "<file>", "", "Output file: theoretical MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("theo_ms2_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("exp_ms2_out", "<file>", "", "Output file: experimental MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("exp_ms2_out", ListUtils::create<String>("mzML"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance (+/- around precursor m/z)", false);

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerIntOption_("precursor:min_charge", "<num>", -1, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", -9, "Maximum precursor charge to be considered.", false, false);

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)", false, false);

    registerTOPPSubsection_("fragment", "Fragment (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerTOPPSubsection_("modifications", "Modifications Options");

    // add modified ribos from database
    vector<String> all_mods;
    for (auto r : *RibonucleotideDB::getInstance())
    {
      if (r->isModified())
      {
        String code = r->getCode();
        // commas aren't allowed in parameter string restrictions:
        all_mods.push_back(code.substitute(',', '_'));
      }
    }

    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_oligo", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate oligo", false, false);

    registerTOPPSubsection_("oligo", "Oligonucleotide Options");
    registerIntOption_("oligo:min_size", "<num>", 5, "Minimum size an oligonucleotide must have after digestion to be considered in the search.", false);
    registerIntOption_("oligo:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

    StringList all_enzymes;
    RNaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("oligo:enzyme", "<choice>", "no cleavage", "The enzyme used for RNA digestion.", false);
    setValidStrings_("oligo:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);
  }


  // slimmer structure to store basic hit information
  struct AnnotatedHit
  {
    String sequence;
    SignedSize mod_index; // enumeration index of the modification
    double score; // the score
    vector<PeptideHit::PeakAnnotation> annotations; // peak/ion annotations

    static bool hasBetterScore(const AnnotatedHit& a, const AnnotatedHit& b)
    {
      return a.score > b.score;
    }
  };


  // query modified residues from database
  vector<ConstRibonucleotidePtr> getModifications_(const set<String>& mod_names)
  {
    vector<ConstRibonucleotidePtr> modifications;
    for (String m : mod_names)
    {
      m.substitute('_', ',');
      ConstRibonucleotidePtr rm = RibonucleotideDB::getInstance()->getRibonucleotide(m);
      modifications.push_back(rm);
    }
    return modifications;
  }


  // check for minimum size
  class HasInvalidLength
  {
    Size min_size_;
  public:
    explicit HasInvalidLength(Size min_size)
      : min_size_(min_size)
    {
    }
    bool operator()(const NASequence& s) { return s.size() < min_size_; }
  };


  // spectrum must not contain 0 intensity peaks and must be sorted by m/z
  void deisotopeAndSingleChargeMSSpectrum_(
    MSSpectrum& in,
    Int min_charge,
    Int max_charge,
    double fragment_tolerance,
    bool fragment_unit_ppm,
    bool keep_only_deisotoped = false,
    Size min_isopeaks = 3,
    Size max_isopeaks = 10,
    bool make_single_charged = true)
  {
    if (in.empty()) return;

    MSSpectrum old_spectrum = in;

    // determine charge seeds and extend them
    vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
    vector<Int> features(old_spectrum.size(), -1);
    Int feature_number = 0;

    bool negative_mode = (max_charge < 0);
    Int step = negative_mode ? -1 : 1;

    for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
    {
      double current_mz = old_spectrum[current_peak].getPosition()[0];

      for (Int q = max_charge; abs(q) >= abs(min_charge); q -= step) // important: test charge hypothesis from high to low (in terms of absolute values)
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1) // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton) // test for missing peak
            {
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            }
            else
            {
/*
              // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
              Size n_extensions = extensions.size();
              if (n_extensions != 0)
              {
                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
                {
                  if (i < min_isopeaks)
                  {
                    has_min_isopeaks = false;
                  }
                  break;
                }
              }

*/
              // averagine check passed
              extensions.push_back(p);
            }
          }

          if (has_min_isopeaks)
          {
            //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
            mono_isotopic_peak[current_peak] = q;
            for (Size i = 0; i != extensions.size(); ++i)
            {
              features[extensions[i]] = feature_number;
            }
            feature_number++;
          }
        }
      }
    }

    in.clear(false);
    for (Size i = 0; i != old_spectrum.size(); ++i)
    {
      Int z = mono_isotopic_peak[i];
      if (keep_only_deisotoped)
      {
        if (z == 0)
        {
          continue;
        }

        // if already single charged or no decharging selected keep peak as it is
        if (!make_single_charged)
        {
          in.push_back(old_spectrum[i]);
        }
        else // make singly charged
        {
          Peak1D p = old_spectrum[i];
          if (negative_mode) // z < 0 in this case
          {
            z = abs(z);
            p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
          }
          else
          {
            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          }
          in.push_back(p);
        }
      }
      else
      {
        // keep all unassigned peaks
        if (features[i] < 0)
        {
          in.push_back(old_spectrum[i]);
          continue;
        }

        // convert mono-isotopic peak with charge assigned by deisotoping
        if (z != 0)
        {
          if (!make_single_charged)
          {
            in.push_back(old_spectrum[i]);
          }
          else // make singly charged
          {
            Peak1D p = old_spectrum[i];
            if (negative_mode) // z < 0 in this case
            {
              z = abs(z);
              p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
            }
            else
            {
              p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            }
            in.push_back(p);
          }
        }
      }
    }

    in.sortByPosition();
  }


  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra, bool negative_mode)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");

    // Note: we expect a higher number for NA than e.g., for peptides
    filter_param.setValue("peakcount", 50, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    // Note: we expect a higher number for NA than e.g., for peptides
    NLargest nlargest_filter = NLargest(1000);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      MSSpectrum& spec = exp[exp_index];

      // sort by mz
      spec.sortByPosition();

      if (spec.getPrecursors().empty()) continue; // this shouldn't happen
      Int precursor_charge = spec.getPrecursors()[0].getCharge();
      // deisotope
      Int coef = negative_mode ? -1 : 1;
      deisotopeAndSingleChargeMSSpectrum_(spec, coef, coef * precursor_charge, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 20, single_charge_spectra);

      // remove noise
      window_mower_filter.filterPeakSpectrum(spec);
      nlargest_filter.filterPeakSpectrum(spec);

      // sort (nlargest changes order)
      spec.sortByPosition();
    }
  }


  struct FragmentAnnotationDetail_
  {
    int charge;
    double mz;
    double intensity;

    bool operator<(const FragmentAnnotationDetail_& other) const
    {
      return std::tie(charge, mz, intensity) <
             std::tie(other.charge, other.mz, other.intensity);
    }

    bool operator==(const FragmentAnnotationDetail_& other) const
    {
      double mz_diff = fabs(mz - other.mz);
      double intensity_diff = fabs(intensity - other.intensity);
      // mz and intensity difference comparison actually not needed but kept for completeness
      return (mz_diff < 1e-6 && intensity_diff < 1e-6);
    }
  };


  void postProcessHits_(const PeakMap& exp,
                        vector<vector<AnnotatedHit>>& annotated_hits,
                        IdentificationData& id_data,
                        IdentificationData::InputFileKey file_key,
                        Size top_hits,
                        const vector<ConstRibonucleotidePtr>& fixed_modifications,
                        const vector<ConstRibonucleotidePtr>& variable_modifications,
                        Size max_variable_mods_per_oligo,
                        const vector<FASTAFile::FASTAEntry>& fasta_db,
                        const map<String, set<Size>>& oligo_map)
  {
    // remove all but top n scoring
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
    }

    IdentificationData::ScoreType score("hyperscore", true);
    IdentificationData::ScoreTypeKey score_key = id_data.registerScoreType(score).first;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spectrum = exp.getSpectrum(scan_index);
        Size charge = spectrum.getPrecursors()[0].getCharge();

        IdentificationData::DataQuery query(spectrum.getNativeID(), file_key, spectrum.getRT(), spectrum.getPrecursors()[0].getMZ());
        query.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        IdentificationData::DataQueryKey query_key = id_data.registerDataQuery(query).first;

        // create full oligo hit structure from annotated hits
        for (const AnnotatedHit& hit : annotated_hits[scan_index])
        {
          // get unmodified string
          LOG_DEBUG << "Hit sequence: " << hit.sequence << endl;
          NASequence seq = NASequence::fromString(hit.sequence);

          // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
          vector<NASequence> all_modified_oligos;
          ModifiedNASequenceGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), seq);
          ModifiedNASequenceGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), seq, max_variable_mods_per_oligo, all_modified_oligos);

          // reannotate NASequence object
          NASequence modified_seq = all_modified_oligos[hit.mod_index];

          IdentificationData::IdentifiedMoleculeKey oligo_key = id_data.registerOligo(modified_seq).first;

          IdentificationData::ScoreList scores;
          scores.push_back(make_pair(score_key, hit.score));
          IdentificationData::MoleculeQueryMatch match(charge, scores);
          match.peak_annotations = hit.annotations;

          id_data.addMoleculeQueryMatch(oligo_key, query_key, match);

          // add parent sequences:
          for (Size index : oligo_map.at(hit.sequence))
          {
            FASTAFile::FASTAEntry fasta = fasta_db[index];
            IdentificationData::ParentMetaData meta(
              IdentificationData::MT_RNA, fasta.sequence, fasta.description);
            IdentificationData::ParentMoleculeKey parent_key =
              id_data.registerParentMolecule(fasta.identifier, meta).first;

            id_data.addMoleculeParentMatch(oligo_key, parent_key);
          }
        }
      }
    }
  }


  IdentificationData::InputFileKey registerIDMetaData_(IdentificationData& id_data, const String& in_mzml, const vector<String>& primary_files,
                           const IdentificationData::DBSearchParameters& search_params)
  {
    IdentificationData::InputFileKey file_key = id_data.registerInputFile(in_mzml).first;
    IdentificationData::DataProcessingSoftware software(toolName_(), version_); // @TODO: add suitable processing action
    IdentificationData::ProcessingSoftwareKey software_key = id_data.registerDataProcessingSoftware(software).first;
    IdentificationData::SearchParamsKey search_key = id_data.registerDBSearchParameters(search_params).first;
    IdentificationData::DataProcessingStep step(software_key, vector<IdentificationData::InputFileKey>(1, file_key), primary_files);
    IdentificationData::ProcessingStepKey step_key = id_data.registerDataProcessingStep(step, search_key).first;
    // reference this step in all following ID data items, if applicable:
    id_data.setCurrentProcessingStep(step_key);
    return file_key;
  }


  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in_mzml = getStringOption_("in");
    String in_db = getStringOption_("database");
    String out = getStringOption_("out");
    String id_out = getStringOption_("id_out");
    String theo_ms2_out = getStringOption_("theo_ms2_out");
    String exp_ms2_out = getStringOption_("exp_ms2_out");

    IdentificationData::DBSearchParameters search_params;
    search_params.molecule_type = IdentificationData::MT_RNA;
    search_params.database = in_db;
    Int min_charge = getIntOption_("precursor:min_charge");
    Int max_charge = getIntOption_("precursor:max_charge");
    // @TODO: allow zero to mean "any charge state in the data"?
    if ((min_charge == 0) || (max_charge == 0))
    {
      LOG_ERROR << "Error: invalid charge state 0" << endl;
      return ILLEGAL_PARAMETERS;
    }
    // charges can be positive or negative, depending on data acquisition mode:
    if (((min_charge < 0) && (max_charge > 0)) ||
        ((min_charge > 0) && (max_charge < 0)))
    {
      LOG_ERROR << "Error: mixing positive and negative charges is not allowed"
                << endl;
      return ILLEGAL_PARAMETERS;
    }
    // min./max. are based on absolute value:
    if (abs(max_charge) < abs(min_charge)) swap(min_charge, max_charge);
    bool negative_mode = (max_charge < 0);
    Int step = negative_mode ? -1 : 1;
    for (Int charge = min_charge; abs(charge) <= abs(max_charge);
         charge += step)
    {
      search_params.charges.insert(charge);
    }
    search_params.precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    search_params.precursor_tolerance_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    search_params.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    search_params.fragment_tolerance_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");
    search_params.min_length = getIntOption_("oligo:min_size");

    StringList fixed_mod_names = getStringList_("modifications:fixed");
    search_params.fixed_mods.insert(fixed_mod_names.begin(), fixed_mod_names.end());

    StringList var_mod_names = getStringList_("modifications:variable");
    search_params.variable_mods.insert(var_mod_names.begin(), var_mod_names.end());

    vector<ConstRibonucleotidePtr> fixed_modifications = getModifications_(search_params.fixed_mods);
    vector<ConstRibonucleotidePtr> variable_modifications = getModifications_(search_params.variable_mods);

    // @TODO: add slots for these to "IdentificationData::DBSearchParameters"?
    IntList precursor_isotopes = getIntList_("precursor:isotopes");
    Size max_variable_mods_per_oligo = getIntOption_("modifications:variable_max_per_oligo");
    Int report_top_hits = getIntOption_("report:top_hits");

    // load MS2 map
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    progresslogger.startProgress(0, 1, "filtering spectra...");
    // @TODO: move this into the loop below (run only when checks pass)
    preprocessSpectra_(spectra, search_params.fragment_mass_tolerance, search_params.fragment_tolerance_ppm, true, negative_mode);
    progresslogger.endProgress();
    LOG_DEBUG << "preprocessed spectra: " << spectra.getNrSpectra() << endl;

    // build multimap of precursor mass to scan index
    multimap<double, Size> multimap_mass_2_scan_index;
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should be only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least one per nucleotide in the chain)
      if (precursor.size() == 1 && s_it->size() >= search_params.min_length)
      {
        int precursor_charge = precursor[0].getCharge();
        // the charge value in mzML seems to be always positive, so compare by
        // absolute value in negative mode:
        if ((negative_mode &&
             ((precursor_charge > abs(*search_params.charges.begin())) ||
              (precursor_charge < abs(*(--search_params.charges.end()))))) ||
            (!negative_mode &&
             ((precursor_charge < *search_params.charges.begin()) ||
              (precursor_charge > *(--search_params.charges.end())))))
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();

        // calculate precursor mass (optionally corrected for misassignment) and map it to MS scan index
        for (int isotope_number : precursor_isotopes)
        {
          double precursor_mass = precursor_mz * precursor_charge;
          if (negative_mode)
          {
            precursor_mass += Constants::PROTON_MASS_U * precursor_charge;
          }
          else
          {
            precursor_mass -= Constants::PROTON_MASS_U * precursor_charge;
          }

          // correct for monoisotopic misassignments of the precursor annotation
          if (isotope_number != 0) { precursor_mass -= isotope_number * Constants::C13C12_MASSDIFF_U; }

          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        }
      }
    }

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param(spectrum_generator.getParameters());
    // set nucleic acid-specific fragmentation pattern:
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "true");
    param.setValue("add_w_ions", "true");
    param.setValue("add_a-B_ions", "true");
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    spectrum_generator.setParameters(param);

    vector<vector<AnnotatedHit>> annotated_hits(spectra.size());
    MSExperiment exp_ms2_spectra, theo_ms2_spectra; // debug output

    progresslogger.startProgress(0, 1, "loading database from FASTA file...");
    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(in_db, fasta_db);
    progresslogger.endProgress();

    search_params.missed_cleavages = getIntOption_("oligo:missed_cleavages");
    String enzyme_name = getStringOption_("oligo:enzyme");
    search_params.digestion_enzyme = RNaseDB::getInstance()->getEnzyme(enzyme_name);
    RNaseDigestion digestor;
    digestor.setEnzyme(search_params.digestion_enzyme);
    digestor.setMissedCleavages(search_params.missed_cleavages);

    progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "scoring oligo models against spectra...");

    // lookup for processed oligos. must be defined outside of omp section and synchronized
    map<String, set<Size>> processed_oligos; // map: sequence -> FASTA index

    // set minimum size of oligo after digestion
    Size min_oligo_length = getIntOption_("oligo:min_size");

    Size hit_counter = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size fasta_index = 0; fasta_index < fasta_db.size(); ++fasta_index)
    {
      IF_MASTERTHREAD
      {
        progresslogger.setProgress((int)fasta_index * NUMBER_OF_THREADS);
      }

      vector<String> current_digest;
      digestor.digest(fasta_db[fasta_index].sequence, current_digest, min_oligo_length);

      for (vector<String>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_oligos_access)
#endif
        {
          auto pos = processed_oligos.find(*cit);
          if (pos != processed_oligos.end())
          {
            // oligo (and all modified variants) already processed so skip it
            pos->second.insert(fasta_index);
            already_processed = true;
          }
        }

        if (already_processed) continue;

#ifdef _OPENMP
#pragma omp critical (processed_oligos_access)
#endif
        {
          processed_oligos[*cit].insert(fasta_index);
        }

        vector<NASequence> all_modified_oligos;

        NASequence ns = NASequence::fromString(*cit);
        ModifiedNASequenceGenerator::applyFixedModifications(
          fixed_modifications.begin(), fixed_modifications.end(), ns);
        ModifiedNASequenceGenerator::applyVariableModifications(
          variable_modifications.begin(), variable_modifications.end(), ns,
          max_variable_mods_per_oligo, all_modified_oligos, true);

        for (SignedSize mod_idx = 0; mod_idx < (SignedSize)all_modified_oligos.size(); ++mod_idx)
        {
          const NASequence& candidate = all_modified_oligos[mod_idx];
          double candidate_mass = candidate.getMonoWeight();
          LOG_DEBUG << "candidate: " << candidate.toString() << " ("
                    << float(candidate_mass) << " Da)" << endl;

          // determine MS2 precursors that match to the current mass
          multimap<double, Size>::const_iterator low_it;
          multimap<double, Size>::const_iterator up_it;

          if (search_params.precursor_tolerance_ppm) // ppm
          {
            low_it = multimap_mass_2_scan_index.lower_bound(candidate_mass - 0.5 * candidate_mass * search_params.precursor_mass_tolerance * 1e-6);
            up_it = multimap_mass_2_scan_index.upper_bound(candidate_mass + 0.5 * candidate_mass * search_params.precursor_mass_tolerance * 1e-6);
          }
          else // Dalton
          {
            low_it = multimap_mass_2_scan_index.lower_bound(candidate_mass - 0.5 * search_params.precursor_mass_tolerance);
            up_it = multimap_mass_2_scan_index.upper_bound(candidate_mass + 0.5 * search_params.precursor_mass_tolerance);
          }

          if (low_it == up_it) continue; // no matching precursor in data

          // create theoretical spectrum
          PeakSpectrum theo_spectrum;

          // add peaks for b and y ions with charge 1
          Int charge = negative_mode ? -1 : 1;
          spectrum_generator.getSpectrum(theo_spectrum, candidate, charge, charge);

          // sort by mz
          theo_spectrum.sortByPosition();

          for (; low_it != up_it; ++low_it)
          {
            LOG_DEBUG << "matching precursor mass: " << float(low_it->first)
                      << endl;

            const Size& scan_index = low_it->second;
            const PeakSpectrum& exp_spectrum = spectra[scan_index];

            vector<PeptideHit::PeakAnnotation> annotations;
            double score = MetaboliteSpectralMatching::computeHyperScore(
              search_params.fragment_mass_tolerance,
              search_params.fragment_tolerance_ppm,
              exp_spectrum, theo_spectrum, annotations);

            if (!exp_ms2_out.empty())
            {
              exp_ms2_spectra.addSpectrum(exp_spectrum);
            }
            if (!theo_ms2_out.empty())
            {
              theo_spectrum.setName(candidate.toString());
              theo_ms2_spectra.addSpectrum(theo_spectrum);
            }

            // no hit
            if (score < 1e-16)
            {
              continue;
            }

            // add oligo hit
            AnnotatedHit ah;
            ah.sequence = *cit;
            ah.mod_index = mod_idx;
            ah.score = score;
            ah.annotations = annotations;

#ifdef _OPENMP
#pragma omp atomic
#endif
            ++hit_counter;

#ifdef DEBUG_NASEARCH
            cout << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
#pragma omp critical (annotated_hits_access)
#endif
            {
              annotated_hits[scan_index].push_back(ah);
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    LOG_INFO << "Undigested nucleic acids: " << fasta_db.size() << endl;
    LOG_INFO << "Oligonucleotides: " << processed_oligos.size() << endl;
    LOG_INFO << "Search hits: " << hit_counter << endl;

    if (!exp_ms2_out.empty())
    {
      MzMLFile().store(exp_ms2_out, exp_ms2_spectra);
    }
    if (!theo_ms2_out.empty())
    {
      MzMLFile().store(theo_ms2_out, theo_ms2_spectra);
    }

    IdentificationData id_data;
    vector<String> primary_files;
    spectra.getPrimaryMSRunPath(primary_files);
    IdentificationData::InputFileKey file_key = registerIDMetaData_(id_data, in_mzml, primary_files, search_params);

    progresslogger.startProgress(0, 1, "post-processing search hits...");
    postProcessHits_(spectra, annotated_hits, id_data, file_key,
                     report_top_hits, fixed_modifications,
                     variable_modifications, max_variable_mods_per_oligo,
                     fasta_db, processed_oligos);
    progresslogger.endProgress();

    // store results
    MzTab results = id_data.exportMzTab();
    LOG_DEBUG << "Nucleic acid rows: "
              << results.getNucleicAcidSectionRows().size()
              << "\nOligonucleotide rows: "
              << results.getOligonucleotideSectionRows().size()
              << "\nOligo-spectrum match rows: "
              << results.getOSMSectionRows().size() << endl;

    MzTabFile().store(out, results);

    // dummy "peptide" results:
    if (!id_out.empty())
    {
      vector<PeptideIdentification> peptides;
      vector<ProteinIdentification> proteins(1);
      proteins[0].setIdentifier("id");
      proteins[0].setDateTime(DateTime::now());
      proteins[0].setSearchEngine(toolName_());
      map<IdentificationData::DataQueryKey, PeptideIdentification> id_map;
      for (const auto& osm : id_data.query_matches)
      {
        IdentificationData::DataQueryKey query_key = osm.first.second;
        IdentificationData::IdentifiedMoleculeKey oligo_key = osm.first.first;
        const IdentificationData::MoleculeQueryMatch& match = osm.second;
        const NASequence& seq = id_data.identified_oligos.left.at(oligo_key);
        PeptideHit hit;
        hit.setMetaValue("label", seq.toString());
        hit.setScore(match.scores.back().second);
        hit.setCharge(match.charge);
        hit.setPeakAnnotations(match.peak_annotations);
        id_map[query_key].insertHit(hit);
      }
      // there should be only one score type:
      const IdentificationData::ScoreType& score_type =
        id_data.score_types.left.begin()->second;
      for (auto& id_pair : id_map)
      {
        const IdentificationData::DataQuery& query =
          id_data.data_queries.left.at(id_pair.first);
        PeptideIdentification& peptide = id_pair.second;
        peptide.setRT(query.rt);
        peptide.setMZ(query.mz);
        peptide.setMetaValue("spectrum_reference", query.data_id);
        peptide.setScoreType(score_type.name);
        peptide.setHigherScoreBetter(score_type.higher_better);
        peptide.setIdentifier("id");
        peptides.push_back(peptide);
      }
      IdXMLFile().store(id_out, proteins, peptides);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  NucleicAcidSearchEngine tool;
  return tool.main(argc, argv);
}
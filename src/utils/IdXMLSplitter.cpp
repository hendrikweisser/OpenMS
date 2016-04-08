// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IdXMLSplitter IdXMLSplitter

    @brief Splits an idXML file into parts corresponding to multiple raw data files

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IdXMLSplitter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IdXMLSplitter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIdXMLSplitter :
  public TOPPBase
{
public:

  TOPPIdXMLSplitter() :
    TOPPBase("IdXMLSplitter", "Splits an idXML file into parts corresponding to multiple raw data files", false)
  {
  }

protected:

  struct MS2Info
  {
    double rt, mz;
    Size origin;
  };


  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file: ID data to split");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerInputFileList_("in_mz", "<files>", vector<String>(), "Input files: raw data");
    setValidFormats_("in_mz", ListUtils::create<String>("mzML"));
    registerOutputFileList_("out", "<files>", vector<String>(), "Output files: split ID data (or an existing directory to create those files in)", false);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
  }


  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    vector<String> in_mz = getStringList_("in_mz"), out = getStringList_("out");

    if (in_mz.size() < 2)
    {
      writeLog_("Error: Need at least two files passed via parameter 'in_mz'.");
      return ILLEGAL_PARAMETERS;
    }
    if (in_mz.size() != out.size())
    {
      if ((out.size() == 1) && (File::isDirectory(out[0])))
      {
        // generate output filenames based on "in_mz":
        String dir = out[0];
        out.resize(in_mz.size());
        writeLog_("Generating " + String(out.size()) + " output filenames...");
        map<String, Size> basenames;
        for (Size i = 0; i < in_mz.size(); ++i)
        {
          String basename = File::basename(in_mz[i]);
          if (basename.hasSuffix(".gz"))
          {
            basename = File::removeExtension(basename);
          }
          basename = File::removeExtension(basename);
          out[i] = dir + String(QDir::separator()) + basename;
          Size& count = ++basenames[basename]; // new basename: count = 1
          if (count > 1)
          {
            out[i] += "_" + String(count);
            ++count;
          }
          out[i] += ".idXML";
        }
      }
      else
      {
        writeLog_("Error: Same number of files must be passed via parameters 'in_mz' and 'out', or 'out' must be an existing directory.");
        return ILLEGAL_PARAMETERS;
      }
    }

    vector<PeptideIdentification> peptides;
    vector<ProteinIdentification> proteins;
    IdXMLFile().load(in, proteins, peptides);

    if (peptides.empty())
    {
      writeLog_("Error: No peptide IDs found.");
      return INPUT_FILE_EMPTY;
    }

    writeLog_("Collecting MS2 locations (RT, m/z)...");
    double rt_min = numeric_limits<double>::max(), rt_max = 0,
      mz_min = numeric_limits<double>::max(), mz_max = 0;
    vector<MS2Info> ms2_infos;
    for (Size i = 0; i < in_mz.size(); ++i)
    {
      MSExperiment<> experiment;
      MzMLFile().load(in_mz[i], experiment);
      for (MSExperiment<>::Iterator it = experiment.begin();
           it != experiment.end(); ++it)
      {
        if (it->getMSLevel() == 2)
        {
          double rt = it->getRT(), mz = it->getPrecursors()[0].getMZ();
          if (rt < rt_min) rt_min = rt;
          if (rt > rt_max) rt_max = rt;
          if (mz < mz_min) mz_min = mz;
          if (mz > mz_max) mz_max = mz;
          MS2Info info = {rt, mz, i};
          ms2_infos.push_back(info);
        }
      }
    }

    writeLog_("Matching " + String(peptides.size()) + " peptide IDs to " +
              String(ms2_infos.size()) + " MS2 spectra...");

    double rt_range = rt_max - rt_min, mz_range = mz_max - mz_min;
    // index into "ms2_infos" for each peptide ID:
    vector<Size> nearest_ms2(peptides.size());
    set<Size> matches;
    ProgressLogger progress;
    progress.setLogType(log_type_);
    progress.startProgress(0, peptides.size(), "matching peptide IDs");
    Size conflicts = 0;
    for (Size i = 0; i < peptides.size(); ++i)
    {
      Size best_index = 0;
      double best_dist = 2.0;
      for (Size j = 0; j < ms2_infos.size(); ++j)
      {
        double rt_dist = abs(peptides[i].getRT() - ms2_infos[j].rt) / rt_range,
          mz_dist = abs(peptides[i].getMZ() - ms2_infos[j].mz) / mz_range;
        double dist = rt_dist * rt_dist + mz_dist * mz_dist;
        if (dist < best_dist)
        {
          best_dist = dist;
          best_index = j;
        }
      }
      nearest_ms2[i] = best_index;
      if (!matches.insert(best_index).second) // MS2 spectrum "claimed" before
      {
        ++conflicts;
        String rt = ms2_infos[best_index].rt, mz = ms2_infos[best_index].mz,
          origin = ms2_infos[best_index].origin;
        writeLog_("Conflict for MS2 spectrum: RT = " + rt + ", m/z = " + mz +
                  ", origin = " + origin);
      }
      progress.setProgress(i + 1);
    }
    progress.endProgress();

    vector<vector<PeptideIdentification> > peptides_by_origin(in_mz.size());
    for (Size i = 0; i < nearest_ms2.size(); ++i)
    {
      Size origin = ms2_infos[nearest_ms2[i]].origin;
      peptides_by_origin[origin].push_back(peptides[i]);
    }

    for (Size i = 0; i < in_mz.size(); ++i)
    {
      writeLog_(String(peptides_by_origin[i].size()) + " peptide IDs in file " +
                       String(i + 1));
      vector<ProteinIdentification> proteins_copy = proteins;
      IDFilter::removeUnreferencedProteins(proteins_copy,
                                           peptides_by_origin[i]);
      IdXMLFile().store(out[i], proteins_copy, peptides_by_origin[i]);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIdXMLSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond

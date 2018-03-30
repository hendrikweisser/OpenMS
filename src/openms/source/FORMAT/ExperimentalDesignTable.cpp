// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/ExperimentalDesignTable.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <QtCore/QString>
#include <QtCore/QFileInfo>

#include <iostream>

using namespace std;

namespace OpenMS
{
    using RunRows = std::vector<ExperimentalDesignTable::RunRow>;

    // Parse Error of filename if test holds
    void parseErrorIf(const bool test, const String &filename, const String &message)
    {
      if (test)
      {
        throw Exception::ParseError(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                filename,
                "Error: " + message);
      }
    }

    void ExperimentalDesignTable::parseHeader_(
                      const StringList &header,
                      const String &filename,
                      std::map <String, Size> &column_map,
                      const std::set <String> &required,
                      const std::set <String> &optional,
                      const bool allow_other_header)
    {
      // Headers as set
      std::set <String> header_set(header.begin(), header.end());
      parseErrorIf(header_set.size() != header.size(), filename,
                   "Some column headers of the table appear multiple times!");

      // Check that all required headers are there
      for (const String &req_header : required)
      {
        parseErrorIf( ! ListUtils::contains(header, req_header), filename,
                     "Missing column header: " + req_header);
      }
      // Assign index in column map and check for weird headers
      for (Size i = 0; i < header.size(); ++i)
      {
        const String &h = header[i];

        // A header is unexpected if it is neither required nor optional and we do not allow other headers
        const bool header_unexpected = (required.find(h) == required.end()) && (optional.find(h) == optional.end());
        parseErrorIf(
                !allow_other_header && header_unexpected,
                filename,
                "Header not allowed in this section of the Experimental Design: " + h
        );
        column_map[h] = i;
      }
    }

    String findSpectraFile(const String &spec_file, const String &tsv_file, const bool require_spectra_file)
    {
      String result;
      QFileInfo spectra_file_info(spec_file.toQString());
      if (spectra_file_info.isRelative())
      {
        // file name is relative so we need to figure out the correct folder

        // first check folder relative to folder of design file
        // to allow, for example, a design in ./design.tsv and spectra in ./spectra/a.mzML
        // where ./ is the same folder
        QFileInfo design_file_info(tsv_file.toQString());
        QString design_file_relative(design_file_info.absolutePath());
        design_file_relative = design_file_relative + "/" + spec_file.toQString();

        if (File::exists(design_file_relative))
        {
          result = design_file_relative.toStdString();
        }
        else
        {
          // check current folder
          String f = File::absolutePath(spec_file);
          if (File::exists(f))
          {
            result = f;
          }
        }

        // if result still empty, just use the provided value
        if (result.empty())
        {
          result = spec_file;
        }
      }
      else
      {
        // set to absolute path
        result = spec_file;
      }

      // Fail if the existence of the spectra files is required but they do not exist
      if (require_spectra_file && !File::exists(result))
      {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
                "Error: Spectra file does not exist: '" + result + "'");
      }
      return result;
    }


    // static
    ExperimentalDesignTable ExperimentalDesignTable::load(const String &tsv_file, const bool require_spectra_file)
    {
      ExperimentalDesignTable design;
      design.run_section_.clear();

      bool has_sample(false);
      bool has_channel(false);

      // Maps the column header string to the column index for
      // the run section
      std::map <String, Size> run_column_header_to_index;

      unsigned line_number(0);

      enum ParseState { RUN_HEADER, RUN_CONTENT, SAMPLE_HEADER, SAMPLE_CONTENT };

      ParseState state(RUN_HEADER);
      Size n_col = 0;

      const TextFile text_file(tsv_file, true);
      for (String s : text_file)
      {
        // skip empty lines (except in state RUN_CONTENT, where the sample table is read)
        const String line(s.trim());

        if (line.empty() && state != RUN_CONTENT)
        {
          continue;
        }

        // Now split the line into individual cells
        StringList cells;
        line.split("\t", cells);

        // Trim whitespace from all cells (so , foo , and  ,foo, is the same)
        std::transform(cells.begin(), cells.end(), cells.begin(),
                       [](String cell) -> String { return cell.trim(); });

        if (state == RUN_HEADER)
        {
          state = RUN_CONTENT;
          parseHeader_(
                  cells,
                  tsv_file,
                  run_column_header_to_index,
                  {"Run", "Fraction", "Path(Spectra File)"},
                  {"Channel", "Sample"}, false
          );
          has_channel = run_column_header_to_index.find("Channel") != run_column_header_to_index.end();
          has_sample = run_column_header_to_index.find("Sample") != run_column_header_to_index.end();
          n_col = run_column_header_to_index.size();
        }
          // End of run section lines, empty line separates run and sample table
        else if (state == RUN_CONTENT && line.empty())
        {
          // Next line is header of Sample table
          state = SAMPLE_HEADER;
        }
          // Line is run line of run section
        else if (state == RUN_CONTENT)
        {
          parseErrorIf(n_col != cells.size(), tsv_file,
                       "Wrong number of records in line");

          RunRow row;

          // Assign run and fraction
          row.run = cells[run_column_header_to_index["Run"]].toInt();
          row.fraction = cells[run_column_header_to_index["Fraction"]].toInt();

          // Assign channel
          row.channel = has_channel ? cells[run_column_header_to_index["Channel"]].toInt() : 1;

          // Assign sample number
          if (has_sample)
          {
            row.sample = cells[run_column_header_to_index["Sample"]].toInt();
          }
          else
          {
            row.sample = has_channel ? row.channel : row.run;
          }

          // Spectra files
          row.path = findSpectraFile(
                  cells[run_column_header_to_index["Path(Spectra File)"]],
                  tsv_file,
                  require_spectra_file);
          design.run_section_.push_back(row);
        }
          // Parse header of the Condition Table
        else if (state == SAMPLE_HEADER)
        {
          state = SAMPLE_CONTENT;
          line_number = 0;
          parseHeader_(
                  cells,
                  tsv_file,
                  design.sample_section_.columnname_to_columnindex_,
                  {"Sample"}, {}, true
          );
          n_col = design.sample_section_.columnname_to_columnindex_.size();
        }
          // Parse Sample Row
        else if (state == SAMPLE_CONTENT)
        {
          // Parse Error if sample appears multiple times
          unsigned sample = cells[design.sample_section_.columnname_to_columnindex_["Sample"]].toInt();
          parseErrorIf(
                  design.sample_section_.sample_to_rowindex_.find(sample) !=
                  design.sample_section_.sample_to_rowindex_.end(),
                  tsv_file,
                  "Sample: " + String(sample) + " appears multiple times in the sample table"
          );
          design.sample_section_.sample_to_rowindex_[sample] = line_number++;
          design.sample_section_.content_.push_back(cells);
        }
      }

      design.sort_();
      design.checkValidRunSection_();
      return design;
    }

    ExperimentalDesignTable ExperimentalDesignTable::fromConsensusMap(const ConsensusMap &cm)
    {
      ExperimentalDesignTable experimental_design;
      // path of the original MS run (mzML / raw file)
      StringList ms_run_paths;
      cm.getPrimaryMSRunPath(ms_run_paths);

      // no fractionation -> as many runs as samples
      // each consensus element corresponds to one sample abundance
      size_t sample(1);
      ExperimentalDesignTable::RunRows rows;
      for (const auto &f : ms_run_paths)
      {
        ExperimentalDesignTable::RunRow r;
        r.path = f;
        r.fraction = 1;
        r.sample = sample;
        r.run = sample;
        r.channel = 1; // TODO MULTIPLEXING: adapt for non-label-free
        rows.push_back(r);
        ++sample;
      }
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (ConsensusMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesignTable ExperimentalDesignTable::fromFeatureMap(const FeatureMap &fm)
    {
      ExperimentalDesignTable experimental_design;
      // path of the original MS run (mzML / raw file)
      StringList ms_paths;
      fm.getPrimaryMSRunPath(ms_paths);

      if (ms_paths.size() != 1)
      {
        throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "FeatureMap annotated with " + String(ms_paths.size()) + " MS files. Must be exactly one.");
      }

      // Feature map is simple. One file, one fraction, one sample, one run
      ExperimentalDesignTable::RunRow r;
      r.path = ms_paths[0];
      r.fraction = 1;
      r.sample = 1;
      r.run = 1;
      r.channel = 1;

      ExperimentalDesignTable::RunRows rows(1, r);
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (FeatureMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesignTable ExperimentalDesignTable::fromIdentifications(const vector <ProteinIdentification> &proteins)
    {
      ExperimentalDesignTable experimental_design;
      // path of the original MS files (mzML / raw file)
      StringList ms_run_paths;
      for (const auto &protein : proteins)
      {
        StringList tmp_ms_run_paths;
        protein.getPrimaryMSRunPath(tmp_ms_run_paths);
        if (tmp_ms_run_paths.size() != 1)
        {
          throw Exception::MissingInformation(
                  __FILE__,
                  __LINE__,
                  OPENMS_PRETTY_FUNCTION,
                  "ProteinIdentification annotated with " + String(tmp_ms_run_paths.size()) +
                  " MS files. Must be exactly one.");
        }
        ms_run_paths.push_back(tmp_ms_run_paths[0]);
      }

      // no fractionation -> as many runs as samples
      // each identification run corresponds to one sample abundance
      unsigned sample(1);
      ExperimentalDesignTable::RunRows rows;
      for (const auto &f : ms_run_paths)
      {
        ExperimentalDesignTable::RunRow r;
        r.path = f;
        r.fraction = 1;
        r.sample = sample;
        r.run = sample;
        r.channel = 1;

        rows.push_back(r);
        ++sample;
      }
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (Identification derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    map<unsigned, vector<String> > ExperimentalDesignTable::getFractionToMSFilesMapping() const
    {
      map<unsigned, vector<String> > ret;

      for (RunRow const &r : run_section_)
      {
        ret[r.fraction].emplace_back(r.path);
      }
      return ret;
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesignTable::pathChannelMapper(
            const bool basename,
            unsigned (*f)(const ExperimentalDesignTable::RunRow &row)) const
    {
      map<pair<String, unsigned>, unsigned> ret;
      for (RunRow const &r : run_section_)
      {
        const String path = String(r.path);
        pair<String, unsigned> tpl = make_pair((basename ? File::basename(path) : path), r.channel);
        ret[tpl] = f(r);
      }
      return ret;
    }


    map<pair<String, unsigned>, unsigned> ExperimentalDesignTable::getPathChannelToSampleMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.sample; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesignTable::getPathChannelToFractionMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.fraction; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesignTable::getPathChannelToRunMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.run; });
    }


    bool ExperimentalDesignTable::sameNrOfMSFilesPerFraction() const
    {
      map<unsigned, vector<String>> frac2files = getFractionToMSFilesMapping();
      if (frac2files.size() <= 1) { return true; }

      Size files_per_fraction(0);
      for (auto const &f : frac2files)
      {
        if (files_per_fraction == 0) // first fraction, initialize
        {
          files_per_fraction = f.second.size();
        }
        else // fraction >= 2
        {
          // different number of associated MS files?
          if (f.second.size() != files_per_fraction)
          {
            return false;
          }
        }
      }
      return true;
    }

    /* Implementation for the Experimental Design and RunSection */

    const RunRows& ExperimentalDesignTable::getRunSection() const
    {
      return run_section_;
    }


    void ExperimentalDesignTable::setRunSection(const RunRows& run_section)
    {
      run_section_ = run_section;
      sort_();
      checkValidRunSection_();
    }

    unsigned ExperimentalDesignTable::getNumberOfSamples() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                return f1.sample < f2.sample;
                              })->sample;
    }

    unsigned ExperimentalDesignTable::getNumberOfFractions() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                  return f1.fraction < f2.fraction;
                              })->fraction;
    }

    // @return the number of channels per file
    unsigned ExperimentalDesignTable::getNumberOfChannels() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                return f1.fraction < f2.fraction;
                              })->channel;
    }

    // @return the number of MS files (= fractions * runs)
    unsigned ExperimentalDesignTable::getNumberOfMSFiles() const
    {
      std::set<std::string> unique_paths;
      for (auto const & r : run_section_) { unique_paths.insert(r.path); }
      return unique_paths.size();
    }

    // @return the number of runs (before fractionation)
    // Allows to group fraction ids and source files
    unsigned ExperimentalDesignTable::getNumberOfPrefractionationRuns() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                  return f1.run < f2.run;
                              })->run;
    }

    // @return sample index (depends on run and channel)
    unsigned ExperimentalDesignTable::getSample(unsigned run, unsigned channel)
    {
      return std::find_if(run_section_.begin(), run_section_.end(),
                          [&run, &channel](const RunRow& r)
                          {
                              return r.run == run && r.channel == channel;
                          })->sample;
    }

    const ExperimentalDesignTable::SampleSection& ExperimentalDesignTable::getSampleSection() const
    {
      return sample_section_;
    }

    std::vector< String > ExperimentalDesignTable::getFileNames(const bool basename) const
    {
      std::vector<String> filenames;
      for (const RunRow &row : run_section_)
      {
        const String path = String(row.path);
        filenames.push_back(basename ? path : File::basename(path));
      }
      return filenames;
    }

    template<typename T>
    void ExperimentalDesignTable::errorIfAlreadyExists(std::set<T> &container, T &item, const String &message)
    {
      if (container.find(item) != container.end())
      {
       throw Exception::MissingInformation(
       __FILE__,
       __LINE__,
        OPENMS_PRETTY_FUNCTION, message);
      }
      container.insert(item);
    }

    void ExperimentalDesignTable::checkValidRunSection_()
    {
      if (getNumberOfMSFiles() == 0)
      {
        throw Exception::MissingInformation(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "No MS files provided.");
      }

    std::set< std::tuple< unsigned, unsigned, unsigned > > run_fraction_sample_set;
    std::set< std::tuple< unsigned, unsigned, unsigned > > run_fraction_channel_set;
    std::set< std::tuple< std::string, unsigned > > path_channel_set;
    std::map< std::tuple< unsigned, unsigned >, std::set< unsigned > > run_channel_to_sample;

    for (const RunRow &row : run_section_)
    {
      // Fail if sample section does not contain the run
      if (sample_section_.hasSample(row.sample) == false)
      {
        throw Exception::MissingInformation(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Sample Section does not contain sample for run " + String(row.run));
      }

      // RUN_FRACTION_SAMPLE TUPLE
      std::tuple<unsigned, unsigned, unsigned> run_fraction_sample = std::make_tuple(row.run, row.fraction, row.sample);
      errorIfAlreadyExists(
        run_fraction_sample_set,
        run_fraction_sample,
      "(Run, Fraction, Sample) combination can only appear once");

      // RUN_FRACTION_CHANNEL TUPLE
      std::tuple<unsigned, unsigned, unsigned> run_fraction_channel = std::make_tuple(row.run, row.fraction, row.channel);
      errorIfAlreadyExists(
        run_fraction_channel_set,
        run_fraction_channel,
      "(Run, Fraction, Channel) combination can only appear once");


      // PATH_CHANNEL_TUPLE
      std::tuple<std::string, unsigned> path_channel = std::make_tuple(row.path, row.channel);
      errorIfAlreadyExists(
        path_channel_set,
        path_channel,
        "(Path, Channel) combination can only appear once");

      // RUN_CHANNEL TUPLE
      std::tuple<unsigned, unsigned> run_channel = std::make_tuple(row.run, row.channel);
      run_channel_to_sample[run_channel].insert(row.sample);

      if (run_channel_to_sample[run_channel].size() > 1)
      {
        throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Multiple Samples encountered for the same Run and the same Channel");
      }
    }
  }

  std::vector<unsigned> ExperimentalDesignTable::getChannels() const
  {
    std::vector<unsigned> channels;
    for (const RunRow &row : run_section_)
    {
      channels.push_back(row.channel);
    }
    return channels;
  }

  std::vector<unsigned> ExperimentalDesignTable::getFractions() const
  {
    std::vector<unsigned> fractions;
    for (const RunRow &row : run_section_)
    {
      fractions.push_back(row.fraction);
    }
    return fractions;
  }

  void ExperimentalDesignTable::sort_()
  {
    std::sort(run_section_.begin(), run_section_.end(),
                [](const RunRow& a, const RunRow& b)
                {
                  return std::tie(a.run, a.fraction, a.channel, a.sample, a.path) <
                         std::tie(b.run, b.fraction, b.channel, b.sample, b.path);
                });
  }

  /* Implementations of SampleSection */

  std::set<unsigned> ExperimentalDesignTable::SampleSection::getSamples() const
  {
    std::set<unsigned> samples;
    for (const auto &kv : sample_to_rowindex_)
    {
      samples.insert(kv.first);
    }
    return samples;
  }

  std::set< String > ExperimentalDesignTable::SampleSection::getFactors() const
  {
    std::set<String> factors;
    for (const auto &kv : columnname_to_columnindex_)
    {
      factors.insert(kv.first);
    }
    return factors;
  }

  bool ExperimentalDesignTable::SampleSection::hasSample(const unsigned sample) const
  {
    return sample_to_rowindex_.find(sample) != sample_to_rowindex_.end();
  }

  bool ExperimentalDesignTable::SampleSection::hasFactor(const String &factor) const
  {
    return columnname_to_columnindex_.find(factor) != columnname_to_columnindex_.end();
  }

  String ExperimentalDesignTable::SampleSection::getFactorValue(const unsigned sample, const String &factor)
  {
   if (! hasSample(sample))
   {
    throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Sample " + String(sample) + " is not present in the Experimental Design");
   }
   if (! hasFactor(factor))
   {
    throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Factor " + factor + " is not present in the Experimental Design");
   }
   StringList sample_row = content_[sample_to_rowindex_[sample]];
   const Size col_index = columnname_to_columnindex_[factor];
   return sample_row[col_index];
  }
}

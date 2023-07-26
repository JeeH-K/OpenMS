//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <chrono>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>


using namespace std;


namespace OpenMS
{
  void FLASHDeconvAlgorithm::filterLowPeaks(MSExperiment& map, Size count)
  {
    for (auto& it : map)
    {
      double threshold;
      if (it.getType(false) == SpectrumSettings::CENTROID)
      {
        if (it.size() <= count)
        {
          return;
        }
        it.sortByIntensity(true);
        threshold = it[count].getIntensity();
      }
      else
      {
        if (it.size() <= count)
        {
          continue;
        }

        it.sortByIntensity(true);
        double max_intensity = log10(it[0].getIntensity());
        double min_intensity = 0;
        for (auto& p : it)
        {
          if (p.getIntensity() <= 0)
          {
            break;
          }
          min_intensity = log10(p.getIntensity());
        }
        Size bin_size = 500;
        std::vector<int> freq(bin_size + 1, 0);
        for (auto& p : it)
        {
          if (p.getIntensity() <= 0)
          {
            break;
          }
          Size bin = round((log10(p.getIntensity()) - min_intensity) / (max_intensity - min_intensity) * bin_size);
          freq[bin]++;
        }

        int mod_bin = std::distance(freq.begin(), std::max_element(freq.begin(), freq.end())); // most frequent intensity is the threshold to distinguish between signal and noise

        threshold =
          3.0 * (pow(10.0, (double)mod_bin / bin_size * (max_intensity - min_intensity) +
                             min_intensity)); // multiply by 3 to the most frequent intensity to make sure more signal component remains. Later this could be determined to use signal-to-noise ratio.
      }
      // pop back the low intensity peaks using threshold
      while (it.size() > 0 && it[it.size() - 1].getIntensity() <= threshold)
      {
        it.pop_back();
      }

      it.sortByPosition();
    }
  }



   int FLASHDeconvAlgorithm::SpectrumDeconvolution(MSExperiment& map, ProgressLogger& progresslogger, uint& current_max_ms_level, const int target_precursor_charge, std::vector<size_t>& spec_cntr,
                                                   std::unordered_map<UInt, std::vector<DeconvolvedSpectrum>>& last_deconvolved_spectra, SpectralDeconvolution fd,
                                                   std::map<int, std::vector<std::vector<float>>> precursor_map_for_real_time_acquisition, double target_precursor_mass,
                                                   std::map<int, PeakGroup>& precursor_peak_groups, double topFD_SNR_threshold, double& expected_identification_count, String out_mzml_file,
                                                   int mzml_charge, uint current_min_ms_level, const DoubleList tols, MSExperiment& exp, String out_anno_mzml_file, MSExperiment& exp_annotated,
                                                   int num_last_deconvolved_spectra, const int merge, std::map<int, double>& scan_rt_map, const bool report_dummy,
                                                   SpectralDeconvolution fd_charge_dummy, SpectralDeconvolution fd_noise_dummy, SpectralDeconvolution fd_iso_dummy,
                                                   std::vector<DeconvolvedSpectrum>& dummy_deconvolved_spectra, std::vector<size_t>& qspec_cntr, std::vector<size_t>& mass_cntr,
                                                   std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<double>& elapsed_deconv_cpu_secs, 
                                                   std::vector<double>& elapsed_deconv_wall_secs)


  {



    for (auto it = map.begin(); it != map.end(); ++it)
    {
      int scan_number = map.getSourceFiles().empty() ? -1 : SpectrumLookup::extractScanNumber(it->getNativeID(), map.getSourceFiles()[0].getNativeIDTypeAccession());

      if (scan_number < 0)
      {
        scan_number = (int)std::distance(map.begin(), it) + 1;
      }
       
      progresslogger.nextProgress();

      if (it->empty())
      {
        continue;
      }

      uint ms_level = it->getMSLevel();
      if (ms_level > current_max_ms_level)
      {
        continue;
      }

      if (ms_level > 1 && target_precursor_charge != 0 && it->getPrecursors().size() == 0)
      {
        OPENMS_LOG_INFO << "Target precursor charge is set but no precursor m/z is found in MS2 spectra. Specify target precursor m/z with -target_precursor_mz option" << std::endl;
        return 9;
      }

      spec_cntr[ms_level - 1]++;
      auto deconv_begin = clock();
      auto deconv_t_start = std::chrono::high_resolution_clock::now();

      // for MS>1 spectrum, register precursor
      std::vector<DeconvolvedSpectrum> precursor_specs;

      if (ms_level > 1 && last_deconvolved_spectra.find(ms_level - 1) != last_deconvolved_spectra.end())
      {
        precursor_specs = (last_deconvolved_spectra[ms_level - 1]);
      }

      fd.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_real_time_acquisition);
      auto& deconvolved_spectrum = fd.getDeconvolvedSpectrum();
      if (deconvolved_spectrum.empty())
      {
        continue;
      }

      if (ms_level > 1 && target_precursor_charge != 0)
      {
        auto precursor = it->getPrecursors()[0];
        target_precursor_mass = (precursor.getMZ() - FLASHDeconvHelperStructs::getChargeMass(target_precursor_charge > 0)) * std::abs(target_precursor_charge);
        // precursor.setCharge(target_precursor_charge);
        PeakGroup precursorPeakGroup(1, std::abs(target_precursor_charge), target_precursor_charge > 0);
        precursorPeakGroup.push_back(FLASHDeconvHelperStructs::LogMzPeak());
        precursorPeakGroup.setMonoisotopicMass(target_precursor_mass);
        precursorPeakGroup.setSNR(1.0);

        precursorPeakGroup.setChargeSNR(std::abs(target_precursor_charge), 1.0);
        precursorPeakGroup.Qscore(1.0);
        deconvolved_spectrum.setPrecursor(precursor);
        deconvolved_spectrum.setPrecursorPeakGroup(precursorPeakGroup);
      }

      if (it->getMSLevel() > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
      {
        precursor_peak_groups[scan_number] = deconvolved_spectrum.getPrecursorPeakGroup();
        if (deconvolved_spectrum.getPrecursorPeakGroup().getChargeSNR(std::abs(deconvolved_spectrum.getPrecursorCharge())) >= topFD_SNR_threshold)
        {
          expected_identification_count += deconvolved_spectrum.getPrecursorPeakGroup().getQscore();
        }
      }
      bool deconved_mzML_written = false;
      if (!out_mzml_file.empty())
      {
        if (!deconvolved_spectrum.empty())
        {
          auto dspec = deconvolved_spectrum.toSpectrum(mzml_charge, current_min_ms_level, tols[ms_level - 1], false);
          if (dspec.size() > 0)
          {
            exp.addSpectrum(dspec);
            deconved_mzML_written = true;
          }
        }
      }

      if (!out_anno_mzml_file.empty())
      {
        if (out_mzml_file.empty() || deconved_mzML_written)
        {
          auto anno_spec = MSSpectrum(*it);

          if (!deconvolved_spectrum.empty())
          {
            std::stringstream val {};

            for (auto& pg : deconvolved_spectrum)
            {
              val << std::to_string(pg.getMonoMass()) << ":";
              for (size_t k = 0; k < pg.size(); k++)
              {
                auto& p = pg[k];
                auto pindex = anno_spec.findNearest(p.mz);
                val << pindex;
                if (k < pg.size() - 1)
                {
                  val << ",";
                }
              }
              val << ";";
            }
            anno_spec.setMetaValue("DeconvMassPeakIndices", val.str());
            exp_annotated.addSpectrum(anno_spec);
          }
        }
      }
      if (ms_level < current_max_ms_level)
      {
        if ((int)last_deconvolved_spectra[ms_level].size() >= num_last_deconvolved_spectra)
        {
          last_deconvolved_spectra.erase(last_deconvolved_spectra.begin());
        }
        last_deconvolved_spectra[ms_level].push_back(deconvolved_spectrum);
      }

      if (merge != 2)
      {
        scan_rt_map[deconvolved_spectrum.getScanNumber()] = it->getRT();
      }

      if (report_dummy)
      {
#pragma omp parallel sections default(none) shared(fd_charge_dummy, fd_noise_dummy, fd_iso_dummy, it, precursor_specs, scan_number, precursor_map_for_real_time_acquisition)
        {
#pragma omp section
          fd_charge_dummy.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_real_time_acquisition);
#pragma omp section
          fd_noise_dummy.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_real_time_acquisition);
#pragma omp section
          fd_iso_dummy.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_real_time_acquisition);
        }
        DeconvolvedSpectrum dummy_deconvolved_spectrum(scan_number);
        deconvolved_spectrum.sortByQscore();
        float qscore_threshold_for_dummy = deconvolved_spectrum[deconvolved_spectrum.size() - 1].getQscore();
        dummy_deconvolved_spectrum.setOriginalSpectrum(*it);
        dummy_deconvolved_spectrum.reserve(fd_iso_dummy.getDeconvolvedSpectrum().size() + fd_charge_dummy.getDeconvolvedSpectrum().size() + fd_noise_dummy.getDeconvolvedSpectrum().size());

        for (auto& pg : fd_charge_dummy.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        for (auto& pg : fd_iso_dummy.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        for (auto& pg : fd_noise_dummy.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        deconvolved_spectrum.sort();
        dummy_deconvolved_spectrum.sort();

        dummy_deconvolved_spectra.push_back(dummy_deconvolved_spectrum);
      }
      qspec_cntr[ms_level - 1]++;
      mass_cntr[ms_level - 1] += deconvolved_spectrum.size();
      deconvolved_spectra.push_back(deconvolved_spectrum);

      elapsed_deconv_cpu_secs[ms_level - 1] += double(clock() - deconv_begin) / CLOCKS_PER_SEC;
      elapsed_deconv_wall_secs[ms_level - 1] += chrono::duration<double>(chrono::high_resolution_clock::now() - deconv_t_start).count();
    
    }




    return 0;
   }



}
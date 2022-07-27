// --------------------------------------------------------------------------
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

namespace OpenMS
{
  DeconvolvedSpectrum::DeconvolvedSpectrum(const MSSpectrum& spectrum, const int scan_number) : scan_number_(scan_number)
  {
    spec_ = spectrum;
  }

  MSSpectrum DeconvolvedSpectrum::toSpectrum(const int to_charge, bool retain_undeconvolved)
  {
    auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);
    if (spec_.getMSLevel() > 1 && precursor_peak_group_.empty())
    {
      return out_spec;
    }
    double charge_mass_offset = abs(to_charge) * FLASHDeconvHelperStructs::getChargeMass(to_charge >= 0);
    std::unordered_set<double> deconvolved_mzs;

    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      out_spec.emplace_back(pg.getMonoMass() + charge_mass_offset, pg.getIntensity());
      if (retain_undeconvolved)
      {
        for (auto& p : pg)
        {
          deconvolved_mzs.insert(p.mz);
        }
      }
    }
    if (retain_undeconvolved)
    {
      for (auto& p : spec_)
      {
        if (deconvolved_mzs.find(p.getMZ()) != deconvolved_mzs.end()) // if p is deconvolved
        {
          continue;
        }
        out_spec.emplace_back(p.getMZ() + charge_mass_offset - FLASHDeconvHelperStructs::getChargeMass(to_charge >= 0), p.getIntensity());
      }
    }
    out_spec.sortByPosition();
    if (!precursor_peak_group_.empty() && !precursor_peak_.empty())
    {
      Precursor precursor(spec_.getPrecursors()[0]);
      // precursor.setCharge((precursor_peak_group_.isPositive() ?
      //                      precursor_peak_group_.getRepAbsCharge() :
      //                      -precursor_peak_group_.getRepAbsCharge()));//getChargeMass
      precursor.setCharge(to_charge);
      precursor.setMZ(precursor_peak_group_.getMonoMass() + charge_mass_offset);
      precursor.setIntensity(precursor_peak_group_.getIntensity());
      out_spec.getPrecursors().clear();
      out_spec.getPrecursors().emplace_back(precursor);
    }
    return out_spec;
  }


  const MSSpectrum& DeconvolvedSpectrum::getOriginalSpectrum() const
  {
    return spec_;
  }

  const PeakGroup& DeconvolvedSpectrum::getPrecursorPeakGroup() const
  {
    if (precursor_peak_group_.empty())
    {
      return *(new PeakGroup());
    }
    return precursor_peak_group_;
  }

  int DeconvolvedSpectrum::getPrecursorCharge() const
  {
    return precursor_peak_.getCharge();
  }

  double DeconvolvedSpectrum::getCurrentMaxMass(const double max_mass) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_mass;
    }
    return precursor_peak_group_.getMonoMass();
  }

  double DeconvolvedSpectrum::getCurrentMinMass(const double min_mass) const
  {
    if (spec_.getMSLevel() == 1)
    {
      return min_mass;
    }
    return 50.0;
  }

  int DeconvolvedSpectrum::getCurrentMaxAbsCharge(const int max_abs_charge) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_abs_charge;
    }
    return abs(precursor_peak_.getCharge());
  }

  const Precursor& DeconvolvedSpectrum::getPrecursor() const
  {
    return precursor_peak_;
  }

  int DeconvolvedSpectrum::getScanNumber() const
  {
    return scan_number_;
  }

  int DeconvolvedSpectrum::getPrecursorScanNumber() const
  {
    return precursor_scan_number_;
  }


  const String& DeconvolvedSpectrum::getActivationMethod() const
  {
    return activation_method_;
  }


  void DeconvolvedSpectrum::setPrecursor(const Precursor& precursor)
  {
    precursor_peak_ = precursor;
  }

  void DeconvolvedSpectrum::setPrecursorIntensity(const double i)
  {
    precursor_peak_.setIntensity(i);
  }

  void DeconvolvedSpectrum::setActivationMethod(const String& method)
  {
    activation_method_ = method;
  }

  void DeconvolvedSpectrum::setPrecursorPeakGroup(const PeakGroup& pg)
  {
    precursor_peak_group_ = pg;
  }

  void DeconvolvedSpectrum::setPrecursorScanNumber(const int scan_number)
  {
    precursor_scan_number_ = scan_number;
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::begin() const noexcept
  {
    return peak_groups.begin();
  }
  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::end() const noexcept
  {
    return peak_groups.end();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::begin() noexcept
  {
    return peak_groups.begin();
  }
  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::end() noexcept
  {
    return peak_groups.end();
  }

  const PeakGroup& DeconvolvedSpectrum::operator[](const Size i) const
  {
    return peak_groups[i];
  }

  void DeconvolvedSpectrum::push_back(const PeakGroup& pg)
  {
    peak_groups.push_back(pg);
  }
  Size DeconvolvedSpectrum::size() const noexcept
  {
    return peak_groups.size();
  }
  void DeconvolvedSpectrum::clear()
  {
    peak_groups.clear();
  }
  void DeconvolvedSpectrum::reserve(Size n)
  {
    peak_groups.reserve(n);
  }
  bool DeconvolvedSpectrum::empty() const
  {
    return peak_groups.empty();
  }
  void DeconvolvedSpectrum::swap(std::vector<PeakGroup>& x)
  {
    peak_groups.swap(x);
  }

  void DeconvolvedSpectrum::sort()
  {
    std::sort(peak_groups.begin(), peak_groups.end());
  }

  void DeconvolvedSpectrum::sortByQScore()
  {
    std::sort(peak_groups.begin(), peak_groups.end(), [](const PeakGroup& p1, const PeakGroup& p2) { return p1.getQScore() > p2.getQScore(); });
  }

  void DeconvolvedSpectrum::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<DeconvolvedSpectrum>& deconvolved_decoy_spectra) // per ms level + precursor update as well.
  {
    std::map<int, std::vector<float>> tscore_map; // per ms level
    std::map<int, std::vector<float>> dscore_map;
    std::map<int, std::map<float, float>> qscore_map; // maps for total qvalues

    std::map<int, std::vector<float>> dscore_with_charge_decoy_only_map;
    std::map<int, std::map<float, float>> qscore_with_charge_decoy_only_map; // maps for charge decoy only qvalues


    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      int ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      if(tscore_map.find(ms_level) == tscore_map.end())
      {
        tscore_map[ms_level] = std::vector<float>();
        dscore_map[ms_level] = std::vector<float>();
        qscore_map[ms_level] = std::map<float, float>();
        dscore_with_charge_decoy_only_map[ms_level] = std::vector<float>();
        qscore_with_charge_decoy_only_map[ms_level] = std::map<float, float>();
      }
      for (auto& pg : deconvolved_spectrum)
      {
        tscore_map[ms_level].push_back(pg.getQScore());
      }
    }
    for (auto& decoy_deconvolved_spectrum : deconvolved_decoy_spectra)
    {
      int ms_level = decoy_deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : decoy_deconvolved_spectrum)
      {
        dscore_map[ms_level].push_back(pg.getQScore());
        if(pg.getDecoyIndex() == 2)
        {
          continue;
        }
        dscore_with_charge_decoy_only_map[ms_level].push_back(pg.getQScore());
      }
    }

    for(auto& titem: tscore_map)
    {
      int ms_level = titem.first;
      auto& tscore = titem.second;
      auto& dscore = dscore_map[ms_level];
      auto& dscore_with_charge_decoy_only = dscore_with_charge_decoy_only_map[ms_level];

      std::sort(tscore.begin(), tscore.end());
      std::sort(dscore.begin(), dscore.end());
      std::sort(dscore_with_charge_decoy_only.begin(), dscore_with_charge_decoy_only.end());

      auto& map = qscore_map[ms_level];
      float tmp_q = 1;
      for (int i = 0; i < tscore.size(); i++)
      {
        float ts = tscore[i];
        int dindex = std::distance(std::upper_bound(dscore.begin(), dscore.end(), ts), dscore.end());
        int tindex = tscore.size() - i;

        tmp_q = std::min(tmp_q, ((float)dindex / tindex));
        map[ts] = tmp_q;
      }

      auto& map_with_charge_decoy_only = qscore_with_charge_decoy_only_map[ms_level];
      float tmp_q_with_charge_decoy_only = 1;
      for (int i = 0; i < tscore.size(); i++)
      {
        float ts = tscore[i];
        int dindex = std::distance(std::upper_bound(dscore_with_charge_decoy_only.begin(), dscore_with_charge_decoy_only.end(), ts), dscore_with_charge_decoy_only.end());
        int tindex = tscore.size() - i;

        tmp_q_with_charge_decoy_only = std::min(tmp_q_with_charge_decoy_only, ((float)dindex / tindex));
        map_with_charge_decoy_only[ts] = tmp_q_with_charge_decoy_only;
      }
    }

    for(auto& titem: qscore_map)
    {
      int ms_level = titem.first;
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if(deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }
        auto& map = qscore_map[ms_level];
        auto& map_with_charge_decoy_only_map = qscore_with_charge_decoy_only_map[ms_level];

        for (auto& pg : deconvolved_spectrum)
        {
          pg.setQvalue(map[pg.getQScore()]);
          pg.setQvalueWithChargeDecoyOnly(map_with_charge_decoy_only_map[pg.getQScore()]);
          if(deconvolved_spectrum.getOriginalSpectrum().getMSLevel() > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
          {
            double qs = deconvolved_spectrum.getPrecursorPeakGroup().getQScore();
            auto& pmap = qscore_map[ms_level - 1];
            auto& pmap_with_charge_decoy_only_map = qscore_with_charge_decoy_only_map[ms_level - 1];
            deconvolved_spectrum.setPrecursorPeakGroupQvalue(pmap[qs], pmap_with_charge_decoy_only_map[qs]);
          }
        }
      }
    }
  }

  void DeconvolvedSpectrum::setPrecursorPeakGroupQvalue(const double qvalue, const double qvalue_with_charge_decoy_only)
  {
    precursor_peak_group_.setQvalue(qvalue);
    precursor_peak_group_.setQvalueWithChargeDecoyOnly(qvalue_with_charge_decoy_only);
  }
} // namespace OpenMS

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



#include<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>

namespace OpenMS
{
  static void filterLowPeaks(MSExperiment& map, Size count)
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




}
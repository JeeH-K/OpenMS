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

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


using namespace OpenMS;

namespace OpenMS
{
  /**
  @brief FLASHDeconv algorithm:ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
  From MSSpectrum, this class outputs DeconvolvedSpectrum.
  */

	class OPENMS_DLLAPI FLASHDeconvAlgorithm
  {
  public:
    static void filterLowPeaks(MSExperiment& map, Size count);

    int SpectrumDeconvolution(MSExperiment& map, ProgressLogger& progresslogger, uint& current_max_ms_level, const int target_precursor_charge, std::vector<size_t>& spec_cntr,
                                                     std::unordered_map<UInt, std::vector<DeconvolvedSpectrum>>& last_deconvolved_spectra, SpectralDeconvolution fd,
                                                     std::map<int, std::vector<std::vector<float>>> precursor_map_for_real_time_acquisition, double target_precursor_mass,
                                                     std::map<int, PeakGroup>& precursor_peak_groups, double topFD_SNR_threshold, double& expected_identification_count, String out_mzml_file,
                                                     int mzml_charge, uint current_min_ms_level, const DoubleList tols, MSExperiment& exp, String out_anno_mzml_file, MSExperiment& exp_annotated,
                                                     int num_last_deconvolved_spectra, const int merge, std::map<int, double>& scan_rt_map, const bool report_dummy,
                                                     SpectralDeconvolution fd_charge_dummy, SpectralDeconvolution fd_noise_dummy, SpectralDeconvolution fd_iso_dummy,
                                                     std::vector<DeconvolvedSpectrum>& dummy_deconvolved_spectra, std::vector<size_t>& qspec_cntr, std::vector<size_t>& mass_cntr,
                                                     std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<double>& elapsed_deconv_cpu_secs,
                                                     std::vector<double>& elapsed_deconv_wall_secs);
  };

    typedef unsigned int uint;
}
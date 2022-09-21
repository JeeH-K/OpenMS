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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <algorithm>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include "boost/dynamic_bitset.hpp"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <fstream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvQuantHelper.h>

using namespace std;
namespace OpenMS
{
  class OPENMS_DLLAPI FLASHDeconvQuantAlgorithm :
      public ProgressLogger,
      public DefaultParamHandler
{
  public:
    typedef FLASHDeconvQuantHelper::FeatureSeed FeatureSeed;
    typedef FLASHDeconvQuantHelper::FeatureGroup FeatureGroup;
    typedef FLASHDeconvQuantHelper::Feature Feature;

    /// Default constructor
    FLASHDeconvQuantAlgorithm();

    /// Default destructor
    ~FLASHDeconvQuantAlgorithm() = default;

    /// copy constructor
    FLASHDeconvQuantAlgorithm(const FLASHDeconvQuantAlgorithm& ) = delete;

    /// move constructor
    FLASHDeconvQuantAlgorithm(FLASHDeconvQuantAlgorithm&& other) = default;

    /// assignment operator
    FLASHDeconvQuantAlgorithm& operator=(const FLASHDeconvQuantAlgorithm& fd);

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace> &input_mtraces, std::vector<FeatureGroup> &out_fgroups);

    // test purpose
    String output_file_path_;

  protected:
    void updateMembers_() override;

  private:
    Param getFLASHDeconvParams_() const;

    // equivalent to FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_
    void getFeatureFromSpectrum_(std::vector<FeatureSeed*> &local_traces, std::vector<FeatureGroup> &local_fgroup, const double &rt);

    void buildMassTraceGroups_(std::vector<FeatureSeed> &in_seeds, std::vector<FeatureGroup> &features);

    bool scoreAndFilterFeatureGroup_(FeatureGroup& fg, double min_iso_score = -1) const;

    void refineFeatureGroups_(std::vector<FeatureGroup>& features);

    bool rescoreFeatureGroup_(FeatureGroup& fg) const;

    void setFeatureGroupScore_(FeatureGroup &fg) const;

    double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const;

    double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const;

    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    bool doFWHMbordersOverlap_(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    bool doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2) const;

    void clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups,
                               std::vector<MassTrace>& input_mtraces);

    void resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                   std::vector<MassTrace> &input_masstraces,
                                   std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   const std::set<Size>& hypo_indices,
                                   std::vector<FeatureGroup>& out_features);

    void setOptionalDetailedOutput_();

    void writeMassTracesOfFeatureGroup_(const FeatureGroup &fgroup,
                                       const Size &fgroup_idx,
                                       const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                       const bool &is_before_resolution);

    void writeTheoreticalShapeForConflictResolution_(const Size &fgroup_idx,
                                                     const FeatureSeed &shared_mt,
                                                     const std::vector<double> &theo_intensities,
                                                     const double &calculated_ratio);

    void resolveConflictRegion_(std::vector<Feature> &conflicting_features,
                                const std::vector<MassTrace> &conflicting_mts,
                                const std::vector<Size> &conflicting_mt_indices);

    MassTrace updateMassTrace_(const MassTrace& ref_trace, const double &ratio) const;

    void fitTraceModelFromUniqueTraces_(Feature const& tmp_feat, EGHTraceFitter* fitted_model) const;

    void runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces, EGHTraceFitter* fitter) const;

    void updateFeatureWithFitModel(std::vector<Feature>& conflicting_features, Size mt_index,
                                   const MassTrace& obs_masstrace, const Size& org_index_of_obs_mt,
                                   Matrix<int>& pointer_to_components, vector<std::vector<double>>& components);

    void getMostAbundantMassTraceFromFeatureGroup_(const FeatureGroup &fgroup,
                                                  const int &ignore_this_charge,
                                                  FeatureSeed* &most_abundant_mt_ptr,
                                                  const std::vector<std::vector<Size>>& shared_m_traces) const;

    void getFLASHDeconvConsensusResult();

    bool isThisMassOneOfTargets_(const double &candi_mass, const double &candi_rt) const;

    void makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const;

    void setFeatureGroupMembersForResultWriting_(std::vector<FeatureGroup> &f_groups) const;

    bool isEligibleFeatureForConflictResolution_(Feature &new_feature, std::vector<std::vector<Size>> &shared_m_traces_indices, FeatureGroup &feat_group) const;

    /// parameter stuff
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    double min_mass_;
    double max_mass_;
    double mz_tolerance_; // ppm
    bool resolving_shared_signal_;

    const double mass_tolerance_da_ = 3; // Da, for feature mass collection
    const double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;

    // advanced parameter?
    Size min_nr_mtraces_ = 3; // minimum number of consecutive bridges among mass traces to support feature
    Size min_nr_peaks_in_mtraces_ = 4; // at least 4 is needed for EGHTraceFitter
    bool use_smoothed_intensities_;
    double rt_window_ = 1; // TODO : remove?

    /// variables for internal use (not for user input)
    FLASHDeconvHelperStructs::PrecalculatedAveragine iso_model_;
    Size max_nr_traces_; // calculated from iso_model_ (setAveragineModel())

    /// cosine threshold between observed and theoretical isotope patterns for MS1
    double min_isotope_cosine_ = 0.90;

    /// FLASHDeconvAlgorithm class for deconvolution
    FLASHDeconvAlgorithm fd_;

    // loop up table
    std::vector<std::pair<double, double>> target_masses_; // mass and rt
    bool with_target_masses_ = false;

    /// for detailed outputs (mostly test purpose)
    bool shared_output_requested_;
    std::fstream shared_out_stream_ = std::fstream();
};
}
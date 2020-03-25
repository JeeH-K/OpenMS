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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <QDirIterator>
#include <QFileInfo>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHDeconvQ TOPP_FLASHDeconvQ
  (Need to be modified)

  @brief  @ref
  @code
  @endcode
  @verbinclude
  @htmlinclude
*/
// We do not want this class to show up in the docu:
// NEED to fill this part later

class TOPPFLASHDeconvQ :
    public TOPPBase
{
public:
  TOPPFLASHDeconvQ() :
      TOPPBase("FLASHDeconvQ", "Proteoform quantification of top-down MS data based on FLASHDeconv",
               false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<input file>", "", "Input file");
    //    //    setValidFormats_("in", ListUtils::create<String>("mzML"),);
    //    registerOutputFile_("out", "<output file prefix/output dir>", "",
    //                        "Output file prefix or output dir (if prefix, [file prefix].tsv will be generated. "
    //                        "if dir, [dir]/[inputfile].tsv is generated per [inputfile])");
  }

  static constexpr double RT_TOLERANCE = 30; // per feature
  static const int CHARGE_TOLERANCE = 3; // per feature


  //  struct OPENMS_DLLAPI AveragineStruct{
  //    IsotopeDistribution isotopes;
  //    double intensityNorm;
  //    double averagine_width;
  //    std::vector<std::pair<double, double>> mz_window_pairs;
  //    std::vector<double> mz_positions; // starting from 0
  //
  //    AveragineStruct(IsotopeDistribution isotope_dist, double pattern_width, double inty_norm,
  //                    std::vector<std::pair<double, double>> mz_windows, std::vector<double> mzs):
  //        isotopes(isotope_dist),
  //        averagine_width(pattern_width),
  //        intensityNorm(inty_norm),
  //        mz_window_pairs(mz_windows),
  //        mz_positions(mzs)
  //    {}
  //  };

  struct OPENMS_DLLAPI FeatureStruct{
    const double average_mass = 0;
    const double rt_start = 0;
    const double rt_end = 0;
    const int min_charge = 0;
    const int max_charge = 0; // +1 than real max_charge (for used in iteration)

    std::vector<double> theoretical_mz_positions; // per charge

    IsotopeDistribution averagine = {}; // mzs are not for use
    double intensityNorm = 0.0;
    // true mz information in averagine
    std::vector<double> iso_peak_deltas; // per charge. Constants::ISOTOPE_MASSDIFF_55K_U/charge
    std::vector<double> averagine_widths; // per charge

    //    std::vector<AveragineStruct> averagines;
    //    std::vector<IsotopeDistribution> theo_averagines;
    //    std::vector<double> theo_averagines_half_width; // half the width

    FeatureStruct();

    ~FeatureStruct()
    {
      std::vector<double>().swap(theoretical_mz_positions);
      std::vector<double>().swap(iso_peak_deltas);
      std::vector<double>().swap(averagine_widths);
    };

    FeatureStruct(double avg_mass, double rt_s, double rt_e, int min_c, int max_c):
        average_mass(avg_mass),
        rt_start(rt_s - RT_TOLERANCE),
        rt_end(rt_e + RT_TOLERANCE),
        min_charge(min_c - CHARGE_TOLERANCE),
        max_charge(max_c + CHARGE_TOLERANCE+1)
    {
      theoretical_mz_positions.reserve(max_charge-min_charge);
      averagine_widths.reserve(max_charge-min_charge);
      iso_peak_deltas.reserve(max_charge-min_charge);
    }

    // for test purpose
    void print_vectors()
    {
      cout << "charge\tios peak deltas\ttheo_avg_peak_pos" << endl;
      for(auto cs_index = 0; cs_index < max_charge-min_charge; cs_index++){
        cout << (max_charge+cs_index) <<"\t"<< to_string(iso_peak_deltas[cs_index]) << "\t" << to_string(theoretical_mz_positions[cs_index]) << endl;
      }

    }
  };

  struct OPENMS_DLLAPI CandidateIsoDist{
    FeatureStruct* fstruct_ptr;
    double score;
    std::vector<Peak1D> found_peaks;

    CandidateIsoDist(FeatureStruct* fs_ptr, double s):
        fstruct_ptr(fs_ptr), score(s)
    {}
  };

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // input file handling (NEED to change later)
    String infilePath = "/Users/jeek/Documents/A4B/FFI_paper/mzml/simple/190226_Cyto_1_FD_500ng.tsv";
    String rawfilePath = "/Users/jeek/Documents/A4B/FFI_paper/mzml/simple/190226_Cyto_1_FD_500ng.mzML";
    Size minCharge = 2;
    Size maxCharge = 100;
    double MIN_RT = std::numeric_limits<double>::max();
    double MAX_RT = 0;
    double mz_tolerance = 100 * 1e-6;
    //    double mz_peak_tolerance = 20 * 1e-6;
    double min_mz_cosine_score = .5;
    //    String outfilePath = ''

    // for checking time
    double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
    double elapsed_deconv_cpu_secs = 0, elapsed_deconv_wall_secs = 0;
    auto begin = clock();
    auto t_start = chrono::high_resolution_clock::now();

    OPENMS_LOG_INFO << "Processing : " << infilePath << endl;

    // ------FLASHDeconv result handling-------
    std::vector<FeatureStruct> fstructs;
    read_FLASHDeconv_result_file(infilePath, fstructs);
    OPENMS_LOG_INFO << fstructs.size() << " features will be used" << endl;

    // calculating theoretical m/z with charge ranges &
    double max_mass = 0; // for calculating Averagine dist. (set the limitation)

    std::vector<std::pair<double, FeatureStruct*>> universal_mzs;
    for (auto& fstr : fstructs){
      if (max_mass < fstr.average_mass){
        max_mass = fstr.average_mass;
      }

      //      cout << "# charge : " << (fstr.max_charge-fstr.min_charge) << endl;
      for (Size cs = fstr.min_charge; cs < fstr.max_charge; cs++){
        auto mz = (fstr.average_mass + cs*Constants::PROTON_MASS_U) / cs;
        fstr.theoretical_mz_positions.push_back(mz);
        fstr.iso_peak_deltas.push_back(Constants::ISOTOPE_MASSDIFF_55K_U/cs);
        universal_mzs.push_back(std::make_pair(mz, &fstr));
        //        cout <<  cs <<  "\t"<< to_string(Constants::ISOTOPE_MASSDIFF_55K_U/cs) <<"\t" << to_string(mz) << endl;
      }
      //      break;

      // calculate rt range
      if (MIN_RT > fstr.rt_start){
        MIN_RT = fstr.rt_start;
      }
      if (MAX_RT < fstr.rt_end){
        MAX_RT = fstr.rt_end;
      }
    }
    sort(universal_mzs.begin(), universal_mzs.end()); // to make a 'pattern' for searching

    //// creating Averagine distribution //// for reducing the calculation time
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(max_mass);
    maxIso.trimRight(0.05 * maxIso.getMostAbundant().getIntensity());
    generator->setMaxIsotope((Size) maxIso.size() - 1);
    auto factor = .02;
    for (auto & fstr : fstructs)
    {
      auto iso = generator->estimateFromPeptideWeight(fstr.average_mass);
      auto inty_cutoff = factor * iso.getMostAbundant().getIntensity();
      iso.trimRight(inty_cutoff);
      iso.trimLeft(inty_cutoff);
      fstr.averagine = iso;

      int charge_size = fstr.max_charge - fstr.min_charge;
      Size avg_peak_delta_size = iso.size() - 1;
      for (auto cs_index = 0; cs_index < charge_size; cs_index++){
        auto cur_cs = fstr.min_charge + cs_index;
        auto isodist_width = fstr.iso_peak_deltas[cs_index] * avg_peak_delta_size;
        fstr.averagine_widths.push_back(isodist_width/cur_cs);
      }
      double intsy_norm = .0;
      for (auto &p : iso)
      {
        intsy_norm += (p.getIntensity() * p.getIntensity());
      }
      fstr.intensityNorm = intsy_norm;
    }
    //      for (auto cur_cs = fstr.min_charge; cur_cs < fstr.max_charge; cur_cs++){
    //      for (auto cur_cs = fstr.min_charge; cur_cs < 300; cur_cs+=50){
    //        auto iso = generator->estimateFromPeptideWeight(Constants::PROTON_MASS_U*cur_cs + fstr.average_mass);
    //      cout << to_string(fstr.average_mass) << endl;

    //        cout << "TITLE:" << to_string(Constants::PROTON_MASS_U*cur_cs + fstr.average_mass) << endl;
    //        auto std = iso[0].getMZ();
    //        for(auto& p : iso){
    //          cout << to_string(p.getMZ()-std) << "\t" << to_string(p.getIntensity()) << endl;
    //        }
    //        continue;

    //        auto inty_cutoff = factor * iso.getMostAbundant().getIntensity();
    //        iso.trimRight(inty_cutoff);
    //        iso.trimLeft(inty_cutoff);
    //
    //        double intsy_norm = .0;
    //        std::vector<std::pair<double, double>> mz_window_pairs; // <window_start, window_end>
    //        std::vector<double> mz_delta_pos;
    //        mz_delta_pos.reserve(iso.size()-1);
    //        mz_window_pairs.reserve(iso.size());

    //        Peak1D* prev_p = nullptr;
    //        double window_size = 0.0;
    //        cout << to_string(fstr.average_mass) << "\t" << cur_cs << endl;
    //        for (auto& p : iso){
    //          intsy_norm += (p.getIntensity() * p.getIntensity());
    //          p.setMZ(p.getMZ()/cur_cs);

    //          if (prev_p==nullptr){
    //            prev_p = &p;
    //            continue;
    //          }
    //          double diff = p.getMZ() - prev_p->getMZ();
    //          mz_delta_pos.push_back(diff);
    //          window_size = (diff/2);
    //          mz_window_pairs.push_back(std::make_pair(prev_p->getMZ()-window_size, prev_p->getMZ()+window_size));
    //
    //          prev_p = &p;
    //        }
    //        mz_window_pairs.push_back(std::make_pair(prev_p->getMZ()-window_size, prev_p->getMZ()+window_size));

    //        auto avg_ptr = new AveragineStruct(iso, iso.averageMass()-iso[0].getMZ(), intsy_norm,
    //                                                        mz_window_pairs, mz_delta_pos);
    //        fstr.averagines.push_back(*avg_ptr);
    //        fstr.intensityNorm = intsy_norm;
    //      }
    //      return EXECUTION_OK;
    //    }

    OPENMS_LOG_INFO << "Processing : " << rawfilePath << endl;

    // -------mzML file handling--------
    MSExperiment map;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.load(rawfilePath, map);
    int ms1Cntr = 0;
    for (auto it = map.begin(); it != map.end(); ++it)
    {
      if (it->getMSLevel() > 1) continue;
      ms1Cntr++;

      if (it->getRT() < MIN_RT || it->getRT() > MAX_RT) continue;

      //// sweep to find existence of theoretical m/z position
      for (auto theo : universal_mzs)
      {
        auto mz_window_size = theo.first * mz_tolerance;
        auto univ_mz_index = it->findNearest(theo.first, mz_window_size);
        if (univ_mz_index == -1) continue; // if no matching mz to DB

        //// mz found, but if cur rt exceed the theo's rt_range, continue
        if(it->getRT() < theo.second->rt_start or it->getRT() > theo.second->rt_end) continue;

        //// change tolerance, if necessary
        auto charge_index = std::distance(theo.second->theoretical_mz_positions.begin(),
                                          find(theo.second->theoretical_mz_positions.begin(),
                                               theo.second->theoretical_mz_positions.end(), theo.first));
        if (theo.second->averagine_widths[charge_index] > mz_window_size){
          //          cout << "iso dist is wider: " <<theo.second->theo_averagines_half_width[theo_avg_index] << "\t" << tol << endl;
          mz_window_size = theo.second->averagine_widths[charge_index];
        }

        //// searching starting & ending indices for mz-window (both should not be 0-intensity)
        auto start_index = it->findNearest(theo.first-mz_window_size, mz_window_size/2);
        if (start_index == -1) continue;
        else{
          MSSpectrum::ConstIterator tmp_it = it->begin();
          std::advance(tmp_it, start_index);
          while( tmp_it->getIntensity() == 0 ){
            start_index++;
            tmp_it++;
          }
        }
        auto end_index = it->findNearest(theo.first+mz_window_size, mz_window_size/2);
        if (end_index == -1) continue;
        else{
          MSSpectrum::ConstIterator tmp_it = it->begin();
          std::advance(tmp_it, end_index);
          while( tmp_it->getIntensity() == 0 ){
            end_index--;
            tmp_it--;
          }
        }
        end_index++;

        //// if not many peaks were found within window, continue (CHECK IF NUMBERS ARE ADEQUATE)
        if ( end_index-start_index < (theo.second->averagine.size()/2) ) continue;

        ///////// test writing part //////
        cout << "feature: " << to_string(theo.first) << "\t" << to_string(theo.second->average_mass) << "\t" << (theo.second->min_charge + charge_index) << endl;
        //        cout << tol << endl;
        cout << "avg iso : " << to_string(theo.second->averagine.getMostAbundant().getMZ()) <<
             "\t"<<  to_string(theo.second->averagine.getMostAbundant().getIntensity()) << endl;
        //        for (auto& peak : theo.second->averagine){
        //          cout << to_string(peak.getMZ()) << "\t" << to_string(peak.getIntensity()) << endl;
        //        }
        cout << "averagine width : " << to_string(theo.second->averagine_widths[charge_index]) << endl;

        cout << "RT: " << to_string(it->getRT()) << endl;
        cout << "indices : " << univ_mz_index << "\t" << start_index << "\t" << end_index << "\t" << theo.second->averagine.size() << endl;
        for (auto tmp_iter = it->begin()+start_index; tmp_iter!= it->begin()+end_index; tmp_iter++){
          cout << to_string(tmp_iter->getMZ()) << "\t" << to_string(tmp_iter->getIntensity()) << endl;
        }
        //////////////////////////////////

        ////// search for the pattern with maximum cosine score

        /// 1. Cosine similarity for mz positions
        /// Find peaks for calculating intensity-level cosine similarity
        for (auto candi_peak_iter = it->begin()+start_index; candi_peak_iter != it->begin()+end_index; candi_peak_iter++)
        {
          if (candi_peak_iter->getIntensity() <= 0 ) continue;

          auto std_mz = candi_peak_iter->getMZ(); // mz for matching first index of Averagine
          cout << "candi peak : " << std_mz << endl;

          std::vector<std::vector<double>> candi_mzs_vec;
          std::vector<double> mz_cosine_scores;

          // sub isodists
          int max_sub_isodists_size = theo.second->averagine.size()-2; // if #peaks of iso dist < 3, no worth matching
          int iso_deltas_size = theo.second->averagine.size()-1;
          for (int sub_iso_starting_index=0; sub_iso_starting_index < max_sub_isodists_size; sub_iso_starting_index++)
          {
            //            vector<double> candi_mz_pos;
            //            vector<double> matched_isos;
            auto tmp_mz = std_mz;
            cout << "standard mz on spec : " << std_mz << endl;
            //            candi_mz_pos.push_back(tmp_mz);
            int size_of_sub_iso = iso_deltas_size-sub_iso_starting_index; // excluding first iso
            // mzs excluding starting point (already matched)
            std::vector<double> candi_mzs;
            std::vector<double> isos_mzs;
            auto iso_delta = theo.second->iso_peak_deltas[charge_index];
            auto peak_tol = iso_delta/2;
            int matched_counter = 1; // including starting point
            cout << "peak_tolearnce : " << peak_tol << endl;

            // iterate through sub iso
            for (int j=0; j < size_of_sub_iso; j++ )
            {
              tmp_mz += iso_delta;
              isos_mzs[j] = tmp_mz;
              //              cout << (sub_iso_starting_index+j+1) <<  "-th iso candi :" << tmp_mz << endl;
              //              auto found_idx = it->findNearest(tmp_mz, peak_tol);
              auto found_idx = it->findHighestInWindow(tmp_mz, peak_tol, peak_tol);
              if (found_idx == -1){
                candi_mzs[j] = 0.0;
                continue;
              }
              //              if (found_idx != found_h){
              //                auto left = it->MZBegin(tmp_mz - peak_tol);
              //                auto right = it->MZEnd(tmp_mz + peak_tol);
              //                cout << "org mz:" << to_string(tmp_mz) << ", found: " << to_string(it->operator[](found_idx).getMZ())
              //                      << ", highest:" << to_string(it->operator[](found_h).getMZ()) << endl;
              //                for (auto k = left; k<right+1; k++)
              //                {
              //                  cout << to_string(k->getMZ()) <<"\t" << to_string(k->getIntensity())<< endl;
              //                }
              //                return EXECUTION_OK;
              //              }


              Peak1D found_peak = it->operator[](found_idx);
              // if int 0 = countinue
              if (found_peak.getIntensity()==0){
                candi_mzs[j] = 0.0;
                continue;
              }
              matched_counter++;
              candi_mzs[j] = found_peak.getMZ();
            }
            //// CHANGE : maybe other than 3...
            if (matched_counter < 3){
              continue;
            }

            //            cout << "found : " << to_string(theo.second->average_mass) << "\t@RT-" <<  to_string(it->getRT()) << endl;
            //            cout << "matched peaks: " << matched_counter << endl;
            //            for (auto x = 0; x< size_of_sub_iso; x++)
            //            {
            //              cout << to_string(candi_mzs[x]) << endl;
            //            }
            //            return EXECUTION_OK;
            //// mz cosine similarity
            candi_mzs_vec.push_back(candi_mzs);
            mz_cosine_scores.push_back(getCosine(candi_mzs, isos_mzs));
          }

          //// matched peak 모아서... '같은 묶음'이 보이는지 확인. matched peaks(candi_mz)의 구성이 비슷하다거나...
          //// check if candi_mzs_vec has duplicates



          //// intensity cosine similarity
        }



        //        // use getCosine.
        //        std::vector<double> candidatePeakMzs;
        //        candidatePeakMzs.reserve(end_index-start_index);
        //        auto tmp_iter = it->begin()+start_index;
        //        double prev_mz = tmp_iter->getMZ();
        //        for (++tmp_iter ; tmp_iter!= it->begin()+end_index; tmp_iter++){
        //          auto diff = tmp_iter->getMZ() - prev_mz;
        //          candidatePeakMzs.push_back(diff);
        //          prev_mz = tmp_iter->getMZ();
        //        }
        //
        //        int offset = 0;
        //        double maxCosine = -1;
        //        int isoSize = avg_iso_ptr->size();
        //        int min_candi_index = 0;
        //        int max_candi_index = end_index - start_index;
        //        for (int f = -isoSize; f <= max_candi_index; f++)
        //        {
        //          double cos;
        //          cos = getCosine(candidatePeakMzs,
        //                          min_candi_index,
        //                          max_candi_index,
        //                          avg_ptr->mz_deltas,
        //                          isoSize,
        //                          avg_ptr->intensityNorm,
        //                          f);
        //
        //          if (maxCosine <= cos)
        //          {
        //            maxCosine = cos;
        //            offset = f;
        //          }
        //        }
        //
        //        cout << to_string(maxCosine) << "\t" << offset << endl;





        //        return EXECUTION_OK;

        ///// 2. Cosine similarity for intensity

        //        //// testing.... if getCosine works
        //
        //        int perIsotopeIntensitiesSize = end_index - start_index;
        //        double *perIsotopeIntensities = new double[perIsotopeIntensitiesSize];
        //        int i = 0;
        //        for (auto tmp_iter = it->begin()+start_index; tmp_iter!= it->begin()+end_index; tmp_iter++){
        //          cout << to_string(tmp_iter->getMZ()) << "\t" << to_string(tmp_iter->getIntensity()) << endl;
        //          perIsotopeIntensities[i] = tmp_iter->getIntensity();
        //          i++;
        //        }
        //
        //        int offset = 0;
        //        double maxCosine = -1;
        //        int maxIsotopeIndex = 0, minIsotopeIndex = -1;
        //        int isoSize = avg_iso_ptr->size();
        //
        //        for (int i = 0; i < perIsotopeIntensitiesSize; i++)
        //        {
        //          if (perIsotopeIntensities[i] <= 0)
        //          {
        //            continue;
        //          }
        //          maxIsotopeIndex = i;
        //          if (minIsotopeIndex < 0)
        //          {
        //            minIsotopeIndex = i;
        //          }
        //        }
        //
        //        for (int f = -isoSize + minIsotopeIndex; f <= maxIsotopeIndex; f++)
        //        {
        //          double cos;
        //          cos = getCosine(perIsotopeIntensities,
        //                                            minIsotopeIndex,
        //                                            maxIsotopeIndex,
        //                                            *avg_iso_ptr,
        //                                            isoSize,
        //                                            avg_ptr->intensityNorm,
        //                                            f);
        //
        //          if (maxCosine <= cos)
        //          {
        //            maxCosine = cos;
        //            offset = f;
        //          }
        //        }
        //
        //        cout << to_string(maxCosine) << "\t" << offset << endl;
      }

    }
    OPENMS_LOG_INFO << "Used "<< ms1Cntr << " MS1 scans" << endl;

    return EXECUTION_OK;
  }

  void read_FLASHDeconv_result_file(String infilePath, vector<FeatureStruct>& fstruct_vec){
    std::ifstream infile(infilePath);
    String line;

    // handling header
    std::getline(infile, line);
    line.trim();
    vector<String> header_vec;
    line.split("\t", header_vec);

    String mass_col_name = "AverageMass";
    String rts_col_name = "StartRetentionTime";
    String rte_col_name = "EndRetentionTime";
    String minC_col_name = "MinCharge";
    String maxC_col_name = "MaxCharge";
    Size index_of_mass = -1;
    Size index_of_rts = -1;
    Size index_of_rte = -1;
    Size index_of_minc = -1;
    Size index_of_maxc = -1;
    for (Size i = 0; i < header_vec.size(); i++){
      auto col = header_vec[i];
      if (col == mass_col_name){
        index_of_mass = i;
      }else if (col == rts_col_name){
        index_of_rts = i;
      }else if (col == rte_col_name){
        index_of_rte = i;
      }else if (col == minC_col_name){
        index_of_minc = i;
      }else if (col == maxC_col_name){
        index_of_maxc = i;
      }
    }

    // read results and store them in fstruct
    while (std::getline(infile, line)) {
      line.trim();
      vector<String> results;
      line.split("\t", results);
      FeatureStruct* tmp_fstruct = new FeatureStruct(
          String(results[index_of_mass]).toDouble(), String(results[index_of_rts]).toDouble(),
          String(results[index_of_rte]).toDouble(), String(results[index_of_minc]).toInt(),
          String(results[index_of_maxc]).toInt());
      fstruct_vec.push_back(*tmp_fstruct);
    }
  }


  double getCosine(const double *a,
                   int &aStart,
                   int &aEnd,
                   IsotopeDistribution &b,
                   int &bSize,
                   double &bNorm,
                   int offset)
  {
    double n = .0, d1 = .0;

    for (int j = aStart; j < aEnd; j++)
    {
      d1 += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= bSize)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }

    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  double getCosine(std::vector<double> &a,
                   int &aStart,
                   int &aEnd,
                   std::vector<double> &b,
                   int &bSize,
                   double &bNorm,
                   int offset)
  {
    double n = .0, d1 = .0;

    for (int j = aStart; j < aEnd; j++)
    {
      d1 += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= bSize)
      {
        continue;
      }
      n += a[j] * b[i]; //
    }

    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  double getCosine(vector<double> &a, vector<double> &b)
  {
    double n = .0, d1 = .0, d2 = .0;
    for (Size j = 0; j < a.size(); j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
    }

    double d = (d1 * d2);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }
};

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPFLASHDeconvQ tool;
  return tool.main(argc, argv);
}

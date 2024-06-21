// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <utility>

namespace OpenMS
{
FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(): DefaultParamHandler("FLASHExtenderAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHExtenderAlgorithm& FLASHExtenderAlgorithm::operator=(const FLASHExtenderAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHExtenderAlgorithm::setDefaultParams_()
{
  defaults_.setValue("max_mass_shift", 10.0, "Maximum mass shift for blind search.");
  defaults_.setValue("num_blind", 1, "Number of blind modification per protein.");

  defaults_.setValue("fdr", 1.0, "Protein FDR threshold.");
  defaults_.setMaxFloat("fdr", 1.0);
  defaults_.setMinFloat("fdr", 0.01);

  defaults_.setValue("keep_decoy", "false", "Keep decoy proteins.");
  defaults_.addTag("keep_decoy", "advanced");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  // registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
  // setValidStrings_("fragment_type", ListUtils::create<String>("full,internal,N-terminal,C-terminal,a-ion,b-ion,c-ion,x-ion,y-ion,z-ion"));

  defaultsToParam_();
}

void FLASHExtenderAlgorithm::updateMembers_()
{
  max_mod_cntr_ = param_.getValue("num_blind");
  max_mod_mass_ = param_.getValue("max_mass_shift");
}

Size FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod) const
{
  return ((node_index * (pro_length_ + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
         + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_);
}

int FLASHExtenderAlgorithm::getNodeIndex_(Size vertex) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1)) / (pro_length_ + 1);
}

int FLASHExtenderAlgorithm::getProIndex_(Size vertex) const
{
  return ((vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % (pro_length_ + 1);
}

int FLASHExtenderAlgorithm::getScore_(Size vertex) const
{
  return vertex % (max_path_score_ - min_path_score_ + 1) + min_path_score_;
}

int FLASHExtenderAlgorithm::getModNumber_(Size vertex) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1)) % (max_mod_cntr_ + 1);
}


void FLASHExtenderAlgorithm::run_(const FLASHTaggerAlgorithm& tagger,
                                  const ProteinHit& hit,
                                  std::vector<std::vector<Size>>& all_paths,
                                  int mode) // per hit
{
  const auto& node_scores = node_score_map_[mode];
  const auto& node_masses = node_mass_map_[mode];

  pro_masses_.clear();
  pro_masses_.reserve(hit.getSequence().size() + 1);
  pro_masses_.push_back(0);

  auto seq = hit.getSequence();
  pro_masses_.reserve(seq.size());

  if (mode == 0) seq = seq.reverse();
  for (const auto& aa : seq)
  {
    pro_masses_.push_back(pro_masses_.back() + AASequence::fromString(aa).getMonoWeight(Residue::Internal));
  }
  pro_length_ = pro_masses_.size();

  std::vector<int> tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends;
  std::vector<FLASHHelperClasses::Tag> tags;
  std::set<Size> sinks;

  tagger.getTagsMatchingTo(hit, tags);

  // std::sort(tags.begin(), tags.end(), [](const FLASHHelperClasses::Tag& v1, const FLASHHelperClasses::Tag& v2)
  //      {return v1.getScore() > v2.getScore();});


  FLASHHelperClasses::DAG dag((1 + node_scores.size()) * (1 + pro_length_) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_));

  for (const auto& tag : tags)
  {
    if (mode == 0 && tag.getCtermMass() <= 0) continue;
    if (mode == 1 && tag.getNtermMass() <= 0) continue;
    if (mode == 2 && precursor_mass_ <= 0) continue;
    std::vector<int> positions;
    std::vector<double> masses;

    tagger.getMatchedPositionsAndFlankingMassDiffs(positions, masses, hit, tag);
    auto tag_masses = tag.getMzs();
    std::sort(tag_masses.begin(), tag_masses.end());

    std::vector<double> start_masses, end_masses;
    if (tag.getCtermMass() >= 0) // suffix
    {
      for (const auto& shift : suffix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;

        if (precursor_mass_ > 0)
        {
          start_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - start_mass;
          end_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - end_mass;
        }

        start_masses.push_back(start_mass);
        end_masses.push_back(end_mass);
      }
    }
    else
    {
      for (const auto& shift : prefix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;
        start_masses.push_back(start_mass);
        end_masses.push_back(end_mass);
      }
    }

    for (int l = 0; l < start_masses.size(); l++)
    {
      double delta_start = start_masses[l] * tol_ * 2;
      auto lower_start = std::lower_bound(node_masses.begin(), node_masses.end(), start_masses[l] - delta_start);

      double delta_end = end_masses[l] * tol_ * 2;
      auto lower_end = std::lower_bound(node_masses.begin(), node_masses.end(), end_masses[l] - delta_end);

      int highest_score = min_path_score_;
      int highest_score_start = -1;
      int highest_score_end = -1;

      while (lower_start != node_masses.end())
      {
        if (std::abs(*lower_start - start_masses[l]) < delta_start)
        {
          int index = std::distance(node_masses.begin(), lower_start);
          if (index < node_scores.size() && highest_score < node_scores[index])
          {
            highest_score_start = index;
            highest_score = node_scores[index];
          }
        }
        lower_start++;
      }

      highest_score = min_path_score_;

      while (lower_end != node_masses.end())
      {
        if (std::abs(*lower_end - end_masses[l]) < delta_end)
        {
          Size index = std::distance(node_masses.begin(), lower_end);
          if (index < node_scores.size() && highest_score < node_scores[index])
          {
            highest_score_end = index;
            highest_score = node_scores[index];
          }
        }
        lower_end++;
      }

      if (highest_score_start < 0 || highest_score_end < 0) continue;
      for (int pos : positions)
      {
        if (tag.getCtermMass() >= 0) // suffix
        {
          if (precursor_mass_ == 0)
          {
            pos = pro_length_ - 1 - pos; // invert pos
          }
        }

        if (mode == 0) // suffix inverted
        {
          if (pos - tag.getLength() >= 0 && pos < node_scores.size())
          {
            tag_pro_starts.push_back(pos - tag.getLength());
            tag_pro_ends.push_back(pos); // this can be much faster...
            tag_node_starts.push_back(highest_score_start);
            tag_node_ends.push_back(highest_score_end);
          }
        }
        else
        {
          if(pos >= 0 && pos + tag.getLength() < node_scores.size())
          {
            tag_pro_starts.push_back(pos);
            tag_pro_ends.push_back(pos + tag.getLength()); // this can be much faster...
            tag_node_starts.push_back(highest_score_start);
            tag_node_ends.push_back(highest_score_end);
          }
        }
      }
    }
  }


  // std::cout<<"start points done"<<std::endl;
  constructDAG_(dag, sinks, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends, mode);
  // std::cout<<"dag construction done"<<std::endl;
  Size src = getVertex_(0, 0, 0, 0);
  std::set<Size> best_sinks;
  int max_score = 0;
  int min_mod = max_mod_cntr_;
  for (Size sink : sinks)
  {
    if (sink == src || getScore_(sink) < max_score) continue; // getModNumber_(sink) == 0 ||
    max_score = getScore_(sink);
    min_mod = std::min(min_mod, getModNumber_(sink));
  }

  for (Size sink : sinks)
  {
    if (sink == src || getScore_(sink) < max_score) continue; // || getModNumber_(sink) > min_mod
    best_sinks.insert(sink);
  }
  // #pragma omp critical
  for (const auto& best_sink : best_sinks)
  {
    std::vector<std::vector<Size>> all_paths_per_sink;
    dag.findAllPaths(best_sink, src, all_paths_per_sink, 0);
    all_paths.insert(all_paths.end(), all_paths_per_sink.begin(), all_paths_per_sink.end());
  }
}

void FLASHExtenderAlgorithm::run(const FLASHTaggerAlgorithm& tagger)
{
  setLogType(CMD);

  ion_types_str_ = std::vector<String>({"b-ion", "c-ion", "y-ion", "z-ion"}); // why ion type order matter? fix TODO
  tol_ = 5e-6;

  std::vector<ProteinHit> hits;
  std::vector<double> mzs;
  std::vector<int> scores;
  tagger.getProteinHits(hits);

  startProgress(0, (SignedSize)std::min(5, (int)hits.size()), "running FLASHExtender ...");

  auto spec = tagger.getSpectrum();

  Residue empty;
  for (const auto& ion_str : ion_types_str_)
  {
    if (ion_str == "a-ion") { prefix_shifts_.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b-ion") { prefix_shifts_.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c-ion") { prefix_shifts_.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x-ion") { suffix_shifts_.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y-ion") { suffix_shifts_.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z-ion") { suffix_shifts_.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else
    {
      continue; // TODO warn
    }
  }

  // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
  for (int mode = 0; mode < 3; mode++) // TODO for each ion type combination
  {
    MSSpectrum merged_spec;
    merged_spec.reserve(spec.size() * ion_types_str_.size());

    for (const auto& p : spec)
    {
      if (mode == 0)
      {
        for (const auto& shift : suffix_shifts_)
        {
          if (p.getMZ() - shift <= 0) continue;
          merged_spec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
      }
      else if (mode == 1)
      {
        for (const auto& shift : prefix_shifts_)
        {
          if (p.getMZ() - shift <= 0) continue;
          merged_spec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
      }
      else if (precursor_mass_ > 0)
      {
        for (const auto& shift : prefix_shifts_)
        {
          if (p.getMZ() - shift <= 0) continue;
          merged_spec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
        for (const auto& shift : suffix_shifts_)
        {
          if (precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - (p.getMZ() - shift) <= 0) continue;
          merged_spec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - (p.getMZ() - shift), p.getIntensity());
          // TODO check after all
        }
      }
    }

    merged_spec.sortByPosition();

    // Assign the sorted values back to the original vectors
    std::vector<int> node_scores;
    std::vector<double> node_masses;

    node_scores.reserve(1 + merged_spec.size());
    node_masses.reserve(1 + merged_spec.size());

    node_scores.push_back(0);
    node_masses.push_back(0);

    for (const auto& p : merged_spec)
    {
      double mass = p.getMZ();
      int score = (int)p.getIntensity();

      // if (score > 0) max_path_score_ += score;
      // else min_path_score_ += score;

      double prev_mass = node_masses.back();
      int prev_score = node_scores.back();
      if (mass - prev_mass < 2 * tol_ * mass) // they are the same
      {
        if (node_masses.size() > 1)
        {
          node_masses.pop_back();
          node_scores.pop_back();
        }
        score = std::max(0, score);
        score += std::max(0, prev_score);
        //score += prev_score;
      }
      node_masses.push_back(mass);
      node_scores.push_back(score);
    }
    node_score_map_[mode] = node_scores;
    node_mass_map_[mode] = node_masses;
  }
  max_path_score_ = 300;
  min_path_score_ = -100;

  // #pragma omp parallel for default(none) shared(tagger, hits, std::cout)
  for (int i = 0; i < std::min(5, (int)hits.size()); i++)
  {
    nextProgress();
    auto hit = hits[i];

    for (int mode = 0; mode < 3; mode++)
    {
      if (mode == 2 && precursor_mass_ <= 0) continue;
      std::vector<std::vector<Size>> all_paths;
      // std::vector<std::vector<Size>> all_paths_per_mode;
      run_(tagger, hit, all_paths, mode);
      // maybe filter out paths

      // std::sort(all_paths.begin(), all_paths.end(), [](const std::vector<Size>& v1, const std::vector<Size>& v2) {return v1.size() > v2.size();

      for (const auto& path : all_paths)
      {
        for (Size v : path)
        {
          std::cout << getNodeIndex_(v) << " " << node_mass_map_[mode][getNodeIndex_(v)] << " " << pro_masses_[getProIndex_(v)] << " "
                    << getModNumber_(v) << " " << getScore_(v) << std::endl;
        }
        std::cout << hit.getAccession() << std::endl;
        break;
      }
    }
  } // add positive proteoforms all?
  endProgress();
}

void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                           std::set<Size>& sinks,
                                           std::vector<int>& tag_node_starts,
                                           std::vector<int>& tag_pro_starts,
                                           std::vector<int>& tag_node_ends,
                                           std::vector<int>& tag_pro_ends,
                                           int mode)
{
  auto visited = boost::dynamic_bitset<>(dag.size());
  if (tag_node_starts.empty()) return;

  Size src = getVertex_(0, 0, 0, 0);

  visited[src] = true;
  std::set<Size> visited_tag_edges;
  connectBetweenTags(dag, visited, visited_tag_edges, sinks, src, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends, mode);
}

void FLASHExtenderAlgorithm::connectBetweenTags(FLASHHelperClasses::DAG& dag,
                                                boost::dynamic_bitset<>& visited,
                                                std::set<Size>& visited_tag_edges,
                                                std::set<Size>& sinks,
                                                Size vertex,
                                                std::vector<int>& tag_node_starts,
                                                std::vector<int>& tag_pro_starts,
                                                std::vector<int>& tag_node_ends,
                                                std::vector<int>& tag_pro_ends,
                                                int mode)
{
  int node_index = getNodeIndex_(vertex);
  int pro_index = getProIndex_(vertex);

  int tag_start_index = -1;
  int tag_end_index = -1;

  for (int i = 0; i < tag_node_starts.size(); i++)
  {
    if (tag_start_index < 0 && tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) { tag_start_index = i; }
    if (tag_end_index < 0 && tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) { tag_end_index = i; }
  }

  Size src = getVertex_(0, 0, 0, 0);

  if (tag_start_index >= 0) // within tag
  {
    int node_end = -1;
    int i = tag_start_index;
    int pro_end = -1;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end < tag_pro_ends[i]) pro_end = tag_pro_ends[i];
      i++;
    }
    i = tag_start_index;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end == tag_pro_ends[i]) { node_end = std::max(tag_node_ends[i], node_end); }
      i++;
    }
    //pro_end = tag_pro_ends[tag_start_index];
    //node_end = tag_node_ends[tag_start_index];

    //int mod_num = getModNumber_(vertex);

    //if (visited_tag_edges.find(getVertex_(node_end, pro_end, 1, mod_num)) != visited_tag_edges.end()) return;

    //visited_tag_edges.insert(getVertex_(node_end, pro_end, 1, mod_num));
    std::set<Size> next_vertices;

    extendBetweenTags(dag, visited, next_vertices, vertex, node_end, pro_end, 1e5, mode);

    for (auto next_vertex : next_vertices)
    {
      if (next_vertex == src) continue;
      connectBetweenTags(dag, visited, visited_tag_edges, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends, mode);
    }
  }
  else if (vertex == src || tag_end_index >= 0) // between tag.
  {
    std::set<Size> reachable_vertices;

    for (int tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
    {
      int node_start = tag_node_starts[tag_index];
      int pro_start = tag_pro_starts[tag_index];

      if (node_index > node_start || pro_index > pro_start) continue;
      bool is_visited_start = false;

      int mod_num = getModNumber_(vertex);
      for (int mn = 0; mn <= mod_num; mn++)
        if (visited_tag_edges.find(getVertex_(node_start, pro_start, 0, mn)) != visited_tag_edges.end()){
          is_visited_start = true;
          break;
        }

      if (is_visited_start) continue;
      visited_tag_edges.insert(getVertex_(node_start, pro_start, 0, mod_num));

      std::set<Size> next_vertices;
      extendBetweenTags(dag, visited, next_vertices, vertex, node_start, pro_start, 0, mode);

      //std::cout <<node_start << " and " << pro_start << ": " << node_mass_map_[mode][node_index] << " to " <<  node_mass_map_[mode][node_start] << " " << pro_masses_[pro_index] << " to " <<  pro_masses_[pro_start] << " " << next_vertices.size() << std::endl;

      for (auto next_vertex : next_vertices)
      {
        if (next_vertex == src) continue;
        connectBetweenTags(dag, visited, visited_tag_edges, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends, mode);
        reachable_vertices.insert(next_vertex);
      }
    }
    if (vertex != src && reachable_vertices.empty())
    {
      extendBetweenTags(dag, visited, sinks, vertex, -1, -1, 1e5, mode); //
    }
  }
}

void FLASHExtenderAlgorithm::extendBetweenTags(FLASHHelperClasses::DAG& dag,
                                               boost::dynamic_bitset<>& visited,
                                               std::set<Size>& sinks,
                                               Size vertex,
                                               int node_index,
                                               int pro_index,
                                               int diagonal_counter,
                                               int mode)
{
  if (! visited[vertex]) return;
  int node_index1 = getNodeIndex_(vertex);
  int pro_index1 = getProIndex_(vertex);
  int score1 = getScore_(vertex);
  int num_mod1 = getModNumber_(vertex);
  if (num_mod1 == max_mod_cntr_) diagonal_counter = 1e5;

  if (vertex == getVertex_(0, 0, 0, 0))
  {
    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      Size vertex2 = getVertex_(node_index1, pro_i, score1, num_mod1);
       dag.addEdge(vertex2, vertex, visited);
      if (! visited[vertex2]) continue;
      extendBetweenTags(dag, visited, sinks, vertex2, node_index, pro_index, diagonal_counter, mode);
    }
  }

  double delta_mass1 = node_mass_map_[mode][node_index1] - pro_masses_[pro_index1];
  double delta_mass;

  if (node_index >= 0)
  {
    double node_mass = node_mass_map_[mode][node_index];
    delta_mass = node_mass - pro_masses_[pro_index];
    double margin = 2 * tol_ * node_mass;
    if (std::abs(delta_mass - delta_mass1) > max_mod_mass_ * (max_mod_cntr_ - num_mod1) + margin) { return; }
    if (diagonal_counter > 0)
    {
      if (std::abs(delta_mass - delta_mass1) > margin) return;
    }
  }

  for (int score = score1 + 1; score <= max_path_score_; score++)
  {
    // if the starting point has already taken by a higher scoring path, don't go further.
    Size higher_score_vertex = getVertex_(node_index1, pro_index1, score, num_mod1);
    if (visited[higher_score_vertex]) return;
  }

  if (node_index < 0 && diagonal_counter > 0)
  {
    sinks.insert(vertex);
    pro_index = std::min(pro_index1 + 10, pro_length_ - 1); // if sink is not specified, stretch up to ten amino acids.
  }
  else if (node_index1 > node_index || pro_index1 > pro_index)
    return;

  if (node_index1 == node_index && pro_index1 == pro_index)
  {
    sinks.insert(vertex);
    return;
  }

  const auto& node_scores = node_score_map_[mode];
  const auto& node_masses = node_mass_map_[mode];

  for (int node_i = node_index1 + 1; node_i <= (node_index < 0 ? (int)node_scores.size() - 1 : node_index); node_i++)
  {
    if (node_index < 0 && node_scores[node_i] < 0) continue;
    int score = score1 + node_scores[node_i];
    double t_node_mass = node_masses[node_i];

    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      double t_delta_mass = t_node_mass - pro_masses_[pro_i];
      double t_margin = 2 * tol_ * t_node_mass;

      if (t_delta_mass - delta_mass1 > max_mod_mass_ + t_margin) continue;
      if (delta_mass1 - t_delta_mass > max_mod_mass_ + t_margin) break;

      int num_mod = num_mod1;
      if (std::abs(t_delta_mass - delta_mass1) > t_margin) { num_mod++; }

      if (diagonal_counter > 0 && num_mod != num_mod1) continue; //
      if (num_mod > max_mod_cntr_) continue;
      if (node_index >= 0 && std::abs(t_delta_mass - delta_mass) > max_mod_mass_ * (max_mod_cntr_ - num_mod) + t_margin) continue;

      Size next_vertex = getVertex_(node_i, pro_i, score, num_mod);

      dag.addEdge(next_vertex, vertex, visited);
      if (! visited[next_vertex]) continue;

      int next_diagonal_counter = diagonal_counter;
      if (diagonal_counter > 0) next_diagonal_counter--;
      else if (num_mod != num_mod1)
        next_diagonal_counter = 1;

      extendBetweenTags(dag, visited, sinks, next_vertex, node_index, pro_index, next_diagonal_counter, mode);
    }
  }
}
} // namespace OpenMS

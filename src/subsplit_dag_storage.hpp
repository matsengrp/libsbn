// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

// The following classes are used internally by SubsplitDAG for storing the nodes and
// edges, and providing convenient views and lookups into the data. Terminology has been
// changed in order to distinguish from the public API of SubsplitDAG - Edge becomes
// DAGLine and Node becomes DAGVertex.
//
// Functionality interface of the individual components has been separated from the
// storage containers to allow for flexible representation of multi-element collections:
// class DAGLineStorage - owns the data representing an Edge (IDs of the two nodes, edge
// ID, clade) class DAGLineView - cheaply copyable wrapper for a reference to
// DAGLineStorage class DAGLine - a pure interface providing accessor methods to
// DAGLine{View,Storage} aliases LineView and ConstLineView - used within SubsplitDAG to
// access const and mutable edges class GenericLinesView - cheaply copyable view into a
// collection of edges
//
// Similar design is applied to nodes - they have storage, view, interface and
// collection classes. The node view class is called GenericSubsplitDAGNode and is
// defined in subsplit_dag_node.hpp. Additionally nodes exposes class
// GenericNeighborsView, which is a view that wraps a reference to the node neighbors
// collection - internally a std::map<VertexId, LineId>.
//
// The SubsplitDAGStorage class glues together all components, and is owned exclusively
// by SubsplitDAG.

#pragma once

#include <limits>
#include <vector>
#include <map>
#include <optional>
#include <functional>
#include <utility>
#include <type_traits>
#include <memory>

#include "bitset.hpp"

enum class Direction { Rootward, Leafward };
enum class Clade { Unspecified, Left, Right };

inline Clade Opposite(Clade clade) {
  switch (clade) {
    case Clade::Left:
      return Clade::Right;
    case Clade::Right:
      return Clade::Left;
    default:
      Failwith("Use of unspecified clade");
  }
}

using VertexId = size_t;
using LineId = size_t;
static constexpr size_t NoId = std::numeric_limits<size_t>::max();

template <typename Derived>
class DAGLine {
 public:
  LineId GetId() const { return storage().id_; }
  VertexId GetParent() const { return storage().parent_; }
  VertexId GetChild() const { return storage().child_; }
  Clade GetClade() const { return storage().clade_; }

  Derived& SetId(LineId id) {
    storage().id_ = id;
    return derived();
  }
  Derived& SetParent(VertexId id) {
    storage().parent_ = id;
    return derived();
  }
  Derived& SetChild(VertexId id) {
    storage().child_ = id;
    return derived();
  }
  Derived& SetClade(Clade clade) {
    storage().clade_ = clade;
    return derived();
  }

  std::pair<VertexId, VertexId> GetVertexIds() const {
    return {GetParent(), GetChild()};
  }

 private:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
  auto& storage() { return derived().storage(); }
  const auto& storage() const { return derived().storage(); }
};

template <typename T>
class DAGLineView;

class DAGLineStorage : public DAGLine<DAGLineStorage> {
 public:
  DAGLineStorage() = default;
  DAGLineStorage(const DAGLineStorage&) = default;
  DAGLineStorage(LineId id, VertexId parent, VertexId child, Clade clade)
      : id_{id}, parent_{parent}, child_{child}, clade_{clade} {}

  template <typename T>
  DAGLineStorage& operator=(DAGLineView<T> other) {
    *this = other.line_;
    return *this;
  }

 private:
  // :: is workaround for GCC bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=52625
  template <typename>
  friend class ::DAGLine;
  DAGLineStorage& storage() { return *this; }
  const DAGLineStorage& storage() const { return *this; }

  LineId id_ = NoId;
  VertexId parent_ = NoId;
  VertexId child_ = NoId;
  Clade clade_ = Clade::Unspecified;
};

template <typename T>
class DAGLineView : public DAGLine<DAGLineView<T>> {
 public:
  DAGLineView(T& line) : line_{line} {
    Assert(line.GetClade() != Clade::Unspecified, "Uninitialized edge");
  }

  template <size_t I>
  auto get() const {
    static_assert(I < 2, "Index out of bounds");
    if constexpr (I == 0) return line_.GetVertexIds();
    if constexpr (I == 1) return line_.GetId();
  }

 private:
  template <typename>
  friend class DAGLine;
  friend class DAGLineStorage;
  DAGLineStorage& storage() { return line_; }
  const DAGLineStorage& storage() const { return line_; }

  T& line_;
};
using LineView = DAGLineView<DAGLineStorage>;
using ConstLineView = DAGLineView<const DAGLineStorage>;

namespace std {
template <>
struct tuple_size<::LineView> : public integral_constant<size_t, 2> {};
template <>
struct tuple_element<0, ::LineView> {
  using type = pair<::VertexId, ::VertexId>;
};
template <>
struct tuple_element<1, ::LineView> {
  using type = ::LineId;
};

template <>
struct tuple_size<::ConstLineView> : public integral_constant<size_t, 2> {};
template <>
struct tuple_element<0, ::ConstLineView> {
  using type = const pair<::VertexId, ::VertexId>;
};
template <>
struct tuple_element<1, ::ConstLineView> {
  using type = const ::LineId;
};
}  // namespace std

template <typename T>
class GenericNeighborsView {
  using iterator_type =
      std::conditional_t<std::is_const_v<T>, std::map<VertexId, LineId>::const_iterator,
                         std::map<VertexId, LineId>::iterator>;
  using map_type =
      std::conditional_t<std::is_const_v<T>, const std::map<VertexId, LineId>,
                         std::map<VertexId, LineId>>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, VertexId> {
   public:
    Iterator(const iterator_type& i, map_type& map) : iter_{i}, map_{map} {}

    Iterator operator++() {
      ++iter_;
      return *this;
    }
    bool operator!=(const Iterator& other) const { return iter_ != other.iter_; }
    bool operator==(const Iterator& other) const { return iter_ == other.iter_; }

    VertexId operator*() { return iter_->first; }

    VertexId GetNodeId() const { return iter_->first; }
    LineId GetEdge() const { return iter_->second; }

   private:
    iterator_type iter_;
    map_type& map_;
  };

  GenericNeighborsView(map_type& neighbors) : neighbors_{neighbors} {}
  template <typename U>
  GenericNeighborsView(const GenericNeighborsView<U>& other)
      : neighbors_{other.neighbors_} {}

  Iterator begin() const { return {neighbors_.begin(), neighbors_}; }
  Iterator end() const { return {neighbors_.end(), neighbors_}; }
  size_t size() const { return neighbors_.size(); }
  bool empty() const { return neighbors_.empty(); }

  void RemapIds(const SizeVector& reindexer) {
    for (auto [vertex_id, line_id] : neighbors_) {
      std::ignore = line_id;
      Assert(vertex_id < reindexer.size(),
             "Neighbors cannot contain an id out of bounds of the reindexer in "
             "GenericNeighborsView::RemapIds.");
    }
    T remapped{};
    for (auto [vertex_id, line_id] : neighbors_) {
      remapped[reindexer.at(vertex_id)] = line_id;
    }
    neighbors_ = remapped;
  }

  operator SizeVector() const {
    SizeVector result;
    for (auto i : neighbors_) result.push_back(i.first);
    return result;
  }

  void SetNeighbors(const T& neighbors) { neighbors_ = neighbors; }

 private:
  friend class Iterator;
  template <typename>
  friend class ::GenericNeighborsView;
  T& neighbors_;
};
using NeighborsView = GenericNeighborsView<std::map<VertexId, LineId>>;
using ConstNeighborsView = GenericNeighborsView<const std::map<VertexId, LineId>>;

class DAGVertex;
template <typename>
class GenericSubsplitDAGNode;
using MutableSubsplitDAGNode = GenericSubsplitDAGNode<DAGVertex>;
using SubsplitDAGNode = GenericSubsplitDAGNode<const DAGVertex>;

class DAGVertex {
 public:
  DAGVertex() = default;
  DAGVertex(const DAGVertex&) = default;
  inline DAGVertex(SubsplitDAGNode node);
  inline DAGVertex(MutableSubsplitDAGNode node);
  DAGVertex(VertexId id, Bitset subsplit) : id_{id}, subsplit_{std::move(subsplit)} {}

  VertexId GetId() const { return id_; }
  const Bitset& GetSubsplit() const { return subsplit_; }

  NeighborsView GetNeighbors(Direction direction, Clade clade) {
    return neighbors_.at({direction, clade});
  }

  ConstNeighborsView GetNeighbors(Direction direction, Clade clade) const {
    return neighbors_.at({direction, clade});
  }

  std::optional<std::tuple<LineId, Direction, Clade>> FindNeighbor(
      VertexId neighbor) const {
    for (auto i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) return {{j->second, i.first.first, i.first.second}};
    }
    return {};
  }

  DAGVertex& SetId(VertexId id) {
    id_ = id;
    return *this;
  }
  DAGVertex& SetSubsplit(Bitset subsplit) {
    subsplit_ = std::move(subsplit);
    return *this;
  }
  DAGVertex& AddNeighbor(Direction direction, Clade clade, VertexId neighbor,
                         LineId line) {
    neighbors_.at({direction, clade}).insert({neighbor, line});
    return *this;
  }
  void RemoveNeighbor(Direction direction, Clade clade, VertexId neighbor) {
    neighbors_.at({direction, clade}).erase(neighbor);
  }

  void SetLineId(VertexId neighbor, LineId line) {
    for (auto& i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) {
        i.second.insert_or_assign(j, neighbor, line);
        return;
      }
    }
    Failwith("Neighbor not found");
  }

  void ClearNeighbors() {
    neighbors_ = {
        {{Direction::Rootward, Clade::Left}, {}},
        {{Direction::Rootward, Clade::Right}, {}},
        {{Direction::Leafward, Clade::Left}, {}},
        {{Direction::Leafward, Clade::Right}, {}},
    };
  }

 private:
  VertexId id_ = NoId;
  Bitset subsplit_ = Bitset{{}};
  std::map<std::pair<Direction, Clade>, std::map<VertexId, LineId>> neighbors_ = {
      {{Direction::Rootward, Clade::Left}, {}},
      {{Direction::Rootward, Clade::Right}, {}},
      {{Direction::Leafward, Clade::Left}, {}},
      {{Direction::Leafward, Clade::Right}, {}},
  };
};

template <typename T>
class GenericLinesView {
  using view_type = std::conditional_t<std::is_const_v<T>, ConstLineView, LineView>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, view_type> {
   public:
    Iterator(const GenericLinesView& view, size_t index) : view_{view}, index_{index} {}

    Iterator operator++() {
      ++index_;
      return *this;
    }
    bool operator!=(const Iterator& other) const { return index_ != other.index_; }
    view_type operator*() const { return view_.storage_.lines_[index_]; }

   private:
    const GenericLinesView& view_;
    size_t index_;
  };

  GenericLinesView(T& storage) : storage_{storage} {}

  Iterator begin() const { return {*this, 0}; }
  Iterator end() const { return {*this, storage_.lines_.size()}; }
  size_t size() const { return storage_.lines_.size(); };
  view_type operator[](size_t i) const { return storage_.lines_[i]; };

 private:
  friend class Iterator;
  T& storage_;
};
class SubsplitDAGStorage;
using LinesView = GenericLinesView<SubsplitDAGStorage>;
using ConstLinesView = GenericLinesView<const SubsplitDAGStorage>;

template <typename T>
class GenericVerticesView {
  using view_type =
      std::conditional_t<std::is_const_v<T>, SubsplitDAGNode, MutableSubsplitDAGNode>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, view_type> {
   public:
    Iterator(T& storage, size_t index) : storage_{storage}, index_{index} {}

    Iterator& operator++() {
      ++index_;
      return *this;
    }
    Iterator operator++(int) { return {storage_, index_++}; }
    Iterator operator+(size_t i) { return {storage_, index_ + i}; }
    Iterator operator-(size_t i) { return {storage_, index_ - i}; }
    bool operator<(const Iterator& other) { return index_ < other.index_; }
    bool operator!=(const Iterator& other) const { return index_ != other.index_; }
    view_type operator*() const;

   private:
    T& storage_;
    size_t index_;
  };

  GenericVerticesView(T& storage) : storage_{storage} {}

  Iterator begin() const { return {storage_, 0}; }
  Iterator end() const { return {storage_, storage_.vertices_.size()}; }
  Iterator cbegin() const { return {storage_, 0}; }
  Iterator cend() const { return {storage_, storage_.vertices_.size()}; }
  size_t size() const { return storage_.vertices_.size(); };
  view_type operator[](size_t i) const;
  view_type at(size_t i) const;

 private:
  friend class Iterator;
  T& storage_;
};
using VerticesView = GenericVerticesView<SubsplitDAGStorage>;
using ConstVerticesView = GenericVerticesView<const SubsplitDAGStorage>;

// A vector that can optionally be prepended with host data for Graft purposes.
// It has no ownership over the host data, so the lifetime of the host data object
// should be greater than the given HostableVector instance.
template <typename T>
class HostableVector {
 public:
  explicit HostableVector(HostableVector<T>* host = nullptr)
      : host_{host ? &host->data_ : nullptr} {}

  T& at(size_t i) {
    if (!host_) {
      return data_.at(i);
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_.at(i - host_->size());
  }

  const T& at(size_t i) const {
    if (!host_) {
      return data_.at(i);
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_.at(i - host_->size());
  }

  T& operator[](size_t i) {
    if (!host_) {
      return data_[i];
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_[i - host_->size()];
  }

  const T& operator[](size_t i) const {
    if (!host_) {
      return data_[i];
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_[i - host_->size()];
  }

  size_t size() const {
    if (!host_) {
      return data_.size();
    }
    return host_->size() + data_.size();
  }

  void resize(size_t new_size) {
    if (!host_) {
      data_.resize(new_size);
    } else {
      data_.resize(new_size - host_->size());
    }
  }

  HostableVector& operator=(const std::vector<T>& data) {
    data_ = data;
    return *this;
  }

  bool HaveHost() const { return host_; }

  size_t HostSize() const {
    if (!host_) {
      return 0;
    }
    return host_->size();
  }

  void ResetHost(HostableVector<T>* host) {
    host_ = host ? &host->data_ : nullptr;
    data_ = {};
  }

 private:
  std::vector<T> data_;
  std::vector<T>* host_;
};

// Tag dispatching type to avoid confusion with copy constructor.
// See https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Tag_Dispatching
struct HostDispatchTag {};

class SubsplitDAGStorage {
 public:
  SubsplitDAGStorage() = default;
  SubsplitDAGStorage(const SubsplitDAGStorage&) = default;
  SubsplitDAGStorage(SubsplitDAGStorage& host, HostDispatchTag)
      : lines_{&host.lines_}, vertices_{&host.vertices_} {}

  ConstVerticesView GetVertices() const { return *this; }
  VerticesView GetVertices() { return *this; }

  ConstLinesView GetLines() const { return *this; }
  LinesView GetLines() { return *this; }

  std::optional<ConstLineView> GetLine(VertexId parent, VertexId child) const {
    if (parent >= vertices_.size() || child >= vertices_.size()) return {};
    auto result = vertices_.at(parent).FindNeighbor(child);
    if (!result.has_value()) return {};
    return lines_.at(std::get<0>(result.value()));
  }

  std::optional<ConstLineView> GetLine(LineId id) const {
    if (id >= lines_.size()) {
      return std::nullopt;
    }
    auto& line = lines_[id];
    if (line.GetId() == NoId) {
      return std::nullopt;
    }
    return line;
  }

  const DAGVertex& GetVertex(VertexId id) const { return vertices_.at(id); }

  bool ContainsVertex(VertexId id) const {
    if (id >= vertices_.size()) return false;
    return vertices_[id].GetId() != NoId;
  }

  std::optional<std::reference_wrapper<DAGVertex>> FindVertex(const Bitset& subsplit) {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      if (vertices_[i].GetSubsplit() == subsplit) return vertices_[i];
    }
    return std::nullopt;
  }

  std::optional<std::reference_wrapper<const DAGVertex>> FindVertex(
      const Bitset& subsplit) const {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      if (vertices_[i].GetSubsplit() == subsplit) return vertices_[i];
    }
    return std::nullopt;
  }

  DAGLineStorage& AddLine(const DAGLineStorage& newLine) {
    if (newLine.GetId() == NoId) Failwith("Set line id before inserting to storage");
    if (newLine.GetClade() == Clade::Unspecified)
      Failwith("Set clade before inserting to storage");
    auto& line = GetOrInsert(lines_, newLine.GetId());
    line = newLine;
    return line;
  }

  DAGVertex& AddVertex(const DAGVertex& newVertex) {
    if (newVertex.GetId() == NoId)
      Failwith("Set vertex id before inserting to storage");
    auto& vertex = GetOrInsert(vertices_, newVertex.GetId());
    vertex = newVertex;
    return vertex;
  }

  void ReindexLine(LineId line, VertexId parent, VertexId child) {
    lines_.at(line).SetParent(parent).SetChild(child);
  }

  void SetLines(const std::vector<DAGLineStorage>& lines) { lines_ = lines; }

  void SetVertices(const std::vector<DAGVertex>& vertices) { vertices_ = vertices; }

  bool HaveHost() const { return lines_.HaveHost(); }

  size_t HostLinesCount() const { return lines_.HostSize(); }

  size_t HostVerticesCount() const { return vertices_.HostSize(); }

  void ResetHost(SubsplitDAGStorage& host) {
    for (size_t i = lines_.HostSize(); i < lines_.size(); ++i) {
      auto& line = lines_[i];
      if (line.GetParent() >= vertices_.HostSize()) {
        if (line.GetChild() < vertices_.HostSize()) {
          vertices_[line.GetChild()].RemoveNeighbor(Direction::Rootward,
                                                    line.GetClade(), line.GetParent());
        }
      }
      if (line.GetChild() >= vertices_.HostSize()) {
        if (line.GetParent() < vertices_.HostSize()) {
          vertices_[line.GetParent()].RemoveNeighbor(Direction::Leafward,
                                                     line.GetClade(), line.GetChild());
        }
      }
    }
    lines_.ResetHost(&host.lines_);
    vertices_.ResetHost(&host.vertices_);
  }

  void ConnectVertices(LineId line_id) {
    auto& line = lines_.at(line_id);
    auto& parent = vertices_.at(line.GetParent());
    auto& child = vertices_.at(line.GetChild());
    parent.AddNeighbor(Direction::Leafward, line.GetClade(), child.GetId(), line_id);
    child.AddNeighbor(Direction::Rootward, line.GetClade(), parent.GetId(), line_id);
  }

  void ConnectAllVertices() {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      vertices_[i].ClearNeighbors();
    }
    for (size_t i = 0; i < lines_.size(); ++i) {
      if (lines_[i].GetId() == NoId) {
        continue;
      }
      ConnectVertices(lines_[i].GetId());
    }
  }

  std::optional<std::reference_wrapper<const DAGVertex>> FindRoot() const {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      auto& vertex = vertices_[i];
      if (vertex.GetId() == NoId) {
        continue;
      }
      if (vertex.GetNeighbors(Direction::Rootward, Clade::Left).empty() &&
          vertex.GetNeighbors(Direction::Rootward, Clade::Right).empty()) {
        return vertex;
      }
    }
    return std::nullopt;
  }

 private:
  template <typename>
  friend class GenericLinesView;
  template <typename>
  friend class GenericVerticesView;

  template <typename T, typename Id>
  static T& GetOrInsert(HostableVector<T>& data, Id id) {
    if (id >= data.size()) data.resize(id + 1);
    return data[id];
  }

  HostableVector<DAGLineStorage> lines_;
  HostableVector<DAGVertex> vertices_;
};

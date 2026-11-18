// Force simulation for Dorling and Demers cartograms.
//
// A direct port of src/cartogram/force-simulation.js. Pairwise repulsion is
// applied sequentially (each pair sees the moves of the previous pairs), so
// the loops cannot be vectorised in R without changing the results.

#include <cpp11.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace cpp11;

namespace {

const double PI = 3.14159265358979323846;

// Linear congruential generator matching createSeededRandom() in JS.
struct Lcg {
  uint32_t s;
  explicit Lcg(double seed) {
    double t = std::fmod(std::trunc(seed), 4294967296.0);
    if (t < 0) t += 4294967296.0;
    s = static_cast<uint32_t>(t);
    if (s == 0) s = 1;
  }
  double next() {
    s = s * 1664525u + 1013904223u;
    return s / 4294967296.0;
  }
};

struct Nodes {
  std::vector<double> x, y, ox, oy, r, w, h;
  size_t n() const { return x.size(); }
};

void unit_vector(Lcg& rng, double& ux, double& uy) {
  double angle = (rng.next() * 2 - 1) * PI;
  ux = std::cos(angle);
  uy = std::sin(angle);
}

void circle_repulsion(Nodes& nd, double strength, Lcg& rng) {
  size_t n = nd.n();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double dx = nd.x[j] - nd.x[i];
      double dy = nd.y[j] - nd.y[i];
      double dist = std::sqrt(dx * dx + dy * dy);
      double min_dist = nd.r[i] + nd.r[j];
      if (dist < min_dist && dist > 0) {
        double overlap = min_dist - dist;
        double nx = dx / dist, ny = dy / dist;
        double ma = nd.r[i] * nd.r[i], mb = nd.r[j] * nd.r[j];
        double total = ma + mb;
        double move_a = (overlap * mb / total) * strength;
        double move_b = (overlap * ma / total) * strength;
        nd.x[i] -= nx * move_a;
        nd.y[i] -= ny * move_a;
        nd.x[j] += nx * move_b;
        nd.y[j] += ny * move_b;
      } else if (dist == 0 && min_dist > 0) {
        double ux, uy;
        unit_vector(rng, ux, uy);
        double push = (min_dist * 0.5) * strength;
        nd.x[i] -= ux * push;
        nd.y[i] -= uy * push;
        nd.x[j] += ux * push;
        nd.y[j] += uy * push;
      }
    }
  }
}

void rect_repulsion(Nodes& nd, double strength, Lcg& rng) {
  size_t n = nd.n();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double dx = nd.x[j] - nd.x[i];
      double dy = nd.y[j] - nd.y[i];
      double ox = (nd.w[i] / 2 + nd.w[j] / 2) - std::fabs(dx);
      double oy = (nd.h[i] / 2 + nd.h[j] / 2) - std::fabs(dy);
      if (ox > 0 && oy > 0) {
        double ma = nd.w[i] * nd.h[i], mb = nd.w[j] * nd.h[j];
        double total = ma + mb;
        if (ox < oy) {
          double sign = dx >= 0 ? 1 : -1;
          nd.x[i] -= sign * (ox * mb / total) * strength;
          nd.x[j] += sign * (ox * ma / total) * strength;
        } else {
          double sign = dy >= 0 ? 1 : -1;
          nd.y[i] -= sign * (oy * mb / total) * strength;
          nd.y[j] += sign * (oy * ma / total) * strength;
        }
      } else if (dx == 0 && dy == 0) {
        double ux, uy;
        unit_vector(rng, ux, uy);
        nd.x[i] -= ux * 0.001 * strength;
        nd.y[i] -= uy * 0.001 * strength;
        nd.x[j] += ux * 0.001 * strength;
        nd.y[j] += uy * 0.001 * strength;
      }
    }
  }
}

bool overlaps(const Nodes& nd, bool rect) {
  size_t n = nd.n();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (!rect) {
        double dx = nd.x[j] - nd.x[i];
        double dy = nd.y[j] - nd.y[i];
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < (nd.r[i] + nd.r[j]) * 0.99) return true;
      } else {
        double ox = (nd.w[i] / 2 + nd.w[j] / 2) - std::fabs(nd.x[j] - nd.x[i]);
        double oy = (nd.h[i] / 2 + nd.h[j] / 2) - std::fabs(nd.y[j] - nd.y[i]);
        if (ox > 0.01 && oy > 0.01) return true;
      }
    }
  }
  return false;
}

void repulsion(Nodes& nd, bool rect, double strength, Lcg& rng) {
  if (rect) {
    rect_repulsion(nd, strength, rng);
  } else {
    circle_repulsion(nd, strength, rng);
  }
}

Nodes make_nodes(doubles x, doubles y, doubles ox, doubles oy, doubles r,
                 doubles w, doubles h) {
  Nodes nd;
  nd.x.assign(x.begin(), x.end());
  nd.y.assign(y.begin(), y.end());
  nd.ox.assign(ox.begin(), ox.end());
  nd.oy.assign(oy.begin(), oy.end());
  nd.r.assign(r.begin(), r.end());
  nd.w.assign(w.begin(), w.end());
  nd.h.assign(h.begin(), h.end());
  return nd;
}

}  // namespace

[[cpp11::register]]
list force_simulation_cpp(doubles x, doubles y, doubles ox, doubles oy,
                          doubles r, doubles w, doubles h, bool rect,
                          integers edge_from, integers edge_to,
                          int iterations, double repulsion_strength,
                          double anchor_strength, double adjacency_strength,
                          double cooling, double min_repulsion, double seed,
                          double initial_jitter) {
  Nodes nd = make_nodes(x, y, ox, oy, r, w, h);
  size_t n = nd.n();
  Lcg rng(seed);
  double current = repulsion_strength;

  if (initial_jitter > 0) {
    for (size_t i = 0; i < n; i++) {
      double jx = (rng.next() * 2 - 1) * initial_jitter;
      double jy = (rng.next() * 2 - 1) * initial_jitter;
      nd.x[i] += jx;
      nd.y[i] += jy;
    }
  }

  R_xlen_t n_edges = edge_from.size();
  for (int it = 0; it < iterations; it++) {
    double progress = static_cast<double>(it) / iterations;
    double anchor = anchor_strength * progress * progress;

    int passes = overlaps(nd, rect) ? 3 : 1;
    for (int p = 0; p < passes; p++) {
      repulsion(nd, rect, current, rng);
    }

    for (R_xlen_t e = 0; e < n_edges; e++) {
      int a = edge_from[e], b = edge_to[e];
      if (a < 0 || b < 0 || a >= (int)n || b >= (int)n) continue;
      double dx = nd.x[b] - nd.x[a];
      double dy = nd.y[b] - nd.y[a];
      double dist = std::sqrt(dx * dx + dy * dy);
      double target = nd.r[a] + nd.r[b];
      if (dist > target) {
        double pull = (dist - target) * adjacency_strength;
        double nx = dx / dist, ny = dy / dist;
        nd.x[a] += nx * pull * 0.5;
        nd.y[a] += ny * pull * 0.5;
        nd.x[b] -= nx * pull * 0.5;
        nd.y[b] -= ny * pull * 0.5;
      }
    }

    for (size_t i = 0; i < n; i++) {
      nd.x[i] += (nd.ox[i] - nd.x[i]) * anchor;
      nd.y[i] += (nd.oy[i] - nd.y[i]) * anchor;
    }

    current = std::max(current * cooling, min_repulsion);
  }

  // Final clean-up: pure repulsion without anchoring.
  double cleanup = std::min(500.0, iterations / 2.0);
  for (int i = 0; i < cleanup && overlaps(nd, rect); i++) {
    repulsion(nd, rect, 1.0, rng);
  }

  writable::doubles out_x(nd.x.begin(), nd.x.end());
  writable::doubles out_y(nd.y.begin(), nd.y.end());
  return writable::list({"x"_nm = out_x, "y"_nm = out_y});
}

[[cpp11::register]]
bool has_overlaps_cpp(doubles x, doubles y, doubles r, doubles w, doubles h,
                      bool rect) {
  Nodes nd = make_nodes(x, y, x, y, r, w, h);
  return overlaps(nd, rect);
}

[[cpp11::register]]
list resolve_rect_overlaps_cpp(doubles x, doubles y, doubles ox, doubles oy,
                               doubles w, doubles h, int max_iterations,
                               double epsilon, double gap, double anchor) {
  Nodes nd = make_nodes(x, y, ox, oy, w, w, h);
  size_t n = nd.n();
  bool before = overlaps(nd, true);
  int moved_pairs = 0;
  int iter = 0;

  if (before) {
    for (iter = 0; iter < max_iterations; iter++) {
      int moved = 0;
      for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
          double dx = nd.x[j] - nd.x[i];
          double dy = nd.y[j] - nd.y[i];
          double ox_ = (nd.w[i] / 2 + nd.w[j] / 2) - std::fabs(dx);
          double oy_ = (nd.h[i] / 2 + nd.h[j] / 2) - std::fabs(dy);
          if (ox_ > epsilon && oy_ > epsilon) {
            double area_a = std::max(epsilon, nd.w[i] * nd.h[i]);
            double area_b = std::max(epsilon, nd.w[j] * nd.h[j]);
            double total = area_a + area_b;
            double push_a = area_b / total, push_b = area_a / total;
            if (ox_ <= oy_) {
              double sign = dx >= 0 ? 1 : -1;
              double delta = ox_ + gap;
              nd.x[i] -= sign * delta * push_a;
              nd.x[j] += sign * delta * push_b;
            } else {
              double sign = dy >= 0 ? 1 : -1;
              double delta = oy_ + gap;
              nd.y[i] -= sign * delta * push_a;
              nd.y[j] += sign * delta * push_b;
            }
            moved_pairs++;
            moved++;
          }
        }
      }
      if (anchor > 0) {
        for (size_t k = 0; k < n; k++) {
          nd.x[k] += (nd.ox[k] - nd.x[k]) * anchor;
          nd.y[k] += (nd.oy[k] - nd.y[k]) * anchor;
        }
      }
      if (moved == 0 || !overlaps(nd, true)) break;
    }
    iter += 1;
  }

  writable::doubles out_x(nd.x.begin(), nd.x.end());
  writable::doubles out_y(nd.y.begin(), nd.y.end());
  return writable::list({"x"_nm = out_x, "y"_nm = out_y,
                         "iterations"_nm = writable::integers({iter}),
                         "moved_pairs"_nm = writable::integers({moved_pairs}),
                         "overlaps_before"_nm = writable::logicals({r_bool(before)}),
                         "overlaps_after"_nm =
                             writable::logicals({r_bool(overlaps(nd, true))})});
}

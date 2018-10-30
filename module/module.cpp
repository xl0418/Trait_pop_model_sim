#include <tbb/tbb.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include "utl/LambertW.h"
#include "rndutils.hpp"


namespace py = pybind11;
using namespace py::literals;


namespace detail {


  inline double ztp_lambda_from_untruncated_mean(double mu)
  {
    if (mu > 38.0) return mu;
    return mu + utl::LambertW<0>(-mu * std::exp(-mu));
  }


  template <typename RENG>
  inline int rztpois_small(double lambda, RENG&& reng)
  {
    if (lambda < 2.22e-16) {
      return 1;
    }
    int k = 1;
    double t = (-std::exp(-lambda) / std::expm1(-lambda)) * lambda;
    auto s = t;
    const auto u = rndutils::uniform01(reng);
    while (s < u) {
      k += 1;
      t *= lambda / k;
      s += t;
    }
    return k;
  }


  template <typename RENG>
  inline int rztpois(double lambda, RENG&& reng)
  {
    if (lambda > 140.0) {
      auto poi = std::poisson_distribution<int>(lambda);
      for (;;) {
        const auto y = poi(reng);
        if (y > 0) return y;
      }
    }
    return rztpois_small(lambda, reng);
  }


  template <typename RENG>
  inline int rztbinom50_small(int n, RENG&& reng)
  {
    if (n == 0) {
      return 0;
    }
    const auto binom50 = rndutils::binomial50_small_distribution<>(n);
    for (;;) {
      const auto y = binom50(reng);
      if (y > 0) return y;
    }
  }


  template <typename RENG>
  inline int rztbinom50(int n, RENG&& reng)
  {
    if (n < 10000) {
      return rztbinom50_small(n, reng);
    }
    auto binom50 = std::binomial_distribution<>(n, 0.5);
    for (;;) {
      const auto y = binom50(reng);
      if (y > 0) return y;
    }
  }


  // double truncated binomial50
  template <typename RENG>
  inline int rtbinom50(int n, RENG&& reng)
  {
    if (n == 1) {
      return 0;
    }
    for (;;) {
      auto y = rztbinom50(n, reng);
      if (y < n) return y;
    }
  }


  // double truncated binomial50
  template <typename RENG>
  inline int ___rtbinom50(int n, RENG&& reng)
  {
    return std::min(rztbinom50(n, reng), n - 1);
  }

}


namespace impl
{

  template <typename T>
  using py_array = py::array_t<T, py::array::c_style | py::array::forcecast>;


  template <typename T, typename I>
  py_array<T> make_py_array(I i)
  {
    return py::array_t<T>(static_cast<ssize_t>(i));
  }


  template <typename T, typename I>
  py_array<T> make_py_array(I i, I j)
  {
    return py::array_t<T>(typename py_array<T>::ShapeContainer{ static_cast<ssize_t>(i), static_cast<ssize_t>(j) });
  }


  struct sim_res
  {
    ssize_t valid_time;
    std::vector<int> N;
    std::vector<double> Z;
    std::vector<double> V;
    const char* msg = "exception";
  };
  using sim_res_vector = std::vector<sim_res, tbb::cache_aligned_allocator<sim_res>>;


  py::object to_py(sim_res&& C)
  {
    const auto n = static_cast<ssize_t>(C.N.size());
    // NaN-ify
    for (auto i = 0; i < n; ++i)
    {
      if (C.N[i] == 0) C.Z[i] = C.V[i] = Py_NAN;
    }
    return py::dict{
      "sim_time"_a = C.valid_time,
      "N"_a = py_array<int>(n, C.N.data()),
      "Z"_a = py_array<double>(n, C.Z.data()),
      "V"_a = py_array<double>(n, C.V.data()),
      "msg"_a = py::str(C.msg)
    };
  }


  py::object to_py(sim_res_vector& R)
  {
    auto py_N = make_py_array<int>(R.size(), R[0].N.size());
    auto py_Z = make_py_array<double>(R.size(), R[0].Z.size());
    auto py_V = make_py_array<double>(R.size(), R[0].V.size());
    auto py_t = make_py_array<ssize_t>(R.size());
    auto py_msg = std::vector<std::string>(R.size());
    auto pN = py_N.mutable_unchecked();
    auto pZ = py_Z.mutable_unchecked();
    auto pV = py_V.mutable_unchecked();
    auto pt = py_t.mutable_unchecked();
    const auto n = static_cast<ssize_t>(R[0].N.size());
    for (size_t j = 0; j < R.size(); ++j)
    {
      pt(j) = R[j].valid_time;
      if (!R[j].N.empty())
      {
        for (auto i = 0; i < n; ++i)
        {
          if (R[j].N[i] == 0) R[j].Z[i] = R[j].V[i] = Py_NAN;
          pN(j, i) = R[j].N[i];
          pZ(j, i) = R[j].Z[i];
          pV(j, i) = R[j].V[i];
        }
      }
      else
      {
        for (auto i = 0; i < n; ++i)
        {
          pN(j, i) = 0; 
          pZ(j, i) = pV(j, i) = Py_NAN;
        }
      }
      py_msg[j] = R[j].msg;
    }
    return py::dict{
      "sim_time"_a = py_t,
      "N"_a = std::move(py_N),
      "Z"_a = std::move(py_Z),
      "V"_a = std::move(py_V),
      "msg"_a = std::move(py_msg)
    };
  }


  enum sim_param : int
  {
    gamma = 0,
    a = 1,
    K = 2,
    nu = 3,
    r = 4,
    theta = 5,
    Vmax = 6,
    init_z = 7,
    init_n = 8,
    init_sigma = 9,
    break_on_mu = 10,
    max_param
  };


  template <typename EVENTS>
  sim_res sim_single(const ssize_t T, const ssize_t S, const EVENTS& events, const double* param)
  {
    const auto gamma = param[sim_param::gamma];
    const auto a = param[sim_param::a];
    const auto K = param[sim_param::K];
    const auto nu = param[sim_param::nu];
    const auto r = param[sim_param::r];
    const auto theta = param[sim_param::theta];
	  const auto Vmax = param[sim_param::Vmax];
    const bool break_on_mu = param[sim_param::break_on_mu] == 1.0;

    std::vector<double> Z(S, param[sim_param::init_z]), nZ(S);
    std::vector<int> N(S, 0), nN(S, 0);
    std::vector<double> V(S, 1.0 / S), nV(S, 0.0);
    auto reng = rndutils::make_random_engine_low_entropy<>();
    {
      auto init_n_dist = std::normal_distribution<double>(param[sim_param::init_n], param[sim_param::init_sigma]);
      N[0] = static_cast<int>(init_n_dist(reng));
      N[1] = static_cast<int>(init_n_dist(reng));
    }
    ssize_t sp = 2;     // existing species
    ssize_t node = 0;
    auto next_event = events.data(node, 0);
    ssize_t sim_time = 0;
    for (; sim_time < T; ++sim_time)
    {
      for (ssize_t i = 0; i < sp; ++i)
      {
        if (N[i] != 0.0)
        {
          double beta = 0.0;
          double sigma = 0.0;
          double sigmasqr = 0.0;
          const double zi = Z[i];
          for (ssize_t j = 0; j < sp; ++j)
          {
            const double zd = zi - Z[j];
            const double t1 = std::exp(-a * zd * zd) * N[j];
            const double t2 = 2.0 * a * zd;
            beta += t1;
            sigma += t2 * t1;
            sigmasqr += t2 * t2 * t1;
          }
          const auto var_trait = V[i] / (2.0 * N[i]);
          const auto dtz = theta - Z[i];
          nZ[i] = Z[i] + V[i] * (2.0 * gamma * dtz + 1.0 / K * sigma) + std::normal_distribution<double>(0.0, var_trait)(reng);
          const auto mu = N[i] * r * std::exp(-gamma * dtz * dtz + (1.0 - beta / K));
          const auto ztp_lambda = detail::ztp_lambda_from_untruncated_mean(mu);
          nN[i] = detail::rztpois(ztp_lambda, reng);
          nV[i] = V[i] / 2.0 + 2.0 * N[i] * nu * Vmax / (1.0 + 4.0 * N[i] * nu)
            + V[i] * V[i] * (
              -2.0 * gamma + 4.0 * gamma * gamma * dtz * dtz +
              (1.0 / K) * (2.0 * a * beta - sigmasqr) + 4.0 * gamma / K *
              dtz * sigma + sigma * sigma / (K * K)
              );
          // sanity checks
          const char* reason = nullptr;
          if (break_on_mu && (mu <= 1.0)) reason = "invalid mu";
          if (nN[i] < 1) reason = "inconsistent zero population";
          if ((nV[i] <= 0.0) || (nV[i] > 100000.0)) reason = "inconsistent variance";
          if (nullptr != reason)
		      {
  			    return { sim_time, std::move(nN), std::move(nZ), std::move(nV), reason };
		      }
        }
      }
      while ((sim_time + 1) == next_event[0])
      {
        const auto parent = next_event[1];
        const auto daughter = next_event[2];
        if (daughter == -1)
        { // extinction
          const auto ext = next_event[1];
          nN[ext] = N[ext] = 0;
        }
        else
        { // speciation
          auto split = detail::rtbinom50(nN[parent], reng);
          nN[daughter] = nN[parent] - split;
          nN[parent] = split;
          nV[parent] *= 0.5;
          nV[daughter] = nV[parent];
          nZ[daughter] = nZ[parent];
          ++sp;
          // sanity checks
          if (nN[parent] == 0 || nN[daughter] == 0) 
          {
            return { sim_time, std::move(nN), std::move(nZ), std::move(nV), "singleton split" };
          }
        }
        ++node;
        next_event = events.data(node, 0);
      }
      N.swap(nN);
      Z.swap(nZ);
      V.swap(nV);
    }
    return { sim_time, std::move(N), std::move(Z), std::move(V), "success" };
  }


  template <typename EVENTS, typename PARAM>
  py::object sim_multiple(ssize_t T, ssize_t S, const EVENTS& events, const PARAM& params)
  {
    const auto sets = params.shape(0);
    auto simres = sim_res_vector(sets);
	  tbb::parallel_for(py::ssize_t{0}, sets, [&](auto s)
    {
      simres[s] = sim_single(T, S, events, params.template unchecked<2>().data(s, 0));
	  });
    return to_py(simres);
  }


  py::object sim_dispatch(const py::object& treedata, py_array<double>& params)
  {
    const auto events = py::cast<py_array<ssize_t>>(treedata.attr("sim_events"));
    const auto T = py::cast<ssize_t>(treedata.attr("sim_evo_time"));
    const auto S = py::cast<ssize_t>(treedata.attr("total_species"));
    if ((events.ndim() == 2) && (events.shape(1) == 3))
    {
      if ((params.ndim() == 1) && (params.shape(0) == sim_param::max_param))
      {
        return impl::to_py(sim_single(T, S, events.unchecked<2>(), params.unchecked<1>().data(0)));
      }
      if ((params.ndim() == 2) && (params.shape(1) == sim_param::max_param))
      {
        return sim_multiple(T, S, events.unchecked<2>(), params);
      }
    }
    throw std::runtime_error("parameter don't match");
    return py::none{};
  }


  py_array<ssize_t> discrete_distribution(const py_array<double>& val, ssize_t num)
  {
    using ddist = rndutils::mutable_discrete_distribution<ssize_t, rndutils::all_zero_policy_uni>;
    auto dist = ddist(val.data(), val.data() + val.shape(0));
    auto res = make_py_array<ssize_t>(num);
    auto p = res.mutable_data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p)
    {
      *p = dist(reng);
    }
    return res;
  }


  py_array<double> cauchy_distribution(double a, double b, ssize_t num)
  {
    auto dist = std::cauchy_distribution<>(a, b);
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p)
    {
      *p = dist(reng);
    }
    return res;
  }


  py_array<double> ztp_lambda_from_untruncated_mean(const py::object& pymu)
  {
    const auto mu = py::cast<py_array<double>>(pymu);
    auto num = mu.size();
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto a = mu.data();
    for (ssize_t i = 0; i < num; ++i, ++p, ++a)
    {
      *p = detail::ztp_lambda_from_untruncated_mean(*a);
    }
    return res;
  }


  py_array<int> ztpoisson(const py_array<double>& pylambda)
  {
    auto num = pylambda.size();
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto a = pylambda.data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p, ++a)
    {
      *p = detail::rztpois(*a, reng);
    }
    return res;
  }


  py_array<int> tbinomial50(const py_array<int>& pyn)
  {
    auto num = pyn.size();
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto a = pyn.data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p, ++a)
    {
      *p = detail::rtbinom50(*a, reng);
    }
    return res;
  }


  py_array<int> ___tbinomial50(const py_array<int>& pyn)
  {
    auto num = pyn.size();
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto a = pyn.data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p, ++a)
    {
      *p = detail::___rtbinom50(*a, reng);
    }
    return res;
  }

}


PYBIND11_MODULE(dvtraitsim_cpp, m) {
  m.doc() = R"pbdoc(
        dvtraitsim_cpp plugin
        ---------------------

        .. currentmodule:: dvtraitsim_cpp

        .. autosummary::
           :toctree: _generate

           DVSim
           discrete_distribution
           cauchy_distribution
           ztp_lambda_from_untruncated_mean
           ztpoisson
           ztbinomial50
    )pbdoc";

  m.def("DVSim", &impl::sim_dispatch, "DVSim");
  m.def("discrete_distribution", &impl::discrete_distribution, "discrete_distribution");
  m.def("cauchy_distribution", &impl::cauchy_distribution, "a"_a = 0.0, "b"_a = 1.0, "num"_a, "cauchy_distribution");

  m.def("ztp_lambda_from_untruncated_mean", 
    [](double mu) {return detail::ztp_lambda_from_untruncated_mean(mu); }, 
    "mu"_a, "lambda from untruncated mean");

  m.def("ztp_lambda_from_untruncated_mean", 
    [](impl::py_array<double> mu) { return impl::ztp_lambda_from_untruncated_mean(mu); }, 
    "mu"_a, "lambda from untruncated mean");

  m.def("ztpoisson",
    [](double lambda) { return detail::rztpois(lambda, rndutils::make_random_engine_low_entropy<>()); },
    "lambda"_a, "zero truncated Poisson");

  m.def("ztpoisson",
    [](impl::py_array<double> lambda) { return impl::ztpoisson(lambda); },
    "lambda"_a, "zero truncated Poisson");

  m.def("split_binomial50",
    [](int n) { return detail::rtbinom50(n, rndutils::make_random_engine_low_entropy<>()); },
    "lambda"_a, "(zero, n) truncated Binomial_p=0.5");

  m.def("split_binomial50",
    [](impl::py_array<int> n) { return impl::tbinomial50(n); },
    "n"_a, "(zero, n) truncated Binomial_p=0.5");

  m.attr("all") = 0;
  m.attr("nodes") = 1;
  m.attr("tips") = 2;

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "0.0.5";
#endif
  m.attr("__author__") = "Hanno Hildenbrandt";
}


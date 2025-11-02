// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2014 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <qlat-utils/qassert.h>
#include <qlat-utils/env.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/show.h>
#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifdef USE_PAPI
#include <papi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#include <qlat-utils/compatible-omp.h>
#endif

#define TIMER(FNAME)                     \
  static std::string fname = FNAME;      \
  static qlat::Timer timer(fname, true); \
  qlat::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE(FNAME)             \
  static std::string fname = FNAME;      \
  static qlat::Timer timer(fname, true); \
  qlat::TimerCtrl timerctrl(timer, true);

#define TIMER_FLOPS(FNAME)                \
  static std::string fname = FNAME;       \
  static qlat::Timer timer(fname, false); \
  qlat::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE_FLOPS(FNAME)        \
  static std::string fname = FNAME;       \
  static qlat::Timer timer(fname, false); \
  qlat::TimerCtrl timerctrl(timer, true);

namespace qlat
{  //

inline RealD get_time()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((RealD)tp.tv_sec + (RealD)tp.tv_usec * 1e-6);
}

API inline RealD& get_actual_start_time()
// not affected by Timer::reset()
{
  static RealD time = get_time();
  return time;
}

inline RealD get_actual_total_time()
{
  return get_time() - get_actual_start_time();
}

API inline RealD& get_start_time()
// will be reset by Timer::reset()
{
  static RealD time = get_actual_start_time();
  return time;
}

inline RealD get_remaining_time()
{
  return get_time_limit() - get_actual_total_time();
}

inline RealD get_total_time() { return get_time() - get_start_time(); }

inline Int get_num_thread()
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline Int get_id_thread()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

API inline Int& get_id_node_in_shuffle_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static Int id_node_in_shuffle = 0;
  return id_node_in_shuffle;
}

inline Int get_id_node_in_shuffle()
{
  return get_id_node_in_shuffle_internal();
}

API inline Int& get_num_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static Int num_node = 1;
  return num_node;
}

API inline Int& get_id_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static Int id_node = 0;
  return id_node;
}

inline Int get_num_node() { return get_num_node_internal(); }

inline Int get_id_node() { return get_id_node_internal(); }

inline void display(const std::string& str)
{
  if (0 == get_id_thread()) {
    get_display_ptr()(str);
  }
}

inline void displayln(const std::string& str)
{
  if (0 == get_id_thread()) {
    DisplayPtr display = get_display_ptr();
    display(str);
    display("\n");
  }
}

inline void display(const Long minimum_verbose_level, const std::string& str)
{
  if (get_verbose_level() >= minimum_verbose_level) {
    display(str);
  }
}

inline void displayln(const Long minimum_verbose_level, const std::string& str)
{
  if (get_verbose_level() >= minimum_verbose_level) {
    displayln(str);
  }
}

inline void display_info(const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(str);
  }
}

inline void display_info(const Long minimum_verbose_level,
                         const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(minimum_verbose_level, str);
  }
}

inline void displayln_info(const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    displayln(str);
  }
}

inline void displayln_info(const Long minimum_verbose_level,
                           const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    displayln(minimum_verbose_level, str);
  }
}

inline Long get_total_flops()
{
  Long flops = 0;
#ifdef USE_PAPI
  const Int n_threads = omp_get_max_threads();
  Long flopses[n_threads];
  memset(flopses, 0, n_threads * sizeof(Long));
#pragma omp parallel
  {
    RealF rtime, ptime, mflops;
    Int i = omp_get_thread_num();
    PAPI_flops(&rtime, &ptime, &flopses[i], &mflops);
  }
  for (Int i = 0; i < n_threads; i++) {
    flops += flopses[i];
  }
#endif
  return flops;
}

template <class K>
bool has(const std::set<K>& m, const K& key)
{
  typename std::set<K>::const_iterator it = m.find(key);
  return it != m.end();
}

template <class K, class M>
bool has(const std::map<K, M>& m, const K& key)
{
  typename std::map<K, M>::const_iterator it = m.find(key);
  return it != m.end();
}

API inline void initialize_papi()
{
#ifdef USE_PAPI
  static bool initialized = false;
  if (initialized) {
    return;
  }
  displayln_info(0, "PAPI::initialize_papi Start.");
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_thread_init((uint64_t (*)(void))(omp_get_thread_num));
  initialized = true;
  displayln_info(0, "PAPI::initialize_papi Finish.");
#endif
}

struct API TimerInfo {
  std::string fname;
  RealD dtime;
  RealD accumulated_time;
  Long dflops;
  Long accumulated_flops;
  Int call_times;
  //
  TimerInfo() { init(); }
  //
  void init()
  {
    fname = "Unknown";
    reset();
  }
  //
  void reset();
  //
  void merge(const TimerInfo& x);
  //
  void show_start(const Int fname_len) const;
  //
  void show_stop(const Int fname_len) const;
  //
  void show_avg_always(const std::string& info, const Int fname_len) const;
  //
  void show_avg(const std::string& info, const Int fname_len) const;
};

struct API Timer {
  const char* cname;
  Long info_index;
  bool is_using_total_flops;
  Long is_running;
  RealD start_time;
  RealD stop_time;
  Long start_flops;
  Long stop_flops;
  Long flops;
  //
  API static std::map<std::string, Long>& get_timer_info_index_map()
  {
    static std::map<std::string, Long> timer_info_index_map;
    return timer_info_index_map;
  }
  //
  API static std::vector<TimerInfo>& get_timer_database()
  {
    static std::vector<TimerInfo> timer_database;
    return timer_database;
  }
  //
  API static std::vector<std::vector<TimerInfo>>& get_timer_database_history()
  {
    static std::vector<std::vector<TimerInfo>> timer_database_history;
    return timer_database_history;
  }
  //
  API static std::vector<RealD>& get_start_time_history()
  {
    static std::vector<RealD> history;
    return history;
  }
  //
  API static std::vector<Long>&
  get_max_call_times_for_always_show_info_history()
  {
    static std::vector<Long> history;
    return history;
  }
  //
  API static std::vector<Long>& get_timer_stack()
  {
    static std::vector<Long> stack;
    return stack;
  }
  //
  API static void reset(const Long max_call_times_for_always_show_info_ = -1);
  //
  API static void fork(const Long max_call_times_for_always_show_info_ = -1);
  //
  API static void merge();
  //
  API static RealD& minimum_autodisplay_interval()
  // qlat parameter
  {
    static RealD time =
        get_env_double_default("q_timer_mini_auto_display", 5.0 * 60.0);
    return time;
  }
  //
  API static RealD& minimum_duration_for_show_info()
  // qlat parameter
  {
    static RealD time = get_env_double_default("q_timer_mini_auto_show", 1.0);
    return time;
  }
  //
  API static RealD& minimum_duration_for_show_stop_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static RealD& minimum_duration_for_show_start_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static Long& max_call_times_for_always_show_info()
  // qlat parameter
  {
    static Long max_call_times =
      get_env_long_default("q_timer_max_always_show", 2);
    return max_call_times;
  }
  //
  API static Long& max_function_name_length_shown()
  // qlat parameter
  {
    static Long max_len = get_env_long_default("q_timer_max_func_name_len", 75);
    return max_len;
  }
  //
  Timer() { init(); }
  Timer(const std::string& fname_str)
  {
    init();
    init(fname_str);
  }
  Timer(const std::string& cname_str, const std::string& fname_str)
  {
    init();
    init(cname_str, fname_str);
  }
  Timer(const std::string& fname_str, const bool is_using_total_flops_)
  {
    init();
    init(fname_str);
    is_using_total_flops = is_using_total_flops_;
  }
  //
  void init();
  //
  void init(const std::string& fname_str);
  //
  void init(const std::string& cname_str, const std::string& fname_str);
  //
  void start(bool verbose = false);
  //
  void stop(bool verbose = false);
  //
  API static void test_timer_time_usage();
  //
  API static void display(const std::string& tag = "");
  //
  API static void autodisplay(const RealD time = get_time());
  //
  API static void display_stack_always();
  //
  API static void display_stack()
  // only display if id_node == 0 and thread_num == 0
  {
    if (0 == get_id_node() && 0 == get_id_thread()) {
      display_stack_always();
    }
  }
};

struct API TimerCtrl {
  Timer* ptimer;
  bool verbose;
  //
  TimerCtrl() { init(); }
  TimerCtrl(Timer& timer, bool verbose_ = false)
  {
    init();
    init(timer, verbose_);
  }
  //
  ~TimerCtrl()
  {
    if (NULL != ptimer) {
      ptimer->stop(verbose);
    }
  }
  //
  void init()
  {
    ptimer = NULL;
    verbose = false;
  }
  void init(Timer& timer, bool verbose_ = false)
  {
    if (get_id_thread() != 0) return;
    ptimer = &timer;
    verbose = verbose_;
    ptimer->start(verbose);
  }
};

///////////////////////////////////////////////////////////////////////

inline void Display(const char* cname, const char* fname, const char* format,
                    ...)
{
  va_list args;
  va_start(args, format);
  char* str;
  const Int ret = vasprintf(&str, format, args);
  Qassert(ret >= 0);
  display(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname,
                        const char* format, ...)
{
  Int rank = get_id_node();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  const Int ret = vasprintf(&str, format, args);
  Qassert(ret >= 0);
  display_info(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void* timer_malloc(size_t size)
{
  TIMER_FLOPS("timer_malloc");
  timer.flops += size;
  void* p = malloc(size);
  memset(p, 0, size);
  return p;
}

inline void timer_free(void* ptr)
{
  TIMER_FLOPS("timer_free");
  free(ptr);
}

inline void* tmalloc(size_t size) { return timer_malloc(size); }

inline void tfree(void* ptr) { timer_free(ptr); }

}  // namespace qlat

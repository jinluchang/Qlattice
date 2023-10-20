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

#include <qlat-utils/assert.h>
#include <qlat-utils/show.h>
#include <qlat-utils/env.h>
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

inline double get_time()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

API inline double& get_actual_start_time()
// not affected by Timer::reset()
{
  static double time = get_time();
  return time;
}

inline double get_actual_total_time()
{
  return get_time() - get_actual_start_time();
}

API inline double& get_start_time()
// will be reset by Timer::reset()
{
  static double time = get_actual_start_time();
  return time;
}

inline double get_remaining_time()
{
  return get_time_limit() - get_actual_total_time();
}

inline double get_total_time() { return get_time() - get_start_time(); }

API inline int& get_num_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int num_node = 1;
  return num_node;
}

inline int get_num_node() { return get_num_node_internal(); }

API inline int& get_id_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int id_node = 0;
  return id_node;
}

inline int get_id_node() { return get_id_node_internal(); }

inline int get_num_thread()
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline int get_id_thread()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

API inline int& get_id_node_in_shuffle_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int id_node_in_shuffle = 0;
  return id_node_in_shuffle;
}

inline int get_id_node_in_shuffle()
{
  return get_id_node_in_shuffle_internal();
}

inline void display_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(str, fp);
  }
}

inline void displayln_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    displayln(str, fp);
  }
}

inline void display(const long minimum_verbose_level, const std::string& str)
{
  if (get_verbose_level() >= minimum_verbose_level) {
    display(str);
  }
}

inline void displayln(const long minimum_verbose_level, const std::string& str)
{
  if (get_verbose_level() >= minimum_verbose_level) {
    displayln(str);
  }
}

inline void display_info(const long minimum_verbose_level,
                         const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(minimum_verbose_level, str);
  }
}

inline void displayln_info(const long minimum_verbose_level,
                           const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    displayln(minimum_verbose_level, str);
  }
}

inline long long get_total_flops()
{
  long long flops = 0;
#ifdef USE_PAPI
  const int n_threads = omp_get_max_threads();
  long long flopses[n_threads];
  std::memset(flopses, 0, n_threads * sizeof(long long));
#pragma omp parallel
  {
    float rtime, ptime, mflops;
    int i = omp_get_thread_num();
    PAPI_flops(&rtime, &ptime, &flopses[i], &mflops);
  }
  for (int i = 0; i < n_threads; i++) {
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
  PAPI_thread_init((unsigned long (*)(void))(omp_get_thread_num));
  initialized = true;
  displayln_info(0, "PAPI::initialize_papi Finish.");
#endif
}

struct API TimerInfo {
  std::string fname;
  double dtime;
  double accumulated_time;
  long long dflops;
  long long accumulated_flops;
  int call_times;
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
  void show_start(const int fname_len) const;
  //
  void show_stop(const int fname_len) const;
  //
  void show_avg_always(const std::string& info, const int fname_len) const;
  //
  void show_avg(const std::string& info, const int fname_len) const
  {
    if (0 == get_id_node() && 0 == get_id_thread()) {
      show_avg_always(info, fname_len);
    }
  }
};

struct API Timer {
  const char* cname;
  long info_index;
  bool is_using_total_flops;
  long is_running;
  double start_time;
  double stop_time;
  long long start_flops;
  long long stop_flops;
  long long flops;
  //
  API static std::map<std::string, long>& get_timer_info_index_map()
  {
    static std::map<std::string, long> timer_info_index_map;
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
  API static std::vector<double>& get_start_time_history()
  {
    static std::vector<double> history;
    return history;
  }
  //
  API static std::vector<long>&
  get_max_call_times_for_always_show_info_history()
  {
    static std::vector<long> history;
    return history;
  }
  //
  API static std::vector<long>& get_timer_stack()
  {
    static std::vector<long> stack;
    return stack;
  }
  //
  API static void reset(const long max_call_times_for_always_show_info_ = -1);
  //
  API static void fork(const long max_call_times_for_always_show_info_ = -1);
  //
  API static void merge();
  //
  API static double& minimum_autodisplay_interval()
  // qlat parameter
  {
    static double time =
        get_env_double_default("q_timer_mini_auto_display", 5.0 * 60.0);
    return time;
  }
  //
  API static double& minimum_duration_for_show_info()
  // qlat parameter
  {
    static double time = get_env_double_default("q_timer_mini_auto_show", 1.0);
    return time;
  }
  //
  API static double& minimum_duration_for_show_stop_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static double& minimum_duration_for_show_start_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static long& max_call_times_for_always_show_info()
  // qlat parameter
  {
    static long max_call_times =
        get_env_long_default("q_timer_max_always_show", 10);
    return max_call_times;
  }
  //
  API static long& max_function_name_length_shown()
  // qlat parameter
  {
    static long max_len = get_env_long_default("q_timer_max_func_name_len", 50);
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
  API static void autodisplay(const double time = get_time());
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
  qassert(vasprintf(&str, format, args) >= 0);
  display(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname,
                        const char* format, ...)
{
  int rank = get_id_node();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  qassert(vasprintf(&str, format, args) >= 0);
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

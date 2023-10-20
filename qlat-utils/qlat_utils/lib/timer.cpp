#include <qlat-utils/timer.h>

namespace qlat
{  //

void TimerInfo::reset()
{
  dtime = 0.0 / 0.0;
  accumulated_time = 0;
  dflops = 0;
  accumulated_flops = 0;
  call_times = 0;
}

void TimerInfo::merge(const TimerInfo& x)
{
  if (std::isnan(dtime)) {
    dtime = x.dtime;
    dflops = x.dflops;
  }
  accumulated_time += x.accumulated_time;
  accumulated_flops += x.accumulated_flops;
  call_times += x.call_times;
}

void TimerInfo::show_start(const int fname_len) const
{
  double total_time = get_total_time();
  std::string fnameCut;
  fnameCut.assign(fname, 0, fname_len);
  displayln_info(ssprintf(
      "Timer::start %s :%5.1f%% %8d calls %.3E,%.3E sec %8.3f Gflops (%.3E "
      "flops)",
      ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str()).c_str(),
      accumulated_time / total_time * 100, call_times, dtime,
      accumulated_time / (call_times - 1), dflops / dtime / 1.0E9,
      (double)dflops));
}

void TimerInfo::show_stop(const int fname_len) const
{
  double total_time = get_total_time();
  std::string fnameCut;
  fnameCut.assign(fname, 0, fname_len);
  displayln_info(ssprintf(
      "Timer::stop  %s :%5.1f%% %8d calls %.3E,%.3E sec %8.3f Gflops (%.3E "
      "flops)",
      ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str()).c_str(),
      accumulated_time / total_time * 100, call_times, dtime,
      accumulated_time / call_times, dflops / dtime / 1.0E9, (double)dflops));
}

void TimerInfo::show_avg_always(const std::string& info,
                                const int fname_len) const
{
  double total_time = get_total_time();
  std::string fnameCut;
  fnameCut.assign(fname, 0, fname_len);
  displayln(ssprintf(
      "Timer::%s %s :%7.3f%% %8d calls; %.2E,%.2E sec; %.2E,%.2E flops; "
      "%5.2f Gflops",
      info.c_str(),
      ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str()).c_str(),
      accumulated_time / total_time * 100, call_times,
      accumulated_time / call_times, accumulated_time,
      (double)accumulated_flops / (double)call_times, (double)accumulated_flops,
      accumulated_flops / accumulated_time / 1.0E9));
}

static bool compare_time_info_p(const TimerInfo* p1, const TimerInfo* p2)
{
  return p1->accumulated_time < p2->accumulated_time;
}

void Timer::init()
{
  cname = "Timer";
  is_using_total_flops = false;
  get_start_time();
  initialize_papi();
  info_index = -1;
  is_running = 0;
}

void Timer::init(const std::string& fname_str)
{
  std::map<std::string, long>& tiim = get_timer_info_index_map();
  std::vector<TimerInfo>& tdb = get_timer_database();
  if (has(tiim, fname_str)) {
    info_index = tiim[fname_str];
  } else {
    info_index = tdb.size();
    tiim[fname_str] = info_index;
    TimerInfo info;
    info.fname = fname_str;
    tdb.push_back(info);
  }
}

void Timer::init(const std::string& cname_str, const std::string& fname_str)
{
  std::string fname = "";
  fname += cname_str;
  fname += "::";
  fname += fname_str;
  init(fname);
}

void Timer::reset(const long max_call_times_for_always_show_info_)
// if max_call_times_for_always_show_info_ <= -1:
// then do not change the current value.
// else update the max_call_times_for_always_show_info.
{
  std::vector<TimerInfo>& tdb = get_timer_database();
  for (long i = 0; i < (long)tdb.size(); ++i) {
    tdb[i].reset();
  }
  get_start_time() = get_time();
  if (max_call_times_for_always_show_info_ >= 0) {
    max_call_times_for_always_show_info() =
        max_call_times_for_always_show_info_;
  }
  std::vector<std::vector<TimerInfo>>& tdb_history =
      get_timer_database_history();
  displayln_info(0,
                 ssprintf("Timer::reset(%ld): Reset all timers! (level = %ld)",
                          max_call_times_for_always_show_info(),
                          (long)tdb_history.size()));
}

void Timer::fork(const long max_call_times_for_always_show_info_)
// if max_call_times_for_always_show_info_ <= -1:
// then do not change the current value.
// else update the max_call_times_for_always_show_info.
{
  get_start_time_history().push_back(get_start_time());
  get_start_time() = get_time();
  get_max_call_times_for_always_show_info_history().push_back(
      max_call_times_for_always_show_info());
  if (max_call_times_for_always_show_info_ >= 0) {
    max_call_times_for_always_show_info() =
        max_call_times_for_always_show_info_;
  }
  std::vector<std::vector<TimerInfo>>& tdb_history =
      get_timer_database_history();
  std::vector<TimerInfo>& tdb = get_timer_database();
  tdb_history.push_back(tdb);
  for (long i = 0; i < (long)tdb.size(); ++i) {
    tdb[i].reset();
  }
  displayln_info(0, ssprintf("Timer::fork(%ld): Fork all timers! (level = %ld)",
                             max_call_times_for_always_show_info(),
                             (long)tdb_history.size()));
}

void Timer::merge()
// call merge only after fork
{
  qassert(get_start_time_history().size() >= 1);
  get_start_time() = get_start_time_history().back();
  get_start_time_history().pop_back();
  qassert(get_max_call_times_for_always_show_info_history().size() >= 1);
  max_call_times_for_always_show_info() =
      get_max_call_times_for_always_show_info_history().back();
  get_max_call_times_for_always_show_info_history().pop_back();
  std::vector<std::vector<TimerInfo>>& tdb_history =
      get_timer_database_history();
  std::vector<TimerInfo>& tdb = get_timer_database();
  qassert(tdb_history.size() >= 1);
  qassert(tdb.size() >= tdb_history.back().size());
  for (long i = 0; i < (long)tdb_history.back().size(); ++i) {
    tdb[i].merge(tdb_history.back()[i]);
  }
  tdb_history.pop_back();
  displayln_info(0, ssprintf("Timer::merge(): Merge all timers! (level = %ld)",
                             (long)tdb_history.size()));
}

void Timer::start(bool verbose)
{
  get_timer_stack().push_back(info_index);
  if (is_running > 0) {
    is_running += 1;
    return;
  } else if (is_running == 0) {
    is_running = 1;
  } else {
    TimerInfo& info = get_timer_database()[info_index];
    info.show_avg_always("debug", max_function_name_length_shown());
    displayln(ssprintf("%s::%s ERROR: is_running=%d", cname, info.fname.c_str(),
                       is_running));
    Timer::display_stack();
    qassert(false);
  }
  TimerInfo& info = get_timer_database()[info_index];
  info.call_times++;
  if (get_verbose_level() > 0) {
    if (verbose ||
        info.accumulated_time >=
            info.call_times * minimum_duration_for_show_info() ||
        info.call_times <= max_call_times_for_always_show_info()) {
      info.show_start(max_function_name_length_shown());
    }
  }
  start_flops = is_using_total_flops ? get_total_flops() : 0;
  flops = 0;
  start_time = get_time();
}

void Timer::stop(bool verbose)
{
  TimerInfo& info = get_timer_database()[info_index];
  std::vector<long>& t_stack = get_timer_stack();
  qassert(not t_stack.empty());
  if (not(t_stack.back() == info_index)) {
    displayln(ssprintf("%s::%s ERROR: stack is corrupted", cname,
                       info.fname.c_str()));
    Timer::display_stack();
    qassert(false);
  }
  t_stack.pop_back();
  if (is_running <= 0) {
    info.show_avg_always("debug", max_function_name_length_shown());
    displayln(ssprintf("%s::%s ERROR: is_running=%d", cname, info.fname.c_str(),
                       is_running));
    Timer::display_stack();
    qassert(false);
  }
  is_running -= 1;
  if (is_running != 0) {
    return;
  }
  stop_time = get_time();
  if (is_using_total_flops) {
    stop_flops = get_total_flops();
  } else {
    stop_flops = start_flops + flops;
  }
  info.dtime = stop_time - start_time;
  info.dflops = stop_flops - start_flops;
  bool is_show = false;
  if (get_verbose_level() > 0) {
    if (verbose ||
        info.accumulated_time >=
            info.call_times * minimum_duration_for_show_info() ||
        info.dtime >= 2.0 * minimum_duration_for_show_info() ||
        info.call_times <= max_call_times_for_always_show_info()) {
      is_show = true;
    }
  }
  info.accumulated_time += info.dtime;
  info.accumulated_flops += info.dflops;
  if (is_show) {
    info.show_stop(max_function_name_length_shown());
  }
}

void Timer::test_timer_time_usage()
{
  {
    static Timer timer("Timer", false);
    static Timer timer_test("Timer-test");
    timer.start();
    timer_test.start();
    timer_test.stop();
    timer.stop();
  }
  {
    static Timer timer_noflop("Timer-noflop", false);
    static Timer timer_test("Timer-test", false);
    timer_noflop.start();
    timer_test.start();
    timer_test.stop();
    timer_noflop.stop();
  }
}

void Timer::display(const std::string& tag)
{
  double total_time = get_total_time();
  const std::vector<TimerInfo>& tdb = get_timer_database();
  const std::vector<std::vector<TimerInfo>>& tdb_history =
      get_timer_database_history();
  std::vector<const TimerInfo*> db;
  const long tdbsize = tdb.size();
  for (long i = 0; i < tdbsize; i++) {
    db.push_back(&tdb[i]);
  }
  std::sort(db.begin(), db.end(), compare_time_info_p);
  displayln_info(ssprintf(
      "Timer::display-start: %s (level=%ld) fname : time%% number of calls; "
      "Avg,Tot sec; Avg,Tot flops; Gflops",
      tag.c_str(), (long)tdb_history.size()));
  const long dbsize = db.size();
  for (long i = 0; i < dbsize; i++) {
    if (db[i]->call_times > 0) {
      db[i]->show_avg("display", max_function_name_length_shown());
    }
  }
  displayln_info(ssprintf(
      "Timer::display-end:   %s (level=%ld) --------------------- total %.4E "
      "sec ----------------------",
      tag.c_str(), (long)tdb_history.size(), total_time));
}

void Timer::autodisplay(const double time)
{
  static double last_time = get_start_time();
  if (time - last_time > minimum_autodisplay_interval()) {
    last_time = time;
    display("autodisplay");
  }
}

void Timer::display_stack_always()
// display on any process or thread
{
  displayln("display_stack start");
  const std::vector<TimerInfo>& tdb = get_timer_database();
  const std::vector<long>& t_stack = get_timer_stack();
  for (long i = (long)t_stack.size() - 1; i >= 0; --i) {
    const long info_index = t_stack[i];
    tdb[info_index].show_avg_always(ssprintf("stack[%3ld]", i),
                                    max_function_name_length_shown());
  }
  displayln("display_stack end");
}

}  // namespace qlat
